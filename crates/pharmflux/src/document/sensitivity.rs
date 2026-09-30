//! Typed document compilation for smooth, constant-binding sensitivities.
mod fit;
mod objective;
mod outputs;
mod population_focei;
mod population_laplace;
mod population_saem;
mod regimen;
mod run;
use super::*;
use crate::runtime::sensitivity::*;
pub use outputs::SensitivitySample;
#[derive(Debug)]
pub struct CompiledSensitivityDocument {
    primal: Arc<CompiledDocument>,
    forward: Arc<UnitCompiledForwardSensitivity>,
    readouts: Arc<UnitCompiledSensitivityReadouts>,
    binding_units: Vec<Unit>,
}
#[derive(Debug)]
pub struct BoundSensitivityDocument {
    compiled: Arc<CompiledSensitivityDocument>,
    primal: BoundDocument,
    forward: UnitBoundForwardSensitivity,
    readouts: UnitBoundSensitivityReadouts,
    checkpoints: Vec<f64>,
}
impl CompiledSensitivityDocument {
    pub fn compile(document: &ModelDocument, parameters: &[String]) -> Result<Arc<Self>, Error> {
        let primal = CompiledDocument::compile(document)?;
        let lowered;
        let document = if document.linear_pk.is_some() {
            lowered = CompiledDocument::lower_linear_pk(document)?;
            &lowered
        } else {
            document
        };
        if !document.switches.is_empty()
            || !document.dose_switches.is_empty()
            || !document.dose_history.is_empty()
        {
            return Err(Error::new(
                ErrorCode::Unsupported,
                "document sensitivities do not yet support history or moving event boundaries",
            ));
        }
        if parameters
            .iter()
            .any(|name| !document.parameters.iter().any(|p| &p.name == name))
        {
            return Err(invalid(
                "select declared model parameters for document sensitivities",
            ));
        }
        let mut definition = runtime::UnitModelDefinition {
            time_unit: Unit::parse(&document.time_unit)?,
            bindings: document
                .parameters
                .iter()
                .map(|p| {
                    Ok(runtime::BindingDefinition {
                        name: p.name.clone(),
                        unit: Unit::parse(&p.default.unit)?,
                        bounds: p
                            .bounds
                            .as_ref()
                            .map(|[a, b]| Ok::<_, Error>((quantity(a)?, quantity(b)?)))
                            .transpose()?,
                    })
                })
                .collect::<Result<_, Error>>()?,
            states: document
                .states
                .iter()
                .map(|s| {
                    Ok(runtime::UnitStateDefinition {
                        name: s.name.clone(),
                        unit: Unit::parse(&s.unit)?,
                        initial: match &s.initial {
                            Initial::Quantity { quantity: q } => {
                                runtime::InitialCondition::Quantity(quantity(q)?)
                            }
                            Initial::Expression { expression } => {
                                runtime::InitialCondition::Expression(expression.clone())
                            }
                        },
                        rhs: s.rhs.clone(),
                        dosing: s
                            .dosing
                            .as_ref()
                            .map(|d| {
                                Ok::<_, Error>((Unit::parse(&d.amount_unit)?, d.scale.clone()))
                            })
                            .transpose()?,
                    })
                })
                .collect::<Result<_, Error>>()?,
        };
        for c in &document.covariates {
            definition.bindings.push(runtime::BindingDefinition {
                name: c.name.clone(),
                unit: Unit::parse(&c.unit)?,
                bounds: None,
            });
        }
        let derived = document
            .individual
            .iter()
            .map(|p| {
                Ok(UnitDerivedBinding {
                    name: p.name.clone(),
                    unit: primal
                        .individual
                        .iter()
                        .find(|(n, _, _)| n == &p.name)
                        .ok_or_else(|| invalid("missing individual unit"))?
                        .1,
                    expression: p.expression.clone(),
                })
            })
            .collect::<Result<Vec<_>, Error>>()?;
        let mut outputs = document
            .outputs
            .iter()
            .map(|o| {
                Ok(UnitSensitivityReadout {
                    expression: o.expression.clone(),
                    unit: Unit::parse(&o.unit)?,
                })
            })
            .collect::<Result<Vec<_>, Error>>()?;
        let mut lags = document
            .states
            .iter()
            .filter_map(|s| {
                s.dosing
                    .as_ref()
                    .and_then(|d| d.lag.as_ref())
                    .map(|e| (s.name.clone(), e.clone()))
            })
            .collect::<Vec<_>>();
        if !document.definitions.is_empty() {
            let mut symbols = definition
                .bindings
                .iter()
                .map(|b| {
                    (
                        Symbol {
                            name: b.name.clone(),
                            kind: SymbolKind::Parameter,
                        },
                        b.unit,
                    )
                })
                .collect::<Vec<_>>();
            symbols.extend(derived.iter().map(|b| {
                (
                    Symbol {
                        name: b.name.clone(),
                        kind: SymbolKind::Parameter,
                    },
                    b.unit,
                )
            }));
            symbols.extend(definition.states.iter().map(|s| {
                (
                    Symbol {
                        name: s.name.clone(),
                        kind: SymbolKind::State,
                    },
                    s.unit,
                )
            }));
            symbols.push((
                Symbol {
                    name: "time".into(),
                    kind: SymbolKind::Time,
                },
                definition.time_unit,
            ));
            for s in &definition.states {
                symbols.push((
                    Symbol {
                        name: format!("input_{}", s.name),
                        kind: SymbolKind::Input,
                    },
                    s.unit.divide(definition.time_unit)?,
                ));
            }
            let mut expanded = definitions::Expanded::new(&document.definitions, &symbols)?;
            for state in &mut definition.states {
                state.rhs = expanded.expand(&state.rhs, "sensitivity.rhs")?;
                if let runtime::InitialCondition::Expression(e) = &mut state.initial {
                    *e = expanded.expand(e, "sensitivity.initial")?;
                }
                if let Some((_, e)) = &mut state.dosing {
                    *e = expanded.expand(e, "sensitivity.dosing")?;
                }
            }
            for o in &mut outputs {
                o.expression = expanded.expand(&o.expression, "sensitivity.output")?;
            }
            for (name, expression) in &mut lags {
                *expression = expanded.expand(expression, &format!("states.{name}.dosing.lag"))?;
            }
        }
        fn depends(expression: &pharmflux_core::expression::Expr, parameter: &str) -> bool {
            matches!(expression,pharmflux_core::expression::Expr::Symbol{name} if name==parameter)
                || expression.children().iter().any(|e| depends(e, parameter))
        }
        for (state, expression) in lags {
            if let Some(parameter) = parameters.iter().find(|p| depends(&expression, p)) {
                let mut error=Error::new(ErrorCode::Unsupported,format!("model lag depends on selected sensitivity parameter {parameter}; moving-event derivatives are not implemented"));
                error.expression = Some(format!("states.{state}.dosing.lag"));
                return Err(error);
            }
        }
        let (definition, outputs) = resolve_unit_bindings(&definition, &derived, &outputs)?;
        Ok(Arc::new(Self {
            forward: UnitCompiledForwardSensitivity::compile(&definition, parameters)?,
            readouts: UnitCompiledSensitivityReadouts::compile(&definition, parameters, &outputs)?,
            binding_units: definition.bindings.iter().map(|b| b.unit).collect(),
            primal,
        }))
    }
    pub fn parameters(&self) -> &[String] {
        self.forward.parameters()
    }
    pub fn output_names(&self) -> &[String] {
        &self.primal.output_names
    }
    pub fn output_units(&self) -> &[String] {
        &self.primal.output_units
    }
    pub fn derivative_units(&self) -> &[Vec<Unit>] {
        &self.readouts.derivative_units
    }
    pub fn augmented_state_units(&self) -> &[Unit] {
        &self.forward.augmented_state_units
    }
    pub fn bind(
        self: &Arc<Self>,
        parameters: &BTreeMap<String, Quantity>,
        covariates: &BTreeMap<String, Quantity>,
    ) -> Result<BoundSensitivityDocument, Error> {
        self.bind_with_initial_states(parameters, covariates, &BTreeMap::new())
    }
    pub fn bind_with_initial_states(
        self: &Arc<Self>,
        parameters: &BTreeMap<String, Quantity>,
        covariates: &BTreeMap<String, Quantity>,
        initial_states: &BTreeMap<String, Quantity>,
    ) -> Result<BoundSensitivityDocument, Error> {
        let primal =
            self.primal
                .bind_with_initial_states(parameters, covariates, initial_states)?;
        let overrides = self
            .primal
            .state_names
            .iter()
            .zip(&self.primal.model.state_units)
            .map(|(name, unit)| {
                initial_states
                    .get(name)
                    .map(|q| quantity(q)?.in_unit(*unit))
                    .transpose()
            })
            .collect::<Result<Vec<_>, Error>>()?;
        let values = primal
            .values
            .iter()
            .zip(&self.binding_units)
            .map(|(v, u)| runtime::Quantity {
                value: *v,
                unit: *u,
            })
            .collect::<Vec<_>>();
        Ok(BoundSensitivityDocument {
            compiled: self.clone(),
            checkpoints: Vec::new(),
            forward: self.forward.bind_initial(&values, Some(&overrides))?,
            readouts: self.readouts.bind_initial(&values, Some(&overrides))?,
            primal,
        })
    }
}
impl BoundSensitivityDocument {
    pub fn rebind(
        &mut self,
        parameters: &BTreeMap<String, Quantity>,
        covariates: &BTreeMap<String, Quantity>,
    ) -> Result<(), Error> {
        *self = self.compiled.bind_with_initial_states(
            parameters,
            covariates,
            &self.primal.initial_overrides,
        )?;
        Ok(())
    }
    pub fn outputs(
        &self,
        time: f64,
        augmented: &[f64],
        raw_rates: &[f64],
    ) -> Result<SensitivityReadoutValues, Error> {
        self.readouts.evaluate(time, augmented, raw_rates)
    }
}
impl Model for BoundSensitivityDocument {
    fn discontinuities(&self) -> &[f64] {
        &self.checkpoints
    }
    fn validate(&self) -> Result<(), Error> {
        self.primal.validate()?;
        self.forward.validate()
    }
    fn initial(&self) -> Vec<f64> {
        self.forward.initial()
    }
    fn check_delivery(&self, event: &pharmflux_core::Event) -> Result<(), Error> {
        self.primal.check_delivery(event)
    }
    fn check_state(&self, time: f64, state: &[f64]) -> Result<(), Error> {
        if state.len() != self.compiled.forward.augmented_state_units.len() {
            return Err(invalid("invalid augmented state dimension"));
        }
        self.primal
            .check_state(time, &state[..self.compiled.primal.state_count])
    }
    fn dose_scale(&self, target: usize) -> Option<f64> {
        self.forward.dose_scale(target)
    }
    fn apply_bolus(
        &self,
        time: f64,
        target: usize,
        amount: f64,
        state: &mut [f64],
    ) -> Result<(), Error> {
        self.forward.apply_bolus(time, target, amount, state)
    }
    fn apply_reset(
        &self,
        time: f64,
        target: usize,
        value: f64,
        state: &mut [f64],
    ) -> Result<(), Error> {
        self.forward.apply_reset(time, target, value, state)
    }
    fn accumulate_infusion_rate(
        &self,
        time: f64,
        target: usize,
        rate: f64,
        rates: &mut [f64],
    ) -> Result<(), Error> {
        self.forward
            .accumulate_infusion_rate(time, target, rate, rates)
    }
    fn rhs(&self, time: f64, state: &[f64], rates: &[f64], out: &mut [f64]) -> Result<(), Error> {
        self.forward.rhs(time, state, rates, out)
    }
    fn jac_mul(
        &self,
        time: f64,
        state: &[f64],
        rates: &[f64],
        v: &[f64],
        out: &mut [f64],
    ) -> Result<(), Error> {
        self.forward.jac_mul(time, state, rates, v, out)
    }
}
