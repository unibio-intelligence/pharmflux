use super::*;
use pharmflux_core::{
    identity::{canonical_hash, SPEC_VERSION},
    run::*,
};
impl CompiledDocument {
    /// Solve a fully contained periodic linear PK cycle, returning pre-event amounts.
    pub fn periodic_steady_state(
        self: &Arc<Self>,
        request: &SteadyStateRequest,
    ) -> Result<SteadyStateResult, Error> {
        let run = &request.run;
        if !matches!(run.solver, Solver::AnalyticLinearPk) || self.linear_pk.is_none() {
            return Err(Error::new(
                ErrorCode::Unsupported,
                "periodic steady state requires analytic_linear_pk",
            ));
        }
        if run.seed.is_some() || !run.initial_states.is_empty() {
            return Err(invalid(
                "steady state prohibits seeds and initial state overrides",
            ));
        }
        if run
            .request_id
            .as_ref()
            .is_some_and(|s| s.is_empty() || s.len() > 128)
        {
            return Err(invalid("request_id must contain 1..128 bytes"));
        }
        if run.budgets.solver_callbacks > u32::MAX as u64
            || run.budgets.events > u32::MAX as usize
            || run.budgets.output_values > u32::MAX as usize
        {
            return Err(invalid(
                "budgets must fit the shared native/WASM unsigned 32-bit range",
            ));
        }
        if !self.history.is_empty()
            || !self.switches.is_empty()
            || !self.dose_switches.is_empty()
            || !run.regimen.covariate_changes.is_empty()
            || !run.regimen.checkpoints.is_empty()
        {
            return Err(Error::new(ErrorCode::Unsupported, "steady state requires a constant-covariate cycle without history, switches or checkpoints"));
        }
        if !run.regimen.samples.is_empty() || run.regimen.observations.is_some() {
            return Err(invalid(
                "steady-state requests do not take trajectory observations",
            ));
        }
        let bound = self.bind_with_covariates(&run.parameters, &run.regimen.covariates)?;
        let kernel = self.linear_pk.as_ref().unwrap().bind(&bound.values)?;
        let lagged =
            self.expand_dose_lags(&run.regimen, &run.parameters, &run.regimen.covariates)?;
        let mut protocol = self.lower_regimen_protocol_internal(&lagged, run.budgets, true)?;
        // The low-level kernel receives state amounts. Preserve the document's
        // input-unit conversion and model dose scaling exactly once.
        for event in &mut protocol.events {
            match event {
                crate::Event::Bolus { target, amount, .. }
                | crate::Event::Infusion { target, amount, .. } => {
                    let scale = bound
                        .dose_scale(*target)
                        .ok_or_else(|| invalid("missing dose scale"))?;
                    *amount *= scale;
                }
                crate::Event::Reset { .. } => {}
            }
        }
        let absolute = Unit::parse(&request.absolute_tolerance.unit)?
            .convert(request.absolute_tolerance.value, self.model.state_units[0])?;
        let result = kernel.periodic_steady_state(
            &protocol,
            run.budgets,
            crate::linear_pk::SteadyStateTolerance {
                absolute,
                relative: request.relative_tolerance,
            },
        )?;
        let make_map = |values: &[f64]| {
            self.state_names
                .iter()
                .zip(values)
                .map(|(name, value)| {
                    (
                        name.clone(),
                        Quantity {
                            value: *value,
                            unit: self.linear_pk.as_ref().unwrap().amount_unit.clone(),
                        },
                    )
                })
                .collect()
        };
        let identity = ExecutionIdentity {
            model_content_hash: self.model_content_hash.clone(),
            spec_version: SPEC_VERSION.into(),
            compiler_version: concat!("pharmflux/", env!("CARGO_PKG_VERSION")).into(),
            compiler_source_sha256: env!("PHARMFLUX_SOURCE_SHA256").into(),
            rustc_version: env!("PHARMFLUX_RUSTC_VERSION").into(),
            target: env!("PHARMFLUX_TARGET").into(),
            build_profile: env!("PHARMFLUX_PROFILE").into(),
            rustflags_sha256: env!("PHARMFLUX_RUSTFLAGS_SHA256").into(),
            bound_value_hash: canonical_hash(&serde_json::json!(bound.values))?,
            run_request_hash: canonical_hash(
                &serde_json::json!({"operation":"periodic_steady_state",
                "schema":run.schema,"regimen":run.regimen,"budgets":run.budgets,
                "absolute_tolerance":request.absolute_tolerance,"relative_tolerance":request.relative_tolerance}),
            )?,
            backend: "pharmflux".into(),
            backend_version: env!("CARGO_PKG_VERSION").into(),
            algorithm: "linear_pk_periodic_affine_taylor20".into(),
        };
        Ok(SteadyStateResult {
            identity,
            initial_states: make_map(&result.amounts),
            residuals: make_map(&result.residual),
            maximum_scaled_residual: result.maximum_scaled_residual,
        })
    }
}

impl CompiledDocument {
    /// Iterate a repeated phase-time protocol using BDF, with a unit-bearing seed.
    pub fn iterate_periodic_state(
        self: &Arc<Self>,
        request: &CycleIterationRequest,
    ) -> Result<CycleIterationResult, Error> {
        let run = &request.run;
        if !matches!(run.solver, Solver::DiffsolBdf) {
            return Err(Error::new(
                ErrorCode::Unsupported,
                "cycle iteration requires diffsol_bdf",
            ));
        }
        if run.seed.is_some() {
            return Err(invalid("stochastic seeds are not supported"));
        }
        if run
            .request_id
            .as_ref()
            .is_some_and(|s| s.is_empty() || s.len() > 128)
        {
            return Err(invalid("request_id must contain 1..128 bytes"));
        }
        if run.budgets.solver_callbacks > u32::MAX as u64
            || run.budgets.events > u32::MAX as usize
            || run.budgets.output_values > u32::MAX as usize
        {
            return Err(invalid(
                "budgets must fit the shared native/WASM unsigned 32-bit range",
            ));
        }
        if !self.history.is_empty()
            || !self.switches.is_empty()
            || !self.dose_switches.is_empty()
            || !run.regimen.covariate_changes.is_empty()
            || !run.regimen.checkpoints.is_empty()
        {
            return Err(Error::new(
                ErrorCode::Unsupported,
                "cycle iteration requires fixed covariates and no history, switches or checkpoints",
            ));
        }
        if !run.regimen.samples.is_empty() || run.regimen.observations.is_some() {
            return Err(invalid(
                "cycle iteration does not take trajectory observations",
            ));
        }
        if request.absolute_tolerances.len() != self.state_count {
            return Err(invalid("provide one convergence tolerance per state"));
        }
        let absolute = self
            .state_names
            .iter()
            .enumerate()
            .map(|(i, name)| {
                let q = request
                    .absolute_tolerances
                    .get(name)
                    .ok_or_else(|| invalid(format!("missing convergence tolerance: {name}")))?;
                Unit::parse(&q.unit)?.convert(q.value, self.model.state_units[i])
            })
            .collect::<Result<Vec<_>, Error>>()?;
        let bound = self.bind_initial_states_internal(
            &run.parameters,
            &run.regimen.covariates,
            &run.initial_states,
        )?;
        let initial = bound.initial();
        let lagged =
            self.expand_dose_lags(&run.regimen, &run.parameters, &run.regimen.covariates)?;
        let protocol = self.lower_regimen_protocol_internal(&lagged, run.budgets, true)?;
        let options = crate::steady_state::CycleIterationOptions {
            absolute_tolerances: absolute,
            integration_absolute_tolerances: self.absolute_tolerances(&run.regimen)?,
            relative_tolerance: request.relative_tolerance,
            maximum_cycles: request.maximum_cycles,
        };
        let result = crate::steady_state::iterate_periodic_cycle(
            |state| {
                let mut next =
                    self.bind_continuing(&run.parameters, &run.regimen.covariates, Some(state))?;
                next.scheduled = true;
                Ok(next)
            },
            &initial,
            &protocol,
            &options,
            run.budgets,
        )?;
        let make_map = |values: &[f64]| {
            self.state_names
                .iter()
                .zip(&self.state_unit_names)
                .zip(values)
                .map(|((name, unit), value)| {
                    (
                        name.clone(),
                        Quantity {
                            value: *value,
                            unit: unit.clone(),
                        },
                    )
                })
                .collect()
        };
        let identity = ExecutionIdentity {
            model_content_hash: self.model_content_hash.clone(),
            spec_version: SPEC_VERSION.into(),
            compiler_version: concat!("pharmflux/", env!("CARGO_PKG_VERSION")).into(),
            compiler_source_sha256: env!("PHARMFLUX_SOURCE_SHA256").into(),
            rustc_version: env!("PHARMFLUX_RUSTC_VERSION").into(),
            target: env!("PHARMFLUX_TARGET").into(),
            build_profile: env!("PHARMFLUX_PROFILE").into(),
            rustflags_sha256: env!("PHARMFLUX_RUSTFLAGS_SHA256").into(),
            bound_value_hash: canonical_hash(
                &serde_json::json!({"values":bound.values,"seed":initial}),
            )?,
            run_request_hash: canonical_hash(
                &serde_json::json!({"operation":"iterate_periodic_state","schema":run.schema,
                "regimen":run.regimen,"budgets":run.budgets,"time_basis":request.time_basis,
                "absolute_tolerances":request.absolute_tolerances,"relative_tolerance":request.relative_tolerance,"maximum_cycles":request.maximum_cycles}),
            )?,
            backend: "diffsol".into(),
            backend_version: crate::DIFFSOL_BACKEND_VERSION.into(),
            algorithm: "bdf_periodic_cycle_iteration".into(),
        };
        Ok(CycleIterationResult {
            steady_state: SteadyStateResult {
                identity,
                initial_states: make_map(&result.initial_states),
                residuals: make_map(&result.residuals),
                maximum_scaled_residual: result.maximum_scaled_residual,
            },
            cycles_completed: result.cycles_completed,
        })
    }
}
