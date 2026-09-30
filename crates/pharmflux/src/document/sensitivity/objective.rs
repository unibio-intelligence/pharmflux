//! Gaussian negative log likelihood with variance derivatives retained.
use super::*;
use pharmflux_core::{fit::*, identity::canonical_hash};
fn normal(
    prediction: f64,
    value: f64,
    additive: f64,
    proportional: f64,
) -> Result<(f64, f64), Error> {
    if ![prediction, value, additive, proportional]
        .iter()
        .all(|x| x.is_finite())
        || additive < 0.
        || proportional < 0.
        || additive == 0. && proportional == 0.
    {
        return Err(invalid(
            "invalid Gaussian prediction, observation or standard deviation",
        ));
    }
    let component = proportional * prediction;
    let sd = libm::hypot(additive, component);
    if !sd.is_finite() || sd <= 0. {
        return Err(Error::new(
            ErrorCode::Domain,
            "Gaussian variance must be finite and positive",
        ));
    }
    let difference = value - prediction;
    let z = if difference.is_finite() {
        difference / sd
    } else {
        value / sd - prediction / sd
    };
    let nll = 0.5 * z * z + libm::log(sd) + 0.9189385332046727;
    let dlogsd = (component / sd) * (proportional / sd);
    let derivative = -z / sd + (1. - z * z) * dlogsd;
    if !nll.is_finite() || !derivative.is_finite() {
        return Err(Error::new(
            ErrorCode::Domain,
            "Gaussian likelihood or derivative overflow",
        ));
    }
    Ok((nll, derivative))
}
impl CompiledSensitivityDocument {
    /// Likelihood-only path for estimators that do not need output derivatives.
    /// Keep the same observation, error and output-budget contract as the
    /// derivative-bearing Gaussian objective.
    pub(super) fn gaussian_primal_nll(
        self: &Arc<Self>,
        request: &GaussianObjectiveRequest,
    ) -> Result<f64, Error> {
        if !request.priors.is_empty() {
            return Err(invalid("population subject priors are unsupported"));
        }
        if request.simulation.with_respect_to != self.parameters()
            || !matches!(
                request.simulation.run.solver,
                pharmflux_core::run::Solver::DiffsolBdf
            )
            || request.simulation.run.seed.is_some()
        {
            return Err(Error::new(
                ErrorCode::Unsupported,
                "population likelihood requires deterministic BDF and the selected parameter order",
            ));
        }
        let observations = &request.observations;
        if observations.is_empty()
            || observations.len() > 1_000_000
            || observations.len() > request.simulation.run.budgets.output_values
        {
            return Err(invalid(
                "Gaussian objective requires a bounded nonempty observation table",
            ));
        }
        let storage = observations
            .len()
            .checked_mul(5)
            .and_then(|n| n.checked_add(2 * self.parameters().len()))
            .ok_or_else(|| invalid("objective storage overflow"))?;
        let mut run = request.simulation.run.clone();
        run.budgets.output_values =
            run.budgets
                .output_values
                .checked_sub(storage)
                .ok_or_else(|| {
                    Error::new(
                        ErrorCode::OutputBudget,
                        "objective observation storage exceeds budget",
                    )
                })?;
        let mut data = Vec::with_capacity(observations.len());
        for (i, observation) in observations.iter().enumerate() {
            let lower = || -> Result<_, Error> {
                let index = self
                    .output_names()
                    .iter()
                    .position(|name| name == &observation.output)
                    .ok_or_else(|| invalid("unknown observation output"))?;
                let unit = Unit::parse(&self.output_units()[index])?;
                let value = quantity(&observation.value)?.in_unit(unit)?;
                let additive = quantity(&observation.error.additive_sd)?.in_unit(unit)?;
                let proportional = observation.error.proportional_sd;
                if !proportional.is_finite()
                    || proportional < 0.
                    || additive < 0.
                    || additive == 0. && proportional == 0.
                {
                    return Err(invalid(
                        "Gaussian standard deviations must be nonnegative and not both zero",
                    ));
                }
                Ok((index, observation.row, value, additive, proportional))
            };
            data.push(lower().map_err(|mut error| {
                error.expression = Some(format!("observations[{i}]"));
                error
            })?);
        }
        let result = self.primal.execute(&run)?;
        let mut nll = 0.0;
        for (i, (output, row, value, additive, proportional)) in data.into_iter().enumerate() {
            let prediction = *result.outputs[output]
                .values
                .get(row)
                .ok_or_else(|| invalid("observation row is outside the simulation result"))?;
            let (term, _) =
                normal(prediction, value, additive, proportional).map_err(|mut error| {
                    error.expression = Some(format!("observations[{i}]"));
                    error
                })?;
            nll += term;
        }
        if !nll.is_finite() {
            return Err(Error::new(
                ErrorCode::Domain,
                "Gaussian objective accumulation overflow",
            ));
        }
        Ok(nll)
    }

    /// Independent Gaussian observations with known error coefficients and
    /// optional independent normal priors on selected natural parameters.
    pub fn gaussian_objective(
        self: &Arc<Self>,
        request: &GaussianObjectiveRequest,
    ) -> Result<GaussianObjectiveResult, Error> {
        if request.observations.is_empty()
            || request.observations.len() > 1_000_000
            || request.observations.len() > request.simulation.run.budgets.output_values
        {
            return Err(invalid(
                "Gaussian objective requires a bounded nonempty observation table",
            ));
        }
        if request.priors.len() > self.parameters().len() {
            return Err(invalid("too many parameter priors"));
        }
        let storage = request
            .observations
            .len()
            .checked_mul(5)
            .and_then(|n| n.checked_add(2 * self.parameters().len() + 3 * request.priors.len()))
            .ok_or_else(|| invalid("objective storage overflow"))?;
        let mut simulation = request.simulation.clone();
        simulation.run.budgets.output_values = simulation
            .run
            .budgets
            .output_values
            .checked_sub(storage)
            .ok_or_else(|| {
                Error::new(
                    ErrorCode::OutputBudget,
                    "objective observation storage exceeds budget",
                )
            })?;
        let mut data = Vec::with_capacity(request.observations.len());
        for (i, observation) in request.observations.iter().enumerate() {
            let lower = || -> Result<_, Error> {
                let index = self
                    .output_names()
                    .iter()
                    .position(|n| n == &observation.output)
                    .ok_or_else(|| invalid("unknown observation output"))?;
                let unit = Unit::parse(&self.output_units()[index])?;
                let value = quantity(&observation.value)?.in_unit(unit)?;
                let additive = quantity(&observation.error.additive_sd)?.in_unit(unit)?;
                let proportional = observation.error.proportional_sd;
                if !proportional.is_finite()
                    || proportional < 0.
                    || additive < 0.
                    || additive == 0. && proportional == 0.
                {
                    return Err(invalid(
                        "Gaussian standard deviations must be nonnegative and not both zero",
                    ));
                }
                Ok((index, observation.row, value, additive, proportional))
            };
            data.push(lower().map_err(|mut e| {
                e.expression = Some(format!("observations[{i}]"));
                e
            })?);
        }
        let bound = self.bind_with_initial_states(
            &request.simulation.run.parameters,
            &request.simulation.run.regimen.covariates,
            &request.simulation.run.initial_states,
        )?;
        let mut gradient = vec![0.; self.parameters().len()];
        let mut prior_nll = 0.;
        let mut names = BTreeSet::new();
        for prior in &request.priors {
            let p = self
                .parameters()
                .iter()
                .position(|n| n == &prior.parameter)
                .ok_or_else(|| invalid("prior must name a selected parameter"))?;
            if !names.insert(&prior.parameter) {
                return Err(invalid("duplicate parameter prior"));
            }
            let index = self
                .primal
                .parameters
                .iter()
                .position(|p| p.name == prior.parameter)
                .ok_or_else(|| invalid("unknown prior parameter"))?;
            let unit = Unit::parse(&self.primal.parameters[index].default.unit)?;
            let mean = quantity(&prior.mean)?.in_unit(unit)?;
            let sd = quantity(&prior.sd)?.in_unit(unit)?;
            if sd <= 0. {
                return Err(invalid("prior standard deviation must be positive"));
            }
            let (nll, derivative) = normal(bound.primal.values[index], mean, sd, 0.)?;
            prior_nll += nll;
            gradient[p] += derivative;
        }
        let result = self.run(&simulation)?;
        let mut nll = 0.;
        let width = self.output_names().len();
        for (i, (output, row, value, additive, proportional)) in data.into_iter().enumerate() {
            let prediction = *result.result.outputs[output]
                .values
                .get(row)
                .ok_or_else(|| invalid("observation row is outside the simulation result"))?;
            let (term, derivative) =
                normal(prediction, value, additive, proportional).map_err(|mut e| {
                    e.expression = Some(format!("observations[{i}]"));
                    e
                })?;
            nll += term;
            for (p, g) in gradient.iter_mut().enumerate() {
                *g += derivative * result.derivatives[p * width + output].values[row];
            }
        }
        let objective = nll + prior_nll;
        if !objective.is_finite() || gradient.iter().any(|g| !g.is_finite()) {
            return Err(Error::new(
                ErrorCode::Domain,
                "Gaussian objective accumulation overflow",
            ));
        }
        let mut identity = result.result.identity;
        identity.run_request_hash = canonical_hash(
            &serde_json::json!({"simulation_request_hash":identity.run_request_hash,"requested_budgets":request.simulation.run.budgets,"observations":request.observations,"priors":request.priors}),
        )?;
        identity.algorithm = "bdf_forward_gaussian_objective".into();
        let gradient = self
            .parameters()
            .iter()
            .zip(gradient)
            .map(|(parameter, value)| {
                let unit = &self
                    .primal
                    .parameters
                    .iter()
                    .find(|p| &p.name == parameter)
                    .expect("selected parameter validated")
                    .default
                    .unit;
                ObjectiveGradient {
                    parameter: parameter.clone(),
                    unit: format!("1/({unit})"),
                    value,
                }
            })
            .collect();
        Ok(GaussianObjectiveResult {
            identity,
            negative_log_likelihood: nll,
            negative_log_prior: prior_nll,
            objective,
            gradient,
            observations: request.observations.len(),
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn primal_likelihood_matches_derivative_objective() {
        let model: pharmflux_core::model::ModelDocument = serde_json::from_str(include_str!(
            "../../../../../conformance/models/synthetic-one-compartment.json"
        ))
        .unwrap();
        let compiled = CompiledSensitivityDocument::compile(&model, &["cl".into()]).unwrap();
        let mut simulation: pharmflux_core::run::SensitivityRequest = serde_json::from_str(
            include_str!("../../../../../conformance/requests/synthetic-sensitivity.json"),
        )
        .unwrap();
        simulation.with_respect_to = vec!["cl".into()];
        simulation.derivative_absolute_tolerances.remove("v");
        let observations = [(2, 1.01), (3, 0.93), (4, 0.79)]
            .into_iter()
            .map(|(row, value)| GaussianObservation {
                row,
                output: "cp".into(),
                value: pharmflux_core::model::Quantity {
                    value,
                    unit: "mg/L".into(),
                },
                error: GaussianError {
                    additive_sd: pharmflux_core::model::Quantity {
                        value: 0.03,
                        unit: "mg/L".into(),
                    },
                    proportional_sd: 0.04,
                },
            })
            .collect();
        let request = GaussianObjectiveRequest {
            simulation,
            observations,
            priors: vec![],
        };
        let expected = compiled.gaussian_objective(&request).unwrap();
        let actual = compiled.gaussian_primal_nll(&request).unwrap();
        assert!(
            (actual - expected.negative_log_likelihood).abs() < 1e-6,
            "primal={actual}, sensitivity={}",
            expected.negative_log_likelihood
        );
    }
}
