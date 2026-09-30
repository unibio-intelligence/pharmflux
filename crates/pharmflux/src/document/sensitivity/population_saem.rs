//! Bounded SAEM for one or two lognormal population parameters.
//!
//! Stochastic approximation tracks the first and second moments of the log
//! individual parameters. A final observed-curvature Laplace calculation
//! reports an approximate marginal OFV when the work budget permits it.
use super::*;
use pharmflux_core::{fit::*, identity::canonical_hash};

const LOG_2PI: f64 = 1.8378770664093453;
type Vector = [f64; 2];
type Matrix = [[f64; 2]; 2];

struct Setup {
    dimensions: usize,
    names: [String; 2],
    units: [String; 2],
    start: Vector,
    lower: Vector,
    upper: Vector,
    omega: Matrix,
    omega_lower: Vector,
    omega_upper: Vector,
}

struct ResidualFit {
    output: String,
    unit: String,
    variance: f64,
    lower_variance: f64,
    upper_variance: f64,
    observations: usize,
}

fn validate_residual_fit(
    model: &CompiledSensitivityDocument,
    request: &PopulationFitRequest,
) -> Result<Option<ResidualFit>, Error> {
    let Some(fit) = &request.additive_error_fit else {
        return Ok(None);
    };
    let output_index = model
        .output_names()
        .iter()
        .position(|name| name == &fit.output)
        .ok_or_else(|| invalid("unknown SAEM additive-error fit output"))?;
    let unit = model.output_units()[output_index].clone();
    let parsed_unit = Unit::parse(&unit)?;
    let initial = quantity(&fit.initial)?.in_unit(parsed_unit)?;
    let lower = quantity(&fit.lower)?.in_unit(parsed_unit)?;
    let upper = quantity(&fit.upper)?.in_unit(parsed_unit)?;
    if !initial.is_finite()
        || !lower.is_finite()
        || !upper.is_finite()
        || lower <= 0.
        || lower >= upper
        || initial < lower
        || initial > upper
        || lower.powi(2) == 0.
        || initial.powi(2) == 0.
        || !upper.powi(2).is_finite()
    {
        return Err(invalid(
            "SAEM additive SD requires positive finite bounds containing its initial value",
        ));
    }
    let mut observations = 0usize;
    for (subject_index, subject) in request.subjects.iter().enumerate() {
        if !subject.priors.is_empty() {
            return Err(invalid("SAEM subject priors are unsupported"));
        }
        for (row, observation) in subject.observations.iter().enumerate() {
            let declared = quantity(&observation.error.additive_sd)?.in_unit(parsed_unit)?;
            if observation.output != fit.output
                || observation.error.proportional_sd != 0.
                || (declared - initial).abs() > 1e-10 * initial
            {
                return Err(invalid(format!(
                    "subjects[{subject_index}].observations[{row}] must use the fitted output, zero proportional error, and initial additive SD"
                )));
            }
            observations += 1;
        }
    }
    if observations == 0 {
        return Err(invalid("SAEM additive SD requires observations"));
    }
    Ok(Some(ResidualFit {
        output: fit.output.clone(),
        unit,
        variance: initial * initial,
        lower_variance: lower * lower,
        upper_variance: upper * upper,
        observations,
    }))
}

fn residual_sum_squares(nll: f64, observations: usize, variance: f64) -> Result<f64, Error> {
    let baseline = observations as f64 * 0.5 * (LOG_2PI + libm::log(variance));
    let excess = nll - baseline;
    if !excess.is_finite() || excess < -1e-9 * (1. + baseline.abs()) {
        return Err(Error::new(
            ErrorCode::Domain,
            "SAEM residual sum of squares is invalid",
        ));
    }
    Ok((2. * variance * excess).max(0.))
}

fn residual_nll(sum_squares: f64, observations: usize, variance: f64) -> f64 {
    0.5 * sum_squares / variance + observations as f64 * 0.5 * (LOG_2PI + libm::log(variance))
}

fn set_residual_sd(subjects: &mut [GaussianObjectiveRequest], fit: &ResidualFit) {
    let sd = libm::sqrt(fit.variance);
    for subject in subjects {
        for observation in &mut subject.observations {
            observation.error.additive_sd = pharmflux_core::model::Quantity {
                value: sd,
                unit: fit.unit.clone(),
            };
        }
    }
}

fn burn_in_iterations(request: &PopulationFitRequest) -> Result<usize, Error> {
    let burn_in = request
        .burn_in_iterations
        .unwrap_or((request.max_iterations / 3).max(4));
    if burn_in < 4 || burn_in >= request.max_iterations {
        return Err(invalid(
            "SAEM burn-in must be at least four iterations and below max_iterations",
        ));
    }
    Ok(burn_in)
}

fn precision(omega: Matrix, dimensions: usize) -> Option<(Matrix, f64)> {
    if dimensions == 1 {
        let variance = omega[0][0];
        return (variance.is_finite() && variance > 0.)
            .then_some(([[1. / variance, 0.], [0., 0.]], variance));
    }
    let det = omega[0][0] * omega[1][1] - omega[0][1] * omega[1][0];
    if !omega.iter().flatten().all(|x| x.is_finite())
        || omega[0][1] != omega[1][0]
        || omega[0][0] <= 0.
        || omega[1][1] <= 0.
        || det <= 1e-10 * omega[0][0] * omega[1][1]
    {
        return None;
    }
    Some((
        [
            [omega[1][1] / det, -omega[0][1] / det],
            [-omega[1][0] / det, omega[0][0] / det],
        ],
        det,
    ))
}

fn prior_nll_vector(eta: Vector, omega: Matrix, dimensions: usize) -> Result<f64, Error> {
    let (inverse, det) = precision(omega, dimensions)
        .ok_or_else(|| invalid("SAEM Omega must remain positive definite"))?;
    let quadratic = (0..dimensions)
        .map(|i| {
            (0..dimensions)
                .map(|j| eta[i] * inverse[i][j] * eta[j])
                .sum::<f64>()
        })
        .sum::<f64>();
    Ok(0.5 * (quadratic + libm::log(det) + dimensions as f64 * LOG_2PI))
}

struct Stream(u64);

impl Stream {
    fn new(seed: u64) -> Self {
        Self(seed)
    }

    fn next_u64(&mut self) -> u64 {
        // SplitMix64 has a fully specified stream across native and WASM.
        self.0 = self.0.wrapping_add(0x9e3779b97f4a7c15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d049bb133111eb);
        z ^ (z >> 31)
    }

    fn uniform(&mut self) -> f64 {
        // In (0, 1), so Box-Muller and log acceptance stay finite.
        ((self.next_u64() >> 11) as f64 + 0.5) * (1.0 / ((1_u64 << 53) as f64))
    }

    fn normal(&mut self) -> f64 {
        libm::sqrt(-2.0 * libm::log(self.uniform()))
            * libm::cos(2.0 * std::f64::consts::PI * self.uniform())
    }
}

fn prior_nll(phi: f64, log_theta: f64, omega: f64) -> f64 {
    let eta = phi - log_theta;
    0.5 * (eta * eta / omega + libm::log(omega)) + 0.9189385332046727
}

// Compare two halves of a post-stabilization history. Consecutive SAEM
// updates retain Monte Carlo noise even when the parameter trajectory is
// stationary, so a single-update threshold is not a useful stopping rule.
fn windowed_drift(values: &[f64], tolerance: f64) -> (f64, f64, bool) {
    let half = values.len() / 2;
    let (older, newer) = values.split_at(half);
    let old_mean = older.iter().sum::<f64>() / older.len() as f64;
    let new_mean = newer.iter().sum::<f64>() / newer.len() as f64;
    let variance_of_mean = |part: &[f64], mean: f64| {
        let centered = part.iter().map(|x| x - mean).collect::<Vec<_>>();
        let sum_squares = centered.iter().map(|x| x * x).sum::<f64>();
        let lag_one = centered
            .windows(2)
            .map(|pair| pair[0] * pair[1])
            .sum::<f64>();
        // SA trajectories are serially correlated. Inflate the variance of
        // each half-window mean by a capped AR(1) effective-sample factor.
        let correlation = if sum_squares > 0. {
            (lag_one / sum_squares).clamp(0., 0.9)
        } else {
            0.
        };
        sum_squares / ((part.len() - 1) * part.len()) as f64 * (1. + correlation)
            / (1. - correlation)
    };
    // Inputs are log parameters, so their difference approximates relative
    // change on the natural scale without a second coordinate-dependent scale.
    let drift = (new_mean - old_mean).abs();
    let noise =
        2.0 * libm::sqrt(variance_of_mean(older, old_mean) + variance_of_mean(newer, new_mean));
    // The noise allowance is capped: a wandering chain must not receive a
    // convergence label just because its Monte Carlo uncertainty is large.
    (drift, noise, drift <= tolerance.max(noise).min(0.005))
}

fn validate_two_effects(
    model: &Arc<CompiledSensitivityDocument>,
    request: &PopulationFitRequest,
) -> Result<Setup, Error> {
    if request.subjects.len() < 2
        || request.subjects.len() > 128
        || request.fixed_effects.len() != 2
        || request.random_effects.len() != 2
        || request.omega.len() != 2
        || request.omega.iter().any(|row| row.len() != 2)
        || request.omega_diagonal_bounds.len() != 2
        || model.parameters().len() != 2
    {
        return Err(Error::new(ErrorCode::Unsupported, "two-effect SAEM requires two selected fixed and random effects across 2..=128 subjects"));
    }
    if !request.eta_bound.is_finite()
        || request.eta_bound <= 0.
        || request.eta_bound > 8.
        || request.max_evaluations == 0
        || request.max_evaluations > 100_000
        || request.max_iterations < 8
        || request.max_iterations > 10_000
        || !request.gradient_tolerance.is_finite()
        || request.gradient_tolerance <= 0.
    {
        return Err(invalid("two-effect SAEM needs bounded positive controls"));
    }
    let mut setup = Setup {
        dimensions: 2,
        names: [model.parameters()[0].clone(), model.parameters()[1].clone()],
        units: [String::new(), String::new()],
        start: [0.; 2],
        lower: [0.; 2],
        upper: [0.; 2],
        omega: [
            [request.omega[0][0], request.omega[0][1]],
            [request.omega[1][0], request.omega[1][1]],
        ],
        omega_lower: [0.; 2],
        omega_upper: [0.; 2],
    };
    for axis in 0..2 {
        let name = &setup.names[axis];
        if request.fixed_effects[axis].parameter != *name
            || request.random_effects[axis].parameter != *name
        {
            return Err(invalid(
                "SAEM fixed and random effects must follow selected parameter order",
            ));
        }
        let parameter = model
            .primal
            .parameters
            .iter()
            .find(|p| p.name == *name)
            .ok_or_else(|| invalid("unknown SAEM parameter"))?;
        if parameter.fixed {
            return Err(invalid("cannot fit a fixed model parameter"));
        }
        let unit = Unit::parse(&parameter.default.unit)?;
        let lower = quantity(&request.fixed_effects[axis].lower)?.in_unit(unit)?;
        let upper = quantity(&request.fixed_effects[axis].upper)?.in_unit(unit)?;
        let start = request.subjects[0]
            .simulation
            .run
            .parameters
            .get(name)
            .map(|q| quantity(q)?.in_unit(unit))
            .transpose()?
            .unwrap_or(quantity(&parameter.default)?.in_unit(unit)?);
        if !lower.is_finite()
            || !upper.is_finite()
            || !start.is_finite()
            || lower <= 0.
            || lower >= upper
            || start < lower
            || start > upper
        {
            return Err(invalid(
                "SAEM lognormal fixed effects need positive finite bounds containing their starts",
            ));
        }
        if let Some(bounds) = &parameter.bounds {
            if lower * libm::exp(-request.eta_bound) < quantity(&bounds[0])?.in_unit(unit)?
                || upper * libm::exp(request.eta_bound) > quantity(&bounds[1])?.in_unit(unit)?
            {
                return Err(invalid(
                    "SAEM fixed-effect and eta boxes exceed model bounds",
                ));
            }
        }
        setup.units[axis] = parameter.default.unit.clone();
        setup.start[axis] = start;
        setup.lower[axis] = lower;
        setup.upper[axis] = upper;
        let [variance_lower, variance_upper] = request.omega_diagonal_bounds[axis];
        if !variance_lower.is_finite()
            || !variance_upper.is_finite()
            || variance_lower <= 0.
            || variance_lower > variance_upper
            || setup.omega[axis][axis] < variance_lower
            || setup.omega[axis][axis] > variance_upper
        {
            return Err(invalid(
                "SAEM Omega variance bounds must be finite and positive",
            ));
        }
        setup.omega_lower[axis] = variance_lower;
        setup.omega_upper[axis] = variance_upper;
    }
    if precision(setup.omega, 2).is_none() {
        return Err(invalid("SAEM Omega must be symmetric positive definite"));
    }
    for (index, subject) in request.subjects.iter().enumerate() {
        if !subject.priors.is_empty() || subject.simulation.with_respect_to != model.parameters() {
            return Err(invalid(format!(
                "subject {index} must share selected parameters and use no priors"
            )));
        }
        let bound = model.bind_with_initial_states(
            &subject.simulation.run.parameters,
            &subject.simulation.run.regimen.covariates,
            &subject.simulation.run.initial_states,
        )?;
        for axis in 0..2 {
            let parameter_index = model
                .primal
                .parameters
                .iter()
                .position(|p| p.name == setup.names[axis])
                .unwrap();
            if (bound.primal.values[parameter_index] - setup.start[axis]).abs()
                > 1e-10 * setup.start[axis].abs().max(1.)
            {
                return Err(invalid(format!(
                    "subject {index} has a conflicting fixed-effect start"
                )));
            }
        }
    }
    Ok(setup)
}

impl CompiledSensitivityDocument {
    fn saem_nll(
        self: &Arc<Self>,
        subject: &GaussianObjectiveRequest,
        setup: &Setup,
        theta: Vector,
        eta: Vector,
    ) -> Result<f64, Error> {
        let mut objective = subject.clone();
        self.saem_nll_in_place(&mut objective, setup, theta, eta)
    }

    // Every selected parameter is overwritten before the solve, so a rejected
    // proposal cannot leave parameter values that affect the next evaluation.
    fn saem_nll_in_place(
        self: &Arc<Self>,
        objective: &mut GaussianObjectiveRequest,
        setup: &Setup,
        theta: Vector,
        eta: Vector,
    ) -> Result<f64, Error> {
        for axis in 0..setup.dimensions {
            let individual = theta[axis] * libm::exp(eta[axis]);
            if !individual.is_finite() || individual <= 0. {
                return Err(Error::new(
                    ErrorCode::Domain,
                    "invalid SAEM individual parameter",
                ));
            }
            objective.simulation.run.parameters.insert(
                setup.names[axis].clone(),
                pharmflux_core::model::Quantity {
                    value: individual,
                    unit: setup.units[axis].clone(),
                },
            );
        }
        self.gaussian_primal_nll(objective)
    }

    fn saem_predictions(
        self: &Arc<Self>,
        subject_index: usize,
        subject: &GaussianObjectiveRequest,
        setup: &Setup,
        theta: Vector,
        eta: Vector,
    ) -> Result<Vec<PopulationFittedObservation>, Error> {
        let mut population = subject.simulation.run.clone();
        let mut individual = subject.simulation.run.clone();
        for axis in 0..setup.dimensions {
            let name = setup.names[axis].clone();
            let unit = setup.units[axis].clone();
            population.parameters.insert(
                name.clone(),
                pharmflux_core::model::Quantity {
                    value: theta[axis],
                    unit: unit.clone(),
                },
            );
            individual.parameters.insert(
                name,
                pharmflux_core::model::Quantity {
                    value: theta[axis] * libm::exp(eta[axis]),
                    unit,
                },
            );
        }
        let pred = self.primal.execute(&population)?;
        let ipred = self.primal.execute(&individual)?;
        if pred.times != ipred.times || pred.sides != ipred.sides {
            return Err(Error::new(
                ErrorCode::Invariant,
                "SAEM prediction grids differ across eta",
            ));
        }
        subject
            .observations
            .iter()
            .enumerate()
            .map(|(observation_index, observation)| {
                let output = self
                    .output_names()
                    .iter()
                    .position(|name| name == &observation.output)
                    .ok_or_else(|| invalid("unknown SAEM observation output"))?;
                let unit = self.output_units()[output].clone();
                let parsed = Unit::parse(&unit)?;
                let value = quantity(&observation.value)?.in_unit(parsed)?;
                let prediction = *pred.outputs[output]
                    .values
                    .get(observation.row)
                    .ok_or_else(|| {
                        invalid("SAEM observation row is outside the simulation result")
                    })?;
                let individual_prediction = *ipred.outputs[output]
                    .values
                    .get(observation.row)
                    .ok_or_else(|| {
                        invalid("SAEM observation row is outside the simulation result")
                    })?;
                let time = *pred.times.get(observation.row).ok_or_else(|| {
                    invalid("SAEM observation row is outside the simulation result")
                })?;
                let side = pred.sides[observation.row];
                let q = |value| pharmflux_core::model::Quantity {
                    value,
                    unit: unit.clone(),
                };
                Ok(PopulationFittedObservation {
                    subject_index,
                    observation_index,
                    row: observation.row,
                    output: observation.output.clone(),
                    time: pharmflux_core::model::Quantity {
                        value: time,
                        unit: self.primal.time_unit.clone(),
                    },
                    side,
                    value: q(value),
                    pred: q(prediction),
                    ipred: q(individual_prediction),
                })
            })
            .collect()
    }

    fn saem_joint(
        self: &Arc<Self>,
        subject: &GaussianObjectiveRequest,
        setup: &Setup,
        theta: Vector,
        omega: Matrix,
        eta: Vector,
        evaluations: &mut usize,
        budget: usize,
    ) -> Result<Option<f64>, Error> {
        if *evaluations >= budget {
            return Ok(None);
        }
        if eta[..setup.dimensions]
            .iter()
            .any(|x| !x.is_finite() || x.abs() >= 0.999 * 12.)
        {
            return Ok(None);
        }
        *evaluations += 1;
        let nll = match self.saem_nll(subject, setup, theta, eta) {
            Ok(value) => value,
            Err(error)
                if matches!(
                    error.code,
                    ErrorCode::Domain | ErrorCode::Solver | ErrorCode::Invariant
                ) =>
            {
                return Ok(None)
            }
            Err(error) => return Err(error),
        };
        Ok(Some(nll + prior_nll_vector(eta, omega, setup.dimensions)?))
    }

    /// Local observed-curvature Laplace approximation at the final SAEM parameters.
    /// `None` means its conditional mode or curvature could not be certified
    /// within the remaining work budget; the caller keeps a labeled surrogate.
    fn saem_marginal_ofv(
        self: &Arc<Self>,
        request: &PopulationFitRequest,
        setup: &Setup,
        theta: Vector,
        omega: Matrix,
        etas: &[Vector],
        evaluations: &mut usize,
        budget: usize,
    ) -> Result<Option<(f64, Vec<Vector>)>, Error> {
        let dimensions = setup.dimensions;
        let mut marginal = 0.;
        let mut modes = Vec::with_capacity(request.subjects.len());
        for (index, subject) in request.subjects.iter().enumerate() {
            let mut eta = etas[index];
            let mut best =
                match self.saem_joint(subject, setup, theta, omega, eta, evaluations, budget)? {
                    Some(value) => value,
                    None => return Ok(None),
                };
            let h = 0.005;
            for _ in 0..16 {
                let mut gradient = [0.; 2];
                let mut hessian = [[0.; 2]; 2];
                for axis in 0..dimensions {
                    let mut trial = eta;
                    trial[axis] += h;
                    if trial[axis] >= request.eta_bound {
                        return Ok(None);
                    }
                    let plus = match self.saem_joint(
                        subject,
                        setup,
                        theta,
                        omega,
                        trial,
                        evaluations,
                        budget,
                    )? {
                        Some(value) => value,
                        None => return Ok(None),
                    };
                    trial[axis] -= 2. * h;
                    if trial[axis] <= -request.eta_bound {
                        return Ok(None);
                    }
                    let minus = match self.saem_joint(
                        subject,
                        setup,
                        theta,
                        omega,
                        trial,
                        evaluations,
                        budget,
                    )? {
                        Some(value) => value,
                        None => return Ok(None),
                    };
                    gradient[axis] = (plus - minus) / (2. * h);
                    hessian[axis][axis] = (plus - 2. * best + minus) / (h * h);
                }
                if dimensions == 2 {
                    let mut corners = [0.; 4];
                    for (i, signs) in [[1., 1.], [1., -1.], [-1., 1.], [-1., -1.]]
                        .iter()
                        .enumerate()
                    {
                        let trial = [eta[0] + h * signs[0], eta[1] + h * signs[1]];
                        corners[i] = match self.saem_joint(
                            subject,
                            setup,
                            theta,
                            omega,
                            trial,
                            evaluations,
                            budget,
                        )? {
                            Some(value) => value,
                            None => return Ok(None),
                        };
                    }
                    hessian[0][1] =
                        (corners[0] - corners[1] - corners[2] + corners[3]) / (4. * h * h);
                    hessian[1][0] = hessian[0][1];
                }
                let determinant = if dimensions == 1 {
                    hessian[0][0]
                } else {
                    hessian[0][0] * hessian[1][1] - hessian[0][1] * hessian[1][0]
                };
                if !determinant.is_finite() || hessian[0][0] <= 0. || determinant <= 0. {
                    return Ok(None);
                }
                let mut next = eta;
                if gradient[..dimensions].iter().all(|x| x.abs() < 1e-4) {
                    break;
                }
                if dimensions == 1 {
                    next[0] -= gradient[0] / determinant;
                } else {
                    next[0] -=
                        (hessian[1][1] * gradient[0] - hessian[0][1] * gradient[1]) / determinant;
                    next[1] -=
                        (-hessian[1][0] * gradient[0] + hessian[0][0] * gradient[1]) / determinant;
                }
                for axis in 0..dimensions {
                    next[axis] =
                        next[axis].clamp(-request.eta_bound + 2. * h, request.eta_bound - 2. * h);
                }
                if next == eta {
                    break;
                }
                let mut accepted = false;
                for _ in 0..12 {
                    let candidate = match self.saem_joint(
                        subject,
                        setup,
                        theta,
                        omega,
                        next,
                        evaluations,
                        budget,
                    )? {
                        Some(value) => value,
                        None => return Ok(None),
                    };
                    if candidate < best {
                        best = candidate;
                        accepted = true;
                        break;
                    }
                    for axis in 0..dimensions {
                        next[axis] = 0.5 * (eta[axis] + next[axis]);
                    }
                }
                if !accepted {
                    break;
                }
                eta = next;
            }
            let mut hessian = [[0.; 2]; 2];
            let mut gradient = [0.; 2];
            for axis in 0..dimensions {
                let mut trial = eta;
                trial[axis] += h;
                let plus = match self.saem_joint(
                    subject,
                    setup,
                    theta,
                    omega,
                    trial,
                    evaluations,
                    budget,
                )? {
                    Some(value) => value,
                    None => return Ok(None),
                };
                trial[axis] -= 2. * h;
                let minus = match self.saem_joint(
                    subject,
                    setup,
                    theta,
                    omega,
                    trial,
                    evaluations,
                    budget,
                )? {
                    Some(value) => value,
                    None => return Ok(None),
                };
                gradient[axis] = (plus - minus) / (2. * h);
                hessian[axis][axis] = (plus - 2. * best + minus) / (h * h);
            }
            if dimensions == 2 {
                let mut corners = [0.; 4];
                for (i, signs) in [[1., 1.], [1., -1.], [-1., 1.], [-1., -1.]]
                    .iter()
                    .enumerate()
                {
                    let trial = [eta[0] + h * signs[0], eta[1] + h * signs[1]];
                    corners[i] = match self.saem_joint(
                        subject,
                        setup,
                        theta,
                        omega,
                        trial,
                        evaluations,
                        budget,
                    )? {
                        Some(value) => value,
                        None => return Ok(None),
                    };
                }
                hessian[0][1] = (corners[0] - corners[1] - corners[2] + corners[3]) / (4. * h * h);
                hessian[1][0] = hessian[0][1];
            }
            let det_hessian = if dimensions == 1 {
                hessian[0][0]
            } else {
                hessian[0][0] * hessian[1][1] - hessian[0][1] * hessian[1][0]
            };
            let mode_step = if dimensions == 1 {
                [gradient[0] / det_hessian, 0.]
            } else {
                [
                    (hessian[1][1] * gradient[0] - hessian[0][1] * gradient[1]) / det_hessian,
                    (-hessian[1][0] * gradient[0] + hessian[0][0] * gradient[1]) / det_hessian,
                ]
            };
            let newton_decrement = gradient[0] * mode_step[0] + gradient[1] * mode_step[1];
            if !det_hessian.is_finite()
                || hessian[0][0] <= 0.
                || det_hessian <= 0.
                || mode_step[..dimensions]
                    .iter()
                    .any(|x| !x.is_finite() || x.abs() > h / 10.)
                || !newton_decrement.is_finite()
                || newton_decrement > 1e-3
                || eta[..dimensions]
                    .iter()
                    .any(|x| x.abs() >= request.eta_bound - h)
            {
                return Ok(None);
            }
            // Integrating the normal eta prior cancels its (2π) factor with
            // the Laplace factor; remove only observation constants.
            let observation_constant = subject.observations.len() as f64 * LOG_2PI;
            marginal += 2. * best + libm::log(det_hessian)
                - dimensions as f64 * LOG_2PI
                - observation_constant;
            modes.push(eta);
        }
        Ok(marginal.is_finite().then_some((marginal, modes)))
    }
}

impl CompiledSensitivityDocument {
    /// Stochastic-approximation EM with a Metropolis conditional simulation.
    /// Estimates one or two positive lognormal effects with fixed Gaussian
    /// error or a common fitted additive Gaussian SD.
    pub fn fit_population_saem(
        self: &Arc<Self>,
        request: &PopulationFitRequest,
    ) -> Result<PopulationFitResult, Error> {
        if request.uncertainty {
            return Err(Error::new(
                ErrorCode::Unsupported,
                "SAEM uncertainty is not implemented",
            ));
        }
        if request.random_effects.len() == 2 {
            return self.fit_population_saem_2d(request);
        }
        let cfg = self.validate_population_1d(request)?;
        let burn_in = burn_in_iterations(request)?;
        let mut residual = validate_residual_fit(self, request)?;
        let mut subjects = request.subjects.clone();
        let setup = Setup {
            dimensions: 1,
            names: [cfg.parameter.clone(), String::new()],
            units: [cfg.unit.clone(), String::new()],
            start: [cfg.start, 0.],
            lower: [cfg.lower, 0.],
            upper: [cfg.upper, 0.],
            omega: [[cfg.omega, 0.], [0., 0.]],
            omega_lower: [cfg.omega_lower, 0.],
            omega_upper: [cfg.omega_upper, 0.],
        };
        let count = request.subjects.len();
        if count < 2 {
            return Err(invalid("SAEM requires at least two subjects"));
        }
        if request.max_iterations < 8 || request.max_evaluations < count.saturating_mul(4) {
            return Err(invalid(
                "SAEM requires at least eight iterations and four evaluations per subject",
            ));
        }
        let mut stream = Stream::new(request.seed);
        let mut theta = cfg.start;
        let mut omega = cfg.omega;
        let mut log_theta = libm::log(theta);
        let mut phi = vec![log_theta; count];
        let mut nll = Vec::with_capacity(count);
        let mut residual_sse = Vec::with_capacity(count);
        let mut evaluations = 0_usize;
        for subject in &subjects {
            if evaluations >= request.max_evaluations {
                return Err(invalid(
                    "SAEM evaluation budget is smaller than subject count",
                ));
            }
            let result =
                self.population_subject_primal_nll(subject, &cfg.parameter, theta, &cfg.unit, 0.0)?;
            evaluations += 1;
            if let Some(fit) = &residual {
                residual_sse.push(residual_sum_squares(
                    result,
                    subject.observations.len(),
                    fit.variance,
                )?);
            }
            nll.push(result);
        }

        // The complete-data family is normal in phi = log(individual parameter).
        // SA updates these sufficient statistics; the M-step is their exact
        // constrained normal-family maximizer. A fitted additive error uses
        // a separate stochastic average of the complete-data residual SSE.
        let mut mean_phi = log_theta;
        let mut second_phi = log_theta * log_theta + omega;
        let mut iterations = 0_usize;
        let mut stable_iterations = 0_usize;
        let mut evaluated_proposals = 0_usize;
        let mut accepted_proposals = 0_usize;
        let mut final_step_size = None;
        let mut final_relative_change = None;
        let mut recent_relative_changes = Vec::with_capacity(20);
        let mut post_burn_history = Vec::with_capacity(50);
        let mut status = FitStatus::IterationLimit;
        for iteration in 0..request.max_iterations {
            let mut moved = 0_usize;
            for (index, subject) in subjects.iter().enumerate() {
                for _ in 0..3 {
                    let proposal = phi[index] + 0.8 * libm::sqrt(omega) * stream.normal();
                    let proposal_eta = proposal - log_theta;
                    if !proposal.is_finite() || proposal_eta.abs() > request.eta_bound {
                        continue;
                    }
                    let individual = libm::exp(proposal);
                    if !individual.is_finite() || individual <= 0.0 {
                        continue;
                    }
                    // Reserve the final likelihood plus population and
                    // individual prediction runs for every subject.
                    if evaluations >= request.max_evaluations - 3 * count {
                        status = FitStatus::EvaluationLimit;
                        break;
                    }
                    let proposal_result = self.population_subject_primal_nll(
                        subject,
                        &cfg.parameter,
                        theta,
                        &cfg.unit,
                        proposal_eta,
                    );
                    evaluations += 1;
                    evaluated_proposals += 1;
                    let proposal_nll = match proposal_result {
                        Ok(value) => value,
                        Err(error)
                            if matches!(
                                error.code,
                                ErrorCode::Domain | ErrorCode::Solver | ErrorCode::Invariant
                            ) =>
                        {
                            continue
                        }
                        Err(error) => return Err(error),
                    };
                    let proposal_sse = if let Some(fit) = &residual {
                        Some(residual_sum_squares(
                            proposal_nll,
                            subject.observations.len(),
                            fit.variance,
                        )?)
                    } else {
                        None
                    };
                    let current = nll[index] + prior_nll(phi[index], log_theta, omega);
                    let candidate = proposal_nll + prior_nll(proposal, log_theta, omega);
                    if libm::log(stream.uniform()) < current - candidate {
                        phi[index] = proposal;
                        nll[index] = proposal_nll;
                        if let Some(sse) = proposal_sse {
                            residual_sse[index] = sse;
                        }
                        moved += 1;
                        accepted_proposals += 1;
                    }
                }
                if status == FitStatus::EvaluationLimit {
                    break;
                }
            }
            if status == FitStatus::EvaluationLimit {
                break;
            }
            let sample_mean = phi.iter().sum::<f64>() / count as f64;
            let sample_second = phi.iter().map(|p| p * p).sum::<f64>() / count as f64;
            let gamma = if iteration < burn_in {
                1.0
            } else {
                1.0 / (iteration - burn_in + 1) as f64
            };
            mean_phi += gamma * (sample_mean - mean_phi);
            second_phi += gamma * (sample_second - second_phi);
            let log_lower = libm::log(cfg.lower)
                .max(phi.iter().copied().fold(f64::NEG_INFINITY, f64::max) - request.eta_bound);
            let log_upper = libm::log(cfg.upper)
                .min(phi.iter().copied().fold(f64::INFINITY, f64::min) + request.eta_bound);
            if log_lower > log_upper {
                return Err(Error::new(
                    ErrorCode::Domain,
                    "SAEM eta constraints have no feasible fixed effect",
                ));
            }
            let next_log_theta = mean_phi.clamp(log_lower, log_upper);
            let next_theta = libm::exp(next_log_theta);
            let next_omega = (second_phi - 2.0 * next_log_theta * mean_phi
                + next_log_theta * next_log_theta)
                .clamp(cfg.omega_lower, cfg.omega_upper);
            if !next_theta.is_finite() || !next_omega.is_finite() || next_omega <= 0.0 {
                return Err(Error::new(
                    ErrorCode::Domain,
                    "SAEM M-step produced invalid population parameters",
                ));
            }
            let residual_change = if let Some(fit) = &mut residual {
                let previous = fit.variance;
                let sample = residual_sse.iter().sum::<f64>() / fit.observations as f64;
                // A single pre-equilibration draw must not overwrite the
                // residual variance and feed that noise into the next E-step.
                let residual_gain = gamma.min(0.2);
                fit.variance = (previous + residual_gain * (sample - previous))
                    .clamp(fit.lower_variance, fit.upper_variance);
                if !fit.variance.is_finite() {
                    return Err(Error::new(
                        ErrorCode::Domain,
                        "SAEM residual M-step overflow",
                    ));
                }
                for (index, subject) in subjects.iter().enumerate() {
                    nll[index] = residual_nll(
                        residual_sse[index],
                        subject.observations.len(),
                        fit.variance,
                    );
                }
                set_residual_sd(&mut subjects, fit);
                ((fit.variance - previous) / previous).abs()
            } else {
                0.
            };
            let relative_change = ((next_theta - theta) / theta)
                .abs()
                .max(((next_omega - omega) / omega).abs())
                .max(residual_change);
            final_step_size = Some(gamma);
            final_relative_change = Some(relative_change);
            if iteration >= burn_in {
                if recent_relative_changes.len() == 20 {
                    recent_relative_changes.remove(0);
                }
                recent_relative_changes.push(relative_change);
                if post_burn_history.len() == 50 {
                    post_burn_history.remove(0);
                }
                post_burn_history.push((
                    next_log_theta,
                    libm::log(next_omega),
                    residual
                        .as_ref()
                        .map(|fit| libm::log(fit.variance))
                        .unwrap_or(0.),
                ));
            }
            theta = next_theta;
            omega = next_omega;
            log_theta = libm::log(theta);
            iterations = iteration + 1;
            let window_stable = if post_burn_history.len() == 50 && moved > 0 {
                let theta_history: Vec<f64> = post_burn_history.iter().map(|x| x.0).collect();
                let omega_history: Vec<f64> = post_burn_history.iter().map(|x| x.1).collect();
                let residual_history: Vec<f64> = post_burn_history.iter().map(|x| x.2).collect();
                let separated_theta: Vec<f64> = theta_history[..15]
                    .iter()
                    .chain(&theta_history[35..])
                    .copied()
                    .collect();
                let separated_omega: Vec<f64> = omega_history[..15]
                    .iter()
                    .chain(&omega_history[35..])
                    .copied()
                    .collect();
                let separated_residual: Vec<f64> = residual_history[..15]
                    .iter()
                    .chain(&residual_history[35..])
                    .copied()
                    .collect();
                windowed_drift(&theta_history[25..], request.gradient_tolerance).2
                    && windowed_drift(&omega_history[25..], request.gradient_tolerance).2
                    && windowed_drift(&residual_history[25..], request.gradient_tolerance).2
                    && windowed_drift(&separated_theta, request.gradient_tolerance).2
                    && windowed_drift(&separated_omega, request.gradient_tolerance).2
                    && windowed_drift(&separated_residual, request.gradient_tolerance).2
            } else {
                false
            };
            stable_iterations = if window_stable {
                stable_iterations + 1
            } else {
                0
            };
            if stable_iterations >= 4 {
                status = FitStatus::Converged;
                break;
            }
        }

        // A converged marginal calculation also supplies conditional eta
        // modes for EBE and IPRED. At a tight work budget the sampled etas
        // remain available with an explicitly labeled complete-data surrogate.
        let final_request = residual.as_ref().map(|_| {
            let mut final_request = request.clone();
            final_request.subjects = subjects.clone();
            final_request
        });
        let fitted_request = final_request.as_ref().unwrap_or(request);
        let sampled_etas: Vec<Vector> = phi.iter().map(|p| [p - log_theta, 0.]).collect();
        let marginal = self.saem_marginal_ofv(
            fitted_request,
            &setup,
            [theta, 0.],
            [[omega, 0.], [0., 0.]],
            &sampled_etas,
            &mut evaluations,
            request.max_evaluations - 3 * count,
        )?;
        if marginal.is_none() && status == FitStatus::Converged {
            status = if evaluations >= request.max_evaluations - 3 * count {
                FitStatus::EvaluationLimit
            } else {
                FitStatus::LineSearchFailed
            };
        }
        let final_etas = marginal
            .as_ref()
            .map(|(_, modes)| modes)
            .unwrap_or(&sampled_etas);
        let mut subject_results = Vec::with_capacity(count);
        let mut fitted_observations = Vec::new();
        let mut objective = 0.0;
        for (index, subject) in fitted_request.subjects.iter().enumerate() {
            let eta = final_etas[index][0];
            let result =
                self.population_subject_primal_nll(subject, &cfg.parameter, theta, &cfg.unit, eta)?;
            evaluations += 1;
            let likelihood = result;
            objective += likelihood + prior_nll_vector([eta, 0.], [[omega, 0.], [0., 0.]], 1)?;
            subject_results.push(PopulationSubjectEstimate {
                eta: vec![eta],
                negative_log_likelihood: likelihood,
            });
            fitted_observations.extend(self.saem_predictions(
                index,
                subject,
                &setup,
                [theta, 0.],
                [eta, 0.],
            )?);
            evaluations += 2;
        }
        if !objective.is_finite() {
            return Err(Error::new(
                ErrorCode::Domain,
                "SAEM complete-data surrogate overflow",
            ));
        }
        let mut fixed_effects = std::collections::BTreeMap::new();
        fixed_effects.insert(
            cfg.parameter,
            pharmflux_core::model::Quantity {
                value: theta,
                unit: cfg.unit,
            },
        );
        let recent_window = recent_relative_changes.len();
        recent_relative_changes.sort_by(f64::total_cmp);
        let recent_relative_change_median = if recent_window == 0 {
            None
        } else if recent_window % 2 == 0 {
            Some(
                (recent_relative_changes[recent_window / 2 - 1]
                    + recent_relative_changes[recent_window / 2])
                    / 2.0,
            )
        } else {
            Some(recent_relative_changes[recent_window / 2])
        };
        Ok(PopulationFitResult {
            status,
            fixed_effects,
            omega: vec![vec![omega]],
            subjects: subject_results,
            fitted_observations: Some(fitted_observations),
            uncertainty: None,
            estimated_residual_error: residual.map(|fit| PopulationResidualErrorEstimate {
                output: fit.output,
                additive_sd: pharmflux_core::model::Quantity {
                    value: libm::sqrt(fit.variance),
                    unit: fit.unit,
                },
            }),
            saem_diagnostics: Some(PopulationSaemDiagnostics {
                burn_in_iterations: burn_in,
                final_step_size,
                final_relative_change,
                recent_relative_change_median,
                recent_window,
                evaluated_proposals,
                accepted_proposals,
                consecutive_stable_iterations: stable_iterations,
            }),
            objective_kind: if marginal.is_some() {
                PopulationObjectiveKind::SaemMarginalObjectiveFunctionValue
            } else {
                PopulationObjectiveKind::SaemCompleteDataSurrogate
            },
            objective: marginal.map(|(value, _)| value).unwrap_or(objective),
            evaluations,
            iterations,
            fit_request_hash: canonical_hash(&serde_json::json!(request))?,
        })
    }

    fn fit_population_saem_2d(
        self: &Arc<Self>,
        request: &PopulationFitRequest,
    ) -> Result<PopulationFitResult, Error> {
        let setup = validate_two_effects(self, request)?;
        let burn_in = burn_in_iterations(request)?;
        let mut residual = validate_residual_fit(self, request)?;
        let mut subjects = request.subjects.clone();
        let count = request.subjects.len();
        if request.max_evaluations < count.saturating_mul(4) {
            return Err(invalid(
                "SAEM requires four evaluations per subject for complete results",
            ));
        }
        let mut stream = Stream::new(request.seed);
        let mut theta = setup.start;
        let mut omega = setup.omega;
        let mut log_theta = [libm::log(theta[0]), libm::log(theta[1])];
        let mut phi = vec![log_theta; count];
        let mut nll = Vec::with_capacity(count);
        let mut residual_sse = Vec::with_capacity(count);
        let mut evaluations = 0;
        for subject in &mut subjects {
            let likelihood = self.saem_nll_in_place(subject, &setup, theta, [0.; 2])?;
            if let Some(fit) = &residual {
                residual_sse.push(residual_sum_squares(
                    likelihood,
                    subject.observations.len(),
                    fit.variance,
                )?);
            }
            nll.push(likelihood);
            evaluations += 1;
        }
        let mut mean_phi = log_theta;
        let mut second_phi = [[0.; 2]; 2];
        for a in 0..2 {
            for b in 0..2 {
                second_phi[a][b] = log_theta[a] * log_theta[b] + omega[a][b];
            }
        }
        let mut iterations = 0;
        let mut status = FitStatus::IterationLimit;
        let mut evaluated_proposals = 0;
        let mut accepted_proposals = 0;
        let mut final_step_size = None;
        let mut final_relative_change = None;
        let mut recent_relative_changes = Vec::with_capacity(20);
        let mut post_burn_history: Vec<[f64; 6]> = Vec::with_capacity(50);
        let mut stable_iterations = 0;
        let mut proposal_scale = [0.8_f64; 2];
        for iteration in 0..request.max_iterations {
            let mut moved = 0;
            let mut proposed_by_axis = [0_usize; 2];
            let mut accepted_by_axis = [0_usize; 2];
            for (index, subject) in subjects.iter_mut().enumerate() {
                let proposals = if residual.is_some() { 5 } else { 3 };
                for step in 0..proposals {
                    let mut proposal = phi[index];
                    // A symmetric coordinate proposal is valid even for
                    // correlated Omega; the acceptance ratio uses its full density.
                    let axis = if residual.is_some() {
                        (iteration + step) % 2
                    } else {
                        (stream.next_u64() % 2) as usize
                    };
                    proposal[axis] +=
                        proposal_scale[axis] * libm::sqrt(omega[axis][axis]) * stream.normal();
                    proposed_by_axis[axis] += 1;
                    let candidate_eta = [proposal[0] - log_theta[0], proposal[1] - log_theta[1]];
                    if candidate_eta
                        .iter()
                        .any(|x| !x.is_finite() || x.abs() > request.eta_bound)
                    {
                        continue;
                    }
                    if evaluations >= request.max_evaluations - 3 * count {
                        status = FitStatus::EvaluationLimit;
                        break;
                    }
                    evaluations += 1;
                    evaluated_proposals += 1;
                    let candidate_nll =
                        match self.saem_nll_in_place(subject, &setup, theta, candidate_eta) {
                            Ok(value) => value,
                            Err(error)
                                if matches!(
                                    error.code,
                                    ErrorCode::Domain | ErrorCode::Solver | ErrorCode::Invariant
                                ) =>
                            {
                                continue
                            }
                            Err(error) => return Err(error),
                        };
                    let candidate_sse = if let Some(fit) = &residual {
                        Some(residual_sum_squares(
                            candidate_nll,
                            subject.observations.len(),
                            fit.variance,
                        )?)
                    } else {
                        None
                    };
                    let current_eta = [phi[index][0] - log_theta[0], phi[index][1] - log_theta[1]];
                    let current = nll[index] + prior_nll_vector(current_eta, omega, 2)?;
                    let candidate = candidate_nll + prior_nll_vector(candidate_eta, omega, 2)?;
                    if libm::log(stream.uniform()) < current - candidate {
                        phi[index] = proposal;
                        nll[index] = candidate_nll;
                        if let Some(sse) = candidate_sse {
                            residual_sse[index] = sse;
                        }
                        accepted_proposals += 1;
                        accepted_by_axis[axis] += 1;
                        moved += 1;
                    }
                }
                if status == FitStatus::EvaluationLimit {
                    break;
                }
            }
            if status == FitStatus::EvaluationLimit {
                break;
            }
            // Tune each coordinate while the gain is one, then freeze the
            // symmetric proposal for the averaging phase. A single scale can
            // leave the tighter posterior coordinate nearly immobile.
            if iteration < burn_in {
                let gain = 1. / libm::pow((iteration + 1) as f64, 0.6);
                for axis in 0..2 {
                    if proposed_by_axis[axis] > 0 {
                        let rate = accepted_by_axis[axis] as f64 / proposed_by_axis[axis] as f64;
                        proposal_scale[axis] =
                            libm::exp(libm::log(proposal_scale[axis]) + gain * (rate - 0.44))
                                .clamp(0.01, 5.);
                    }
                }
            }
            let gamma = if iteration < burn_in {
                1.
            } else {
                1. / (iteration - burn_in + 1) as f64
            };
            for axis in 0..2 {
                let sample_mean = phi.iter().map(|p| p[axis]).sum::<f64>() / count as f64;
                mean_phi[axis] += gamma * (sample_mean - mean_phi[axis]);
                for other in 0..2 {
                    let sample_second =
                        phi.iter().map(|p| p[axis] * p[other]).sum::<f64>() / count as f64;
                    second_phi[axis][other] += gamma * (sample_second - second_phi[axis][other]);
                }
            }
            let mut next_log_theta = [0.; 2];
            let mut next_theta = [0.; 2];
            for axis in 0..2 {
                let lower = libm::log(setup.lower[axis]).max(
                    phi.iter()
                        .map(|p| p[axis])
                        .fold(f64::NEG_INFINITY, f64::max)
                        - request.eta_bound,
                );
                let upper = libm::log(setup.upper[axis]).min(
                    phi.iter().map(|p| p[axis]).fold(f64::INFINITY, f64::min) + request.eta_bound,
                );
                if lower > upper {
                    return Err(Error::new(
                        ErrorCode::Domain,
                        "SAEM eta constraints have no feasible fixed effect",
                    ));
                }
                next_log_theta[axis] = mean_phi[axis].clamp(lower, upper);
                next_theta[axis] = libm::exp(next_log_theta[axis]);
            }
            let mut next_omega = [[0.; 2]; 2];
            for a in 0..2 {
                for b in 0..2 {
                    next_omega[a][b] = second_phi[a][b]
                        - next_log_theta[a] * mean_phi[b]
                        - next_log_theta[b] * mean_phi[a]
                        + next_log_theta[a] * next_log_theta[b];
                }
            }
            for axis in 0..2 {
                next_omega[axis][axis] =
                    next_omega[axis][axis].clamp(setup.omega_lower[axis], setup.omega_upper[axis]);
            }
            let correlation_limit = if residual.is_some() && iteration < burn_in {
                0.75
            } else {
                0.95
            };
            let covariance_limit =
                correlation_limit * libm::sqrt(next_omega[0][0] * next_omega[1][1]);
            let covariance = (0.5 * (next_omega[0][1] + next_omega[1][0]))
                .clamp(-covariance_limit, covariance_limit);
            next_omega[0][1] = covariance;
            next_omega[1][0] = covariance;
            if next_theta.iter().any(|x| !x.is_finite()) || precision(next_omega, 2).is_none() {
                return Err(Error::new(
                    ErrorCode::Domain,
                    "SAEM M-step produced invalid population parameters",
                ));
            }
            let mut change: f64 = 0.;
            for axis in 0..2 {
                change = change
                    .max(((next_theta[axis] - theta[axis]) / theta[axis]).abs())
                    .max(((next_omega[axis][axis] - omega[axis][axis]) / omega[axis][axis]).abs());
            }
            change = change.max(
                (next_omega[0][1] - omega[0][1]).abs() / libm::sqrt(omega[0][0] * omega[1][1]),
            );
            if let Some(fit) = &mut residual {
                let previous = fit.variance;
                let sample = residual_sse.iter().sum::<f64>() / fit.observations as f64;
                let residual_gain = gamma.min(0.2);
                fit.variance = (previous + residual_gain * (sample - previous))
                    .clamp(fit.lower_variance, fit.upper_variance);
                if !fit.variance.is_finite() {
                    return Err(Error::new(
                        ErrorCode::Domain,
                        "SAEM residual M-step overflow",
                    ));
                }
                change = change.max(((fit.variance - previous) / previous).abs());
                for (index, subject) in subjects.iter().enumerate() {
                    nll[index] = residual_nll(
                        residual_sse[index],
                        subject.observations.len(),
                        fit.variance,
                    );
                }
                set_residual_sd(&mut subjects, fit);
            }
            final_step_size = Some(gamma);
            final_relative_change = Some(change);
            if iteration >= burn_in {
                if recent_relative_changes.len() == 20 {
                    recent_relative_changes.remove(0);
                }
                recent_relative_changes.push(change);
                if post_burn_history.len() == 50 {
                    post_burn_history.remove(0);
                }
                post_burn_history.push([
                    next_log_theta[0],
                    next_log_theta[1],
                    libm::log(next_omega[0][0]),
                    libm::log(next_omega[1][1]),
                    next_omega[0][1] / libm::sqrt(next_omega[0][0] * next_omega[1][1]),
                    residual
                        .as_ref()
                        .map(|fit| libm::log(fit.variance))
                        .unwrap_or(0.),
                ]);
            }
            theta = next_theta;
            omega = next_omega;
            log_theta = next_log_theta;
            iterations = iteration + 1;
            let window_stable = if post_burn_history.len() == 50 && moved > 0 {
                (0..6).all(|axis| {
                    let history: Vec<f64> = post_burn_history.iter().map(|x| x[axis]).collect();
                    let separated: Vec<f64> = history[..15]
                        .iter()
                        .chain(&history[35..])
                        .copied()
                        .collect();
                    windowed_drift(&history[25..], request.gradient_tolerance).2
                        && windowed_drift(&separated, request.gradient_tolerance).2
                })
            } else {
                false
            };
            stable_iterations = if window_stable {
                stable_iterations + 1
            } else {
                0
            };
            if stable_iterations >= 4 {
                status = FitStatus::Converged;
                break;
            }
        }
        let sampled_etas: Vec<Vector> = phi
            .iter()
            .map(|p| [p[0] - log_theta[0], p[1] - log_theta[1]])
            .collect();
        let final_request = residual.as_ref().map(|_| {
            let mut updated = request.clone();
            updated.subjects = subjects.clone();
            updated
        });
        let fitted_request = final_request.as_ref().unwrap_or(request);
        let marginal = self.saem_marginal_ofv(
            fitted_request,
            &setup,
            theta,
            omega,
            &sampled_etas,
            &mut evaluations,
            request.max_evaluations - 3 * count,
        )?;
        if marginal.is_none() && status == FitStatus::Converged {
            status = if evaluations >= request.max_evaluations - 3 * count {
                FitStatus::EvaluationLimit
            } else {
                FitStatus::LineSearchFailed
            };
        }
        let final_etas = marginal
            .as_ref()
            .map(|(_, modes)| modes)
            .unwrap_or(&sampled_etas);
        let mut subject_results = Vec::with_capacity(count);
        let mut fitted_observations = Vec::new();
        let mut surrogate = 0.;
        for (index, subject) in fitted_request.subjects.iter().enumerate() {
            let eta = final_etas[index];
            let likelihood = self.saem_nll(subject, &setup, theta, eta)?;
            evaluations += 1;
            surrogate += likelihood + prior_nll_vector(eta, omega, 2)?;
            fitted_observations.extend(self.saem_predictions(index, subject, &setup, theta, eta)?);
            evaluations += 2;
            subject_results.push(PopulationSubjectEstimate {
                eta: eta.to_vec(),
                negative_log_likelihood: likelihood,
            });
        }
        if !surrogate.is_finite() {
            return Err(Error::new(ErrorCode::Domain, "SAEM surrogate overflow"));
        }
        let mut fixed_effects = std::collections::BTreeMap::new();
        for axis in 0..2 {
            fixed_effects.insert(
                setup.names[axis].clone(),
                pharmflux_core::model::Quantity {
                    value: theta[axis],
                    unit: setup.units[axis].clone(),
                },
            );
        }
        let recent_window = recent_relative_changes.len();
        recent_relative_changes.sort_by(f64::total_cmp);
        let recent_relative_change_median = if recent_window == 0 {
            None
        } else if recent_window % 2 == 0 {
            Some(
                0.5 * (recent_relative_changes[recent_window / 2 - 1]
                    + recent_relative_changes[recent_window / 2]),
            )
        } else {
            Some(recent_relative_changes[recent_window / 2])
        };
        Ok(PopulationFitResult {
            status,
            fixed_effects,
            omega: vec![omega[0].to_vec(), omega[1].to_vec()],
            subjects: subject_results,
            fitted_observations: Some(fitted_observations),
            uncertainty: None,
            estimated_residual_error: residual.map(|fit| PopulationResidualErrorEstimate {
                output: fit.output,
                additive_sd: pharmflux_core::model::Quantity {
                    value: libm::sqrt(fit.variance),
                    unit: fit.unit,
                },
            }),
            saem_diagnostics: Some(PopulationSaemDiagnostics {
                burn_in_iterations: burn_in,
                final_step_size,
                final_relative_change,
                recent_relative_change_median,
                recent_window,
                evaluated_proposals,
                accepted_proposals,
                consecutive_stable_iterations: stable_iterations,
            }),
            objective_kind: if marginal.is_some() {
                PopulationObjectiveKind::SaemMarginalObjectiveFunctionValue
            } else {
                PopulationObjectiveKind::SaemCompleteDataSurrogate
            },
            objective: marginal.map(|(value, _)| value).unwrap_or(surrogate),
            evaluations,
            iterations,
            fit_request_hash: canonical_hash(&serde_json::json!(request))?,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::{windowed_drift, Stream};

    #[test]
    fn stream_is_reproducible_and_finite() {
        let mut a = Stream::new(42);
        let mut b = Stream::new(42);
        for _ in 0..1000 {
            let x = a.normal();
            assert_eq!(x, b.normal());
            assert!(x.is_finite());
        }
    }

    #[test]
    fn windowed_drift_rejects_trend_and_high_noise() {
        let flat: Vec<f64> = (0..50).map(|i| 1.0 + (i % 2) as f64 * 0.0001).collect();
        assert!(windowed_drift(&flat, 1e-5).2);
        let trending: Vec<f64> = (0..50).map(|i| 1.0 + i as f64 * 0.001).collect();
        assert!(!windowed_drift(&trending, 1e-5).2);
        assert!(!windowed_drift(&trending, 0.1).2);
        let noisy: Vec<f64> = (0..50)
            .map(|i| 1.0 + (i % 2) as f64 * 0.1 + i as f64 * 0.001)
            .collect();
        assert!(!windowed_drift(&noisy, 1e-5).2);
        let correlated_drift: Vec<f64> = (0..50)
            .map(|i| 1.0 + 0.0008 * i as f64 + 0.002 * libm::sin(i as f64 / 6.))
            .collect();
        assert!(!windowed_drift(&correlated_drift, 1e-5).2);
    }
}
