//! Bounded one-random-effect Laplace approximation to marginal likelihood.
//! The conditional Gaussian likelihood includes prediction-dependent variance.
use super::*;
use pharmflux_core::{fit::*, identity::canonical_hash};

struct Eval {
    value: f64,
    subjects: Vec<PopulationSubjectEstimate>,
}

struct Evaluator<'a> {
    model: &'a Arc<CompiledSensitivityDocument>,
    request: &'a PopulationFitRequest,
    parameter: &'a str,
    unit: &'a str,
    calls: usize,
    warm_etas: Vec<f64>,
}

impl Evaluator<'_> {
    fn subject(
        &mut self,
        subject: &GaussianObjectiveRequest,
        theta: f64,
        eta: f64,
    ) -> Result<(f64, f64), Error> {
        if self.calls >= self.request.max_evaluations {
            return Err(Error::new(
                ErrorCode::WorkBudget,
                "Laplace evaluation limit reached",
            ));
        }
        self.calls += 1;
        let result = self.model.population_subject_likelihood(
            subject,
            self.parameter,
            theta,
            self.unit,
            eta,
        )?;
        let score = result.gradient[0].value * theta * libm::exp(eta);
        if !score.is_finite() {
            return Err(Error::new(ErrorCode::Domain, "Laplace score overflow"));
        }
        Ok((result.negative_log_likelihood, score))
    }

    fn conditional(
        &mut self,
        subject: &GaussianObjectiveRequest,
        theta: f64,
        omega: f64,
        warm_eta: f64,
    ) -> Result<(f64, f64, f64), Error> {
        let limit = self.request.eta_bound;
        let score = |this: &mut Self, eta: f64| -> Result<(f64, f64), Error> {
            let (nll, derivative) = this.subject(subject, theta, eta)?;
            Ok((nll, derivative + eta / omega))
        };
        // Sample the full declared eta domain to detect extra score roots. The
        // samples also furnish a tight bracket for the safeguarded secant solve.
        let mut grid = Vec::with_capacity(17);
        for point in 0..=16 {
            let eta = -limit + limit * point as f64 / 8.;
            let (nll, g) = score(self, eta)?;
            grid.push((eta, nll, g));
        }
        if grid[0].2 >= 0. || grid[16].2 <= 0. {
            return Err(Error::new(
                ErrorCode::Unsupported,
                "Laplace eta mode reaches its bound",
            ));
        }
        let brackets = grid
            .windows(2)
            .filter(|pair| pair[0].2 <= 0. && pair[1].2 > 0.)
            .collect::<Vec<_>>();
        if brackets.len() != 1
            || grid
                .windows(2)
                .any(|pair| pair[0].2 > 0. && pair[1].2 <= 0.)
        {
            return Err(Error::new(
                ErrorCode::Unsupported,
                "Laplace conditional score has multiple modes in the eta domain",
            ));
        }
        let mut left = brackets[0][0];
        let mut right = brackets[0][1];
        let mut best = if left.2.abs() < right.2.abs() {
            left
        } else {
            right
        };
        if warm_eta > left.0 && warm_eta < right.0 {
            let (nll, g) = score(self, warm_eta)?;
            let trial = (warm_eta, nll, g);
            if g > 0. {
                right = trial;
            } else {
                left = trial;
            }
            if g.abs() < best.2.abs() {
                best = trial;
            }
        }
        for _ in 0..16 {
            if best.2.abs() <= 1e-4 || right.0 - left.0 <= 1e-6 {
                break;
            }
            let width = right.0 - left.0;
            let eta = (left.0 - left.2 * width / (right.2 - left.2))
                .clamp(left.0 + 0.1 * width, right.0 - 0.1 * width);
            let (nll, g) = score(self, eta)?;
            let trial = (eta, nll, g);
            if g > 0. {
                right = trial;
            } else {
                left = trial;
            }
            if g.abs() < best.2.abs() {
                best = trial;
            }
        }
        // Flat scores can slow secant interpolation. Finish those subjects
        // with bisection while retaining the same root bracket.
        if best.2.abs() > 1e-4 {
            for _ in 0..18 {
                if best.2.abs() <= 1e-4 || right.0 - left.0 <= 1e-6 {
                    break;
                }
                let eta = (left.0 + right.0) / 2.;
                let (nll, g) = score(self, eta)?;
                let trial = (eta, nll, g);
                if g > 0. {
                    right = trial;
                } else {
                    left = trial;
                }
                if g.abs() < best.2.abs() {
                    best = trial;
                }
            }
        }
        if best.2.abs() > 1e-2 {
            return Err(Error::new(
                ErrorCode::Unsupported,
                "Laplace conditional mode did not converge",
            ));
        }
        let (eta, nll, _) = best;
        let h = 1e-3_f64.min((limit - eta.abs()) / 4.);
        if h <= 0. {
            return Err(Error::new(
                ErrorCode::Unsupported,
                "Laplace eta mode has no curvature domain",
            ));
        }
        let curvature = (score(self, eta + h)?.1 - score(self, eta - h)?.1) / (2. * h);
        if !curvature.is_finite() || curvature <= 0. {
            return Err(Error::new(
                ErrorCode::Unsupported,
                "Laplace conditional mode lacks positive curvature",
            ));
        }
        let laplace = nll + 0.5 * eta * eta / omega + 0.5 * libm::log(omega * curvature);
        Ok((eta, nll, laplace))
    }

    fn population(&mut self, theta: f64, omega: f64) -> Result<Eval, Error> {
        let mut value = 0.;
        let mut subjects = Vec::with_capacity(self.request.subjects.len());
        for (i, subject) in self.request.subjects.iter().enumerate() {
            let warm_eta = self.warm_etas[i];
            let (eta, nll, term) =
                self.conditional(subject, theta, omega, warm_eta)
                    .map_err(|mut e| {
                        e.expression = Some(format!("subjects[{i}]"));
                        e
                    })?;
            self.warm_etas[i] = eta;
            value += term;
            subjects.push(PopulationSubjectEstimate {
                eta: vec![eta],
                negative_log_likelihood: nll,
            });
        }
        if !value.is_finite() {
            return Err(Error::new(ErrorCode::Domain, "Laplace objective overflow"));
        }
        Ok(Eval { value, subjects })
    }
}

impl CompiledSensitivityDocument {
    /// Estimate the population mean and, when its bounds differ, eta variance.
    /// The returned objective is a Laplace approximation to marginal NLL.
    /// `Converged` means the normalized coordinate-search step reached the
    /// request tolerance; it does not certify a projected-gradient threshold.
    pub fn fit_population_laplace(
        self: &Arc<Self>,
        request: &PopulationFitRequest,
    ) -> Result<PopulationFitResult, Error> {
        if request.uncertainty {
            return Err(Error::new(
                ErrorCode::Unsupported,
                "Laplace uncertainty is not implemented",
            ));
        }
        if request.additive_error_fit.is_some() {
            return Err(Error::new(
                ErrorCode::Unsupported,
                "Laplace additive-error estimation is not implemented",
            ));
        }
        let setup = self.validate_population_1d(request)?;
        let mut evaluator = Evaluator {
            model: self,
            request,
            parameter: &setup.parameter,
            unit: &setup.unit,
            calls: 0,
            warm_etas: vec![0.; request.subjects.len()],
        };
        let omega_free = setup.omega_lower < setup.omega_upper;
        let mut x = [
            (setup.start - setup.lower) / (setup.upper - setup.lower),
            if omega_free {
                (setup.omega - setup.omega_lower) / (setup.omega_upper - setup.omega_lower)
            } else {
                0.
            },
        ];
        let decode = |x: &[f64; 2]| {
            (
                setup.lower + x[0] * (setup.upper - setup.lower),
                if omega_free {
                    setup.omega_lower + x[1] * (setup.omega_upper - setup.omega_lower)
                } else {
                    setup.omega
                },
            )
        };
        let (theta, omega) = decode(&x);
        let mut best = evaluator.population(theta, omega)?;
        let mut steps: [f64; 2] = [0.2, if omega_free { 0.2 } else { 0. }];
        let mut iterations = 0;
        let status = loop {
            if steps[0].max(steps[1]) <= request.gradient_tolerance {
                break FitStatus::Converged;
            }
            if iterations >= request.max_iterations {
                break FitStatus::IterationLimit;
            }
            if evaluator.calls >= request.max_evaluations {
                break FitStatus::EvaluationLimit;
            }
            let mut improved = false;
            for dimension in 0..(if omega_free { 2 } else { 1 }) {
                for direction in [-1., 1.] {
                    let mut trial = x;
                    trial[dimension] =
                        (trial[dimension] + direction * steps[dimension]).clamp(0., 1.);
                    if trial == x {
                        continue;
                    }
                    let (theta, omega) = decode(&trial);
                    match evaluator.population(theta, omega) {
                        Ok(candidate) if candidate.value < best.value => {
                            x = trial;
                            best = candidate;
                            improved = true;
                        }
                        Ok(_) => {}
                        Err(e) if e.code == ErrorCode::WorkBudget => break,
                        Err(e)
                            if matches!(
                                e.code,
                                ErrorCode::Domain
                                    | ErrorCode::Solver
                                    | ErrorCode::Invariant
                                    | ErrorCode::Unsupported
                            ) => {}
                        Err(e) => return Err(e),
                    }
                    if evaluator.calls >= request.max_evaluations {
                        break;
                    }
                }
                if evaluator.calls >= request.max_evaluations {
                    break;
                }
            }
            iterations += 1;
            if !improved {
                for step in &mut steps {
                    *step *= 0.5;
                }
            }
        };
        let (theta, omega) = decode(&x);
        Ok(PopulationFitResult {
            status,
            fixed_effects: BTreeMap::from([(
                setup.parameter.clone(),
                pharmflux_core::model::Quantity {
                    value: theta,
                    unit: setup.unit.clone(),
                },
            )]),
            omega: vec![vec![omega]],
            subjects: best.subjects,
            fitted_observations: None,
            uncertainty: None,
            estimated_residual_error: None,
            saem_diagnostics: None,
            objective_kind: PopulationObjectiveKind::LaplaceNegativeLogLikelihood,
            objective: best.value,
            evaluations: evaluator.calls,
            iterations,
            fit_request_hash: canonical_hash(
                &serde_json::json!({"algorithm":"laplace_1d_v1","request":request}),
            )?,
        })
    }
}
