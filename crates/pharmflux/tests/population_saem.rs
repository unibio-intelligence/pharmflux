use pharmflux::document::sensitivity::CompiledSensitivityDocument;
use pharmflux_core::{
    expression::{Binary, Expr, Function},
    fit::*,
    model::{Initial, ModelDocument, Output, Parameter, Quantity, Transform},
    run::SensitivityRequest,
};

fn quantity(value: f64, unit: &str) -> Quantity {
    Quantity {
        value,
        unit: unit.into(),
    }
}

fn fixture() -> (
    std::sync::Arc<CompiledSensitivityDocument>,
    PopulationFitRequest,
) {
    let model: ModelDocument = serde_json::from_str(include_str!(
        "../../../conformance/models/synthetic-one-compartment.json"
    ))
    .unwrap();
    let mut simulation: SensitivityRequest = serde_json::from_str(include_str!(
        "../../../conformance/requests/synthetic-sensitivity.json"
    ))
    .unwrap();
    simulation.with_respect_to = vec!["cl".into()];
    simulation.derivative_absolute_tolerances.remove("v");
    simulation
        .run
        .parameters
        .insert("cl".into(), quantity(0.5, "L/h"));
    simulation
        .run
        .parameters
        .insert("v".into(), quantity(4.5, "L"));
    let compiled = CompiledSensitivityDocument::compile(&model, &["cl".into()]).unwrap();
    let eta: [f64; 8] = [-0.45, -0.32, -0.21, -0.12, 0.12, 0.21, 0.32, 0.45];
    let times = compiled.run(&simulation).unwrap().result.times;
    let subjects = eta
        .iter()
        .map(|eta| GaussianObjectiveRequest {
            simulation: simulation.clone(),
            observations: times
                .iter()
                .enumerate()
                .filter(|(row, _)| *row > 0)
                .map(|(row, time)| GaussianObservation {
                    row,
                    output: "cp".into(),
                    value: quantity(6.0 / 4.5 * (-0.6 * (*eta).exp() * time / 4.5).exp(), "mg/L"),
                    error: GaussianError {
                        additive_sd: quantity(0.02, "mg/L"),
                        proportional_sd: 0.0,
                    },
                })
                .collect(),
            priors: Vec::new(),
        })
        .collect();
    let request = PopulationFitRequest {
        subjects,
        fixed_effects: vec![FitBound {
            parameter: "cl".into(),
            lower: quantity(0.2, "L/h"),
            upper: quantity(1.0, "L/h"),
        }],
        random_effects: vec![PopulationRandomEffect {
            parameter: "cl".into(),
        }],
        omega: vec![vec![0.1]],
        omega_diagonal_bounds: vec![[0.01, 0.5]],
        eta_bound: 2.0,
        max_evaluations: 20_000,
        max_iterations: 800,
        gradient_tolerance: 0.001,
        seed: 20260929,
        burn_in_iterations: None,
        uncertainty: false,
        additive_error_fit: None,
    };
    (compiled, request)
}

fn two_effect_fixture() -> (
    std::sync::Arc<CompiledSensitivityDocument>,
    PopulationFitRequest,
) {
    let mut model: ModelDocument = serde_json::from_str(include_str!(
        "../../../conformance/models/synthetic-one-compartment.json"
    ))
    .unwrap();
    model.model_id = "saem_two_log_effects".into();
    model.parameters[0].default = quantity(1., "1");
    model.parameters[1].default = quantity(1., "1");
    model.parameters.push(Parameter {
        name: "k".into(),
        default: quantity(1., "1/h"),
        bounds: None,
        transform: Transform::Identity,
        fixed: true,
        description: String::new(),
    });
    model.states[0].unit = "1".into();
    model.states[0].initial = Initial::Quantity {
        quantity: quantity(0., "1"),
    };
    model.states[0].rhs = Expr::binary(
        Binary::Multiply,
        Expr::number(-1.),
        Expr::binary(Binary::Multiply, Expr::symbol("k"), Expr::symbol("central")),
    );
    model.states[0].dosing = None;
    model.outputs[0].unit = "1".into();
    model.outputs[0].expression = Expr::call(Function::Log, vec![Expr::symbol("cl")]);
    model.outputs.push(Output {
        name: "log_v".into(),
        unit: "1".into(),
        expression: Expr::call(Function::Log, vec![Expr::symbol("v")]),
    });
    let compiled =
        CompiledSensitivityDocument::compile(&model, &["cl".into(), "v".into()]).unwrap();
    let mut simulation: SensitivityRequest = serde_json::from_str(include_str!(
        "../../../conformance/requests/synthetic-sensitivity.json"
    ))
    .unwrap();
    simulation.with_respect_to = vec!["cl".into(), "v".into()];
    simulation.derivative_absolute_tolerances.clear();
    for name in ["cl", "v"] {
        simulation.derivative_absolute_tolerances.insert(
            name.into(),
            std::collections::BTreeMap::from([("central".into(), quantity(1e-12, "1"))]),
        );
        simulation
            .run
            .parameters
            .insert(name.into(), quantity(1., "1"));
    }
    simulation.run.regimen.administrations.clear();
    let effects = [
        [-0.24, -0.19],
        [-0.18, -0.13],
        [-0.11, -0.05],
        [-0.03, 0.02],
        [0.03, -0.02],
        [0.11, 0.05],
        [0.18, 0.13],
        [0.24, 0.19],
    ];
    let subjects = effects
        .iter()
        .map(|eta| GaussianObjectiveRequest {
            simulation: simulation.clone(),
            observations: [
                GaussianObservation {
                    row: 0,
                    output: "cp".into(),
                    value: quantity(0.9_f64.ln() + eta[0], "1"),
                    error: GaussianError {
                        additive_sd: quantity(0.06, "1"),
                        proportional_sd: 0.,
                    },
                },
                GaussianObservation {
                    row: 0,
                    output: "log_v".into(),
                    value: quantity(1.1_f64.ln() + eta[1], "1"),
                    error: GaussianError {
                        additive_sd: quantity(0.06, "1"),
                        proportional_sd: 0.,
                    },
                },
            ]
            .into(),
            priors: vec![],
        })
        .collect();
    (
        compiled,
        PopulationFitRequest {
            subjects,
            fixed_effects: vec![
                FitBound {
                    parameter: "cl".into(),
                    lower: quantity(0.7, "1"),
                    upper: quantity(1.3, "1"),
                },
                FitBound {
                    parameter: "v".into(),
                    lower: quantity(0.7, "1"),
                    upper: quantity(1.3, "1"),
                },
            ],
            random_effects: vec![
                PopulationRandomEffect {
                    parameter: "cl".into(),
                },
                PopulationRandomEffect {
                    parameter: "v".into(),
                },
            ],
            omega: vec![vec![0.04, 0.012], vec![0.012, 0.04]],
            omega_diagonal_bounds: vec![[0.005, 0.12], [0.005, 0.12]],
            eta_bound: 1.5,
            max_evaluations: 40_000,
            max_iterations: 800,
            gradient_tolerance: 0.001,
            seed: 20260930,
            burn_in_iterations: None,
            uncertainty: false,
            additive_error_fit: None,
        },
    )
}

fn two_effect_dosed_pk_fixture() -> (
    std::sync::Arc<CompiledSensitivityDocument>,
    PopulationFitRequest,
) {
    let model: ModelDocument = serde_json::from_str(include_str!(
        "../../../conformance/models/synthetic-one-compartment.json"
    ))
    .unwrap();
    let mut simulation: SensitivityRequest = serde_json::from_str(include_str!(
        "../../../conformance/requests/synthetic-sensitivity.json"
    ))
    .unwrap();
    simulation
        .run
        .parameters
        .insert("cl".into(), quantity(0.5, "L/h"));
    simulation
        .run
        .parameters
        .insert("v".into(), quantity(4.0, "L"));
    let compiled =
        CompiledSensitivityDocument::compile(&model, &["cl".into(), "v".into()]).unwrap();
    let grid = compiled.run(&simulation).unwrap().result;
    let effects = [
        [-0.35, -0.22],
        [-0.28, 0.12],
        [-0.21, -0.08],
        [-0.14, 0.25],
        [-0.07, -0.18],
        [0.07, 0.18],
        [0.14, -0.25],
        [0.21, 0.08],
        [0.28, -0.12],
        [0.35, 0.22],
    ];
    let subjects = effects
        .iter()
        .map(|eta| {
            let cl = 0.6 * libm::exp(eta[0]);
            let v = 4.5 * libm::exp(eta[1]);
            GaussianObjectiveRequest {
                simulation: simulation.clone(),
                observations: grid
                    .times
                    .iter()
                    .enumerate()
                    .filter(|(row, time)| *row > 0 && **time > 0.)
                    .map(|(row, time)| GaussianObservation {
                        row,
                        output: "cp".into(),
                        value: quantity(6. / v * libm::exp(-cl * time / v), "mg/L"),
                        error: GaussianError {
                            additive_sd: quantity(0.02, "mg/L"),
                            proportional_sd: 0.,
                        },
                    })
                    .collect(),
                priors: Vec::new(),
            }
        })
        .collect();
    (
        compiled,
        PopulationFitRequest {
            subjects,
            fixed_effects: vec![
                FitBound {
                    parameter: "cl".into(),
                    lower: quantity(0.3, "L/h"),
                    upper: quantity(1., "L/h"),
                },
                FitBound {
                    parameter: "v".into(),
                    lower: quantity(2.5, "L"),
                    upper: quantity(7., "L"),
                },
            ],
            random_effects: vec![
                PopulationRandomEffect {
                    parameter: "cl".into(),
                },
                PopulationRandomEffect {
                    parameter: "v".into(),
                },
            ],
            omega: vec![vec![0.05, 0.], vec![0., 0.05]],
            omega_diagonal_bounds: vec![[0.005, 0.2], [0.005, 0.2]],
            eta_bound: 1.5,
            max_evaluations: 100_000,
            max_iterations: 1200,
            gradient_tolerance: 0.001,
            seed: 20260930,
            burn_in_iterations: None,
            uncertainty: false,
            additive_error_fit: None,
        },
    )
}

#[test]
fn saem_two_effects_recovers_dosed_pk_and_publishes_fitted_rows() {
    let (compiled, request) = two_effect_dosed_pk_fixture();
    let fitted = compiled.fit_population_saem(&request).unwrap();
    assert_eq!(fitted.status, FitStatus::Converged, "{fitted:?}");
    assert_eq!(
        fitted.objective_kind,
        PopulationObjectiveKind::SaemMarginalObjectiveFunctionValue
    );
    assert!((fitted.fixed_effects["cl"].value - 0.6).abs() < 0.08);
    assert!((fitted.fixed_effects["v"].value - 4.5).abs() < 0.45);
    assert!(fitted.objective.is_finite());
    let rows = fitted.fitted_observations.unwrap();
    assert_eq!(rows.len(), 30);
    for row in rows {
        let cl = fitted.fixed_effects["cl"].value;
        let v = fitted.fixed_effects["v"].value;
        let eta = &fitted.subjects[row.subject_index].eta;
        let pred = 6. / v * libm::exp(-cl * row.time.value / v);
        let individual_cl = cl * libm::exp(eta[0]);
        let individual_v = v * libm::exp(eta[1]);
        let ipred = 6. / individual_v * libm::exp(-individual_cl * row.time.value / individual_v);
        assert!((row.pred.value - pred).abs() < 1e-6, "{row:?}");
        assert!((row.ipred.value - ipred).abs() < 1e-6, "{row:?}");
    }
}

#[test]
fn saem_two_effects_dosed_pk_converges_across_seeds() {
    let (compiled, mut request) = two_effect_dosed_pk_fixture();
    for seed in [17029, 32321] {
        request.seed = seed;
        let fitted = compiled.fit_population_saem(&request).unwrap();
        assert_eq!(
            fitted.status,
            FitStatus::Converged,
            "seed={seed}, {fitted:?}"
        );
        assert_eq!(
            fitted.objective_kind,
            PopulationObjectiveKind::SaemMarginalObjectiveFunctionValue,
            "seed={seed}"
        );
        assert!(
            (fitted.fixed_effects["cl"].value - 0.6).abs() < 0.08,
            "seed={seed}"
        );
        assert!(
            (fitted.fixed_effects["v"].value - 4.5).abs() < 0.45,
            "seed={seed}"
        );
        assert_eq!(fitted.fitted_observations.unwrap().len(), 30);
    }
}

#[test]
fn saem_one_effect_estimates_additive_error_and_preserves_fitted_rows() {
    let (compiled, mut request) = fixture();
    for (subject_index, subject) in request.subjects.iter_mut().enumerate() {
        for (row, observation) in subject.observations.iter_mut().enumerate() {
            observation.value.value += if (subject_index + row) % 2 == 0 {
                0.015
            } else {
                -0.015
            };
            observation.error.additive_sd = quantity(0.03, "mg/L");
        }
    }
    request.additive_error_fit = Some(PopulationAdditiveErrorFit {
        output: "cp".into(),
        initial: quantity(0.03, "mg/L"),
        lower: quantity(0.005, "mg/L"),
        upper: quantity(0.05, "mg/L"),
    });
    request.max_iterations = 1200;
    request.max_evaluations = 100_000;
    for seed in [20260929, 17029, 32321] {
        request.seed = seed;
        let fitted = compiled.fit_population_saem(&request).unwrap();
        assert_eq!(
            fitted.status,
            FitStatus::Converged,
            "seed={seed}, {fitted:?}"
        );
        assert_eq!(
            fitted.objective_kind,
            PopulationObjectiveKind::SaemMarginalObjectiveFunctionValue
        );
        let error = fitted.estimated_residual_error.as_ref().unwrap();
        assert_eq!(error.output, "cp");
        assert_eq!(error.additive_sd.unit, "mg/L");
        // Synthetic alternating residuals imply an interior residual scale.
        assert!(
            (error.additive_sd.value - 0.015).abs() < 0.005,
            "seed={seed}, {fitted:?}"
        );
        assert!((fitted.fixed_effects["cl"].value - 0.6).abs() < 0.02);
        assert!((fitted.omega[0][0] - 0.088).abs() < 0.02);
        assert_eq!(
            fitted.fitted_observations.as_ref().unwrap().len(),
            request
                .subjects
                .iter()
                .map(|s| s.observations.len())
                .sum::<usize>()
        );
    }
}

#[test]
fn saem_burn_in_is_explicit_and_bounded() {
    let (compiled, mut request) = fixture();
    request.burn_in_iterations = Some(200);
    let fitted = compiled.fit_population_saem(&request).unwrap();
    assert_eq!(fitted.saem_diagnostics.unwrap().burn_in_iterations, 200);
    request.burn_in_iterations = Some(request.max_iterations);
    assert_eq!(
        compiled.fit_population_saem(&request).unwrap_err().code,
        pharmflux_core::ErrorCode::InvalidInput
    );
}

#[test]
fn saem_two_effects_free_error_preserves_evaluation_budget() {
    let (compiled, mut request) = two_effect_dosed_pk_fixture();
    for subject in &mut request.subjects {
        for observation in &mut subject.observations {
            observation.error.additive_sd = quantity(0.04, "mg/L");
        }
    }
    request.additive_error_fit = Some(PopulationAdditiveErrorFit {
        output: "cp".into(),
        initial: quantity(0.04, "mg/L"),
        lower: quantity(0.005, "mg/L"),
        upper: quantity(0.1, "mg/L"),
    });
    request.max_evaluations = request.subjects.len() * 4;
    let fit = compiled.fit_population_saem(&request).unwrap();
    assert_eq!(fit.status, FitStatus::EvaluationLimit);
    assert!(fit.evaluations <= request.max_evaluations);
    assert_eq!(fit.fitted_observations.as_ref().unwrap().len(), 30);
}

#[test]
fn saem_two_effects_matches_exact_gaussian_marginal_at_fitted_parameters() {
    let (compiled, request) = two_effect_fixture();
    let fitted = compiled.fit_population_saem(&request).unwrap();
    assert_eq!(fitted.status, FitStatus::Converged, "{fitted:?}");
    assert_eq!(
        fitted.objective_kind,
        PopulationObjectiveKind::SaemMarginalObjectiveFunctionValue,
        "{fitted:?}"
    );
    assert!((fitted.fixed_effects["cl"].value - 0.9).abs() < 0.12);
    assert!((fitted.fixed_effects["v"].value - 1.1).abs() < 0.12);
    assert_eq!(fitted.subjects.len(), request.subjects.len());
    assert_eq!(fitted.fitted_observations.as_ref().unwrap().len(), 16);
    let covariance = [
        [fitted.omega[0][0] + 0.06_f64.powi(2), fitted.omega[0][1]],
        [fitted.omega[1][0], fitted.omega[1][1] + 0.06_f64.powi(2)],
    ];
    let determinant = covariance[0][0] * covariance[1][1] - covariance[0][1].powi(2);
    assert!(determinant > 0.);
    let mut exact = 0.;
    for subject in &request.subjects {
        let residual = [
            subject.observations[0].value.value - fitted.fixed_effects["cl"].value.ln(),
            subject.observations[1].value.value - fitted.fixed_effects["v"].value.ln(),
        ];
        exact += (covariance[1][1] * residual[0].powi(2) + covariance[0][0] * residual[1].powi(2)
            - 2. * covariance[0][1] * residual[0] * residual[1])
            / determinant
            + determinant.ln();
    }
    assert!(
        (fitted.objective - exact).abs() < 0.02,
        "Laplace OFV {} vs exact {}",
        fitted.objective,
        exact
    );
    let grid = compiled
        .run(&request.subjects[0].simulation)
        .unwrap()
        .result;
    for row in fitted.fitted_observations.unwrap() {
        let axis = row.observation_index;
        let theta = if axis == 0 {
            fitted.fixed_effects["cl"].value
        } else {
            fitted.fixed_effects["v"].value
        };
        assert!((row.pred.value - theta.ln()).abs() < 1e-10);
        assert!(
            (row.ipred.value - theta.ln() - fitted.subjects[row.subject_index].eta[axis]).abs()
                < 1e-10
        );
        assert_eq!(row.side, grid.sides[row.row]);
        assert_eq!(row.time.value, grid.times[row.row]);
    }
}

#[test]
fn tiny_positive_residual_variance_has_continuous_likelihood_and_score() {
    let (compiled, request) = two_effect_fixture();
    let mut subject = request.subjects[0].clone();
    subject.observations.truncate(1);
    subject.observations[0].value.value = 0.;
    subject.observations[0].error.additive_sd.value = 1e-5;
    let narrow = compiled.gaussian_objective(&subject).unwrap();
    assert!(narrow.objective.is_finite());
    subject.observations[0].error.additive_sd.value = 2e-5;
    let wider = compiled.gaussian_objective(&subject).unwrap();
    assert!((wider.objective - narrow.objective - 2_f64.ln()).abs() < 1e-9);

    subject.observations[0].error.additive_sd.value = 1e-5;
    subject.observations[0].value.value = 1e-6;
    let center = compiled.gaussian_objective(&subject).unwrap();
    let step = 1e-7;
    subject
        .simulation
        .run
        .parameters
        .get_mut("cl")
        .unwrap()
        .value = 1. + step;
    let plus = compiled.gaussian_objective(&subject).unwrap().objective;
    subject
        .simulation
        .run
        .parameters
        .get_mut("cl")
        .unwrap()
        .value = 1. - step;
    let minus = compiled.gaussian_objective(&subject).unwrap().objective;
    let numeric = (plus - minus) / (2. * step);
    assert!(
        (center.gradient[0].value - numeric).abs() < 0.01,
        "analytic {} vs finite difference {numeric}",
        center.gradient[0].value
    );
}

#[test]
fn saem_recovers_lognormal_population_on_synthetic_pk() {
    let (compiled, request) = fixture();
    let result = compiled.fit_population_saem(&request).unwrap();
    assert_eq!(result.status, FitStatus::Converged, "{result:?}");
    assert_eq!(
        result.objective_kind,
        PopulationObjectiveKind::SaemMarginalObjectiveFunctionValue
    );
    assert_eq!(result.subjects.len(), 8);
    assert!(result.evaluations <= request.max_evaluations);
    assert!(result.objective.is_finite());
    let rows = result.fitted_observations.as_ref().unwrap();
    assert_eq!(
        rows.len(),
        request
            .subjects
            .iter()
            .map(|s| s.observations.len())
            .sum::<usize>()
    );
    for row in rows {
        assert_eq!(row.value.unit, row.pred.unit);
        assert_eq!(row.pred.unit, row.ipred.unit);
        assert_eq!(
            row.output,
            request.subjects[row.subject_index].observations[row.observation_index].output
        );
        assert_eq!(
            row.row,
            request.subjects[row.subject_index].observations[row.observation_index].row
        );
    }
    let theta = result.fixed_effects["cl"].value;
    let omega = result.omega[0][0];
    assert!((theta - 0.6).abs() < 0.12, "theta={theta}");
    assert!((omega - 0.083).abs() < 0.06, "omega={omega}");
    let diagnostics = result.saem_diagnostics.unwrap();
    assert_eq!(diagnostics.burn_in_iterations, request.max_iterations / 3);
    assert!(diagnostics.accepted_proposals <= diagnostics.evaluated_proposals);
    assert!(diagnostics.evaluated_proposals <= result.evaluations - 4 * request.subjects.len());
    if result.iterations > 0 {
        let last = result.iterations - 1;
        let expected_gain = if last < diagnostics.burn_in_iterations {
            1.0
        } else {
            1.0 / (last - diagnostics.burn_in_iterations + 1) as f64
        };
        assert!((diagnostics.final_step_size.unwrap() - expected_gain).abs() < 1e-12);
        assert!(diagnostics.final_relative_change.unwrap().is_finite());
        assert!(diagnostics.recent_window <= 20);
        if diagnostics.recent_window > 0 {
            assert!(diagnostics
                .recent_relative_change_median
                .unwrap()
                .is_finite());
        }
    }
}

#[test]
fn saem_converges_on_twenty_four_subject_pk_cohort() {
    let (compiled, mut request) = fixture();
    let simulation = request.subjects[0].simulation.clone();
    let times = compiled.run(&simulation).unwrap().result.times;
    request.subjects = (0..24)
        .map(|index| {
            let eta = -0.45 + 0.9 * index as f64 / 23.;
            GaussianObjectiveRequest {
                simulation: simulation.clone(),
                observations: times
                    .iter()
                    .enumerate()
                    .filter(|(row, _)| *row > 0)
                    .map(|(row, time)| GaussianObservation {
                        row,
                        output: "cp".into(),
                        value: quantity(6.0 / 4.5 * (-0.6 * eta.exp() * time / 4.5).exp(), "mg/L"),
                        error: GaussianError {
                            additive_sd: quantity(0.02, "mg/L"),
                            proportional_sd: 0.,
                        },
                    })
                    .collect(),
                priors: Vec::new(),
            }
        })
        .collect();
    request.max_iterations = 600;
    request.max_evaluations = 60_000;
    let result = compiled.fit_population_saem(&request).unwrap();
    assert_eq!(
        result.status,
        FitStatus::Converged,
        "{:?}",
        result.saem_diagnostics
    );
    assert_eq!(
        result.objective_kind,
        PopulationObjectiveKind::SaemMarginalObjectiveFunctionValue
    );
    assert!((result.fixed_effects["cl"].value - 0.6).abs() < 0.05);
    assert!((result.omega[0][0] - 0.073).abs() < 0.035);
    assert_eq!(result.fitted_observations.unwrap().len(), 24 * 4);
}

#[test]
fn saem_rejects_unsupported_population_shape() {
    let (compiled, mut request) = fixture();
    request.random_effects.push(PopulationRandomEffect {
        parameter: "v".into(),
    });
    assert_eq!(
        compiled.fit_population_saem(&request).unwrap_err().code,
        pharmflux_core::ErrorCode::Unsupported
    );
}

#[test]
fn saem_preserves_complete_results_at_evaluation_limit() {
    let (compiled, mut request) = fixture();
    request.max_evaluations = request.subjects.len() * 4;
    let result = compiled.fit_population_saem(&request).unwrap();
    assert!(matches!(result.status, FitStatus::EvaluationLimit));
    assert_eq!(result.evaluations, request.max_evaluations);
    assert_eq!(result.subjects.len(), request.subjects.len());
    assert_eq!(
        result.objective_kind,
        PopulationObjectiveKind::SaemCompleteDataSurrogate
    );
    assert_eq!(
        result.fitted_observations.as_ref().unwrap().len(),
        request
            .subjects
            .iter()
            .map(|s| s.observations.len())
            .sum::<usize>()
    );
    assert!(result
        .subjects
        .iter()
        .all(|s| s.negative_log_likelihood.is_finite()));
    let diagnostics = result.saem_diagnostics.unwrap();
    assert_eq!(diagnostics.evaluated_proposals, 0);
    assert_eq!(diagnostics.accepted_proposals, 0);
    assert_eq!(diagnostics.final_step_size, None);
    assert_eq!(diagnostics.final_relative_change, None);
    assert_eq!(diagnostics.recent_relative_change_median, None);
    assert_eq!(diagnostics.recent_window, 0);
}

#[test]
fn saem_diagnostics_preserve_exact_gaussian_random_intercept_oracle() {
    // log(theta * exp(eta)) = log(theta) + eta, so the marginal Gaussian
    // MLEs are theta=exp(mean(y)) and Omega=mean((y-mean)^2)-sigma^2.
    let mut model: ModelDocument = serde_json::from_str(include_str!(
        "../../../conformance/models/synthetic-one-compartment.json"
    ))
    .unwrap();
    model.parameters[0].default = quantity(1., "1");
    model.parameters[1].default = quantity(1., "h");
    model.states[0].unit = "1".into();
    model.states[0].initial = pharmflux_core::model::Initial::Quantity {
        quantity: quantity(0., "1"),
    };
    model.states[0].dosing = None;
    model.outputs[0].unit = "1".into();
    model.outputs[0].expression = pharmflux_core::expression::Expr::call(
        pharmflux_core::expression::Function::Log,
        vec![pharmflux_core::expression::Expr::symbol("cl")],
    );
    let compiled = CompiledSensitivityDocument::compile(&model, &["cl".into()]).unwrap();
    let mut simulation: SensitivityRequest = serde_json::from_str(include_str!(
        "../../../conformance/requests/synthetic-sensitivity.json"
    ))
    .unwrap();
    simulation.with_respect_to = vec!["cl".into()];
    simulation.derivative_absolute_tolerances.clear();
    simulation.derivative_absolute_tolerances.insert(
        "cl".into(),
        std::collections::BTreeMap::from([("central".into(), quantity(1e-12, "1"))]),
    );
    simulation
        .run
        .parameters
        .insert("cl".into(), quantity(1., "1"));
    simulation
        .run
        .parameters
        .insert("v".into(), quantity(1., "h"));
    simulation.run.regimen.administrations.clear();
    let values = [-0.3, -0.2, -0.1, 0., 0., 0.1, 0.2, 0.3];
    let subjects = values
        .into_iter()
        .map(|value| GaussianObjectiveRequest {
            simulation: simulation.clone(),
            observations: vec![GaussianObservation {
                row: 0,
                output: "cp".into(),
                value: quantity(value, "1"),
                error: GaussianError {
                    additive_sd: quantity(0.05, "1"),
                    proportional_sd: 0.,
                },
            }],
            priors: vec![],
        })
        .collect();
    let request = PopulationFitRequest {
        subjects,
        fixed_effects: vec![FitBound {
            parameter: "cl".into(),
            lower: quantity(0.7, "1"),
            upper: quantity(1.3, "1"),
        }],
        random_effects: vec![PopulationRandomEffect {
            parameter: "cl".into(),
        }],
        omega: vec![vec![0.03]],
        omega_diagonal_bounds: vec![[0.005, 0.08]],
        eta_bound: 1.5,
        max_evaluations: 30_000,
        max_iterations: 800,
        gradient_tolerance: 1e-5,
        seed: 17029,
        burn_in_iterations: None,
        uncertainty: false,
        additive_error_fit: None,
    };
    let fitted = compiled.fit_population_saem(&request).unwrap();
    let exact_omega = values.iter().map(|value| value * value).sum::<f64>() / 8. - 0.05_f64.powi(2);
    assert!(
        (fitted.fixed_effects["cl"].value - 1.).abs() < 0.05,
        "{fitted:?}"
    );
    assert!(
        (fitted.omega[0][0] - exact_omega).abs() < 0.02,
        "{fitted:?}"
    );
    assert_eq!(
        fitted.status,
        FitStatus::Converged,
        "{:?}",
        fitted.saem_diagnostics
    );
    let diagnostics = fitted.saem_diagnostics.unwrap();
    assert_eq!(diagnostics.burn_in_iterations, 800 / 3);
    assert_eq!(diagnostics.recent_window, 20);
    assert!(diagnostics.recent_relative_change_median.unwrap() > request.gradient_tolerance);
    assert!(diagnostics.accepted_proposals > 0);
    assert!(diagnostics.accepted_proposals < diagnostics.evaluated_proposals);
}
