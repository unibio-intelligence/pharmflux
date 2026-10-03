# Population SAEM fitting

This page describes the native Rust fitting contract. The runnable quickstart deliberately uses a short iteration budget. Inspect fit status and method-specific objective values before interpreting estimates.

`CompiledSensitivityDocument::fit_population_saem` fits one or two positive
model parameters with lognormal subject effects. Each subject supplies its own
regimen, observations, and fixed model inputs. The request binds the subject
order, initial fixed effects, eta covariance, parameter and variance bounds,
random seed, and work budget. A successful result has `status: "converged"`
and `objective_kind: "saem_marginal_objective_function_value"`. An evaluation
limit may return a complete-data surrogate; that value is not interchangeable
with the marginal objective.

## Observation models

The normal likelihood supports fixed additive and proportional standard
deviations, including both terms together. `error.model: "lognormal"` uses
zero additive SD and interprets `proportional_sd` as a dimensionless log-scale
SD. Predictions and values must then be positive. All coefficients are fixed
unless `additive_error_fit` explicitly estimates one common additive SD.

For a below-quantification observation, set `censoring: "left"` and put the
quantification limit in `value`. The likelihood uses the probability below
that limit. Censored observations can use fixed normal or lognormal error.
`additive_error_fit` can estimate a common additive normal SD with mixed
observed and left-censored rows. Its SAEM variance step uses the conditional
expected squared residual under the current prediction and censoring limit.
This estimate requires at least one uncensored observation and a recoverable
cohort; all-censored or weakly identified designs may converge at a bound.
Priors and correlated observation errors are outside this contract.

## Covariate effect

One-effect SAEM can estimate up to four exponential effects through
`covariate_effect`. Supply `parameter`, `covariate`, `subject_values` in the
same order as `subjects`, a `reference`, and finite `initial` and `bounds` for
each coefficient. Put effects after the first in `additional`. The subject
parameter is

`theta_subject = theta_reference * exp(sum_j coefficient_j * (value_j - reference_j) + eta)`.

Values must be finite and subject-constant. The centered design must have full
rank and more subjects than coefficients. The implementation requires fixed
residual error and does not combine fitted covariates with uncertainty.
Two-random-effect covariate inference, linear and power forms, and estimation
of coefficients on time-varying inputs remain outside this contract.

A model can already contain a stepwise time-varying covariate in its equations.
Each subject may then supply `regimen.covariates` and
`regimen.covariate_changes`; the SAEM fit estimates the selected population
parameter and eta covariance while the model-authored covariate coefficient
remains fixed. Changes take effect on the post-event side of their time point.
This path has a synthetic recovery test with different change times and values
across subjects. Linear interpolation remains unsupported by the native model.
Other population fitters reject the SAEM-specific fitted-covariate field.

## Uncertainty and limits

For a converged one- or two-effect fit with fixed residual error and no fitted
covariate, `uncertainty: true` estimates a local covariance for the fixed
effects and free eta covariance coordinates. It finite-differences the final
Laplace-approximated marginal objective and inverts its observed curvature.
Boundary estimates, singular curvature, and uncertified perturbed objectives
return an error. The covariance is an asymptotic local approximation; it does
not include stochastic-approximation variability, a sandwich correction, or
interval calibration. The two-effect calculation includes the eta covariance
between effects. Strongly correlated cohorts can produce singular or
boundary curvature; these return an error instead of standard errors.

The native tests cover independent likelihood and gradient calculations,
synthetic parameter recovery, censored rows with fitted additive SD,
lognormal residuals, two independent exponential covariates, model-authored
step covariates, and finite symmetric one- and two-effect uncertainty.
These checks do not establish recovery for every model equation,
sampling design, dose history, or parameter boundary.
