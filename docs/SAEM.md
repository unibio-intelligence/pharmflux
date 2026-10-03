# Population SAEM fitting

The public native SAEM contract fits one or two positive model parameters with
lognormal subject effects and Gaussian observation error. The runnable guides
use a deliberately short iteration budget; inspect status and objective kind
before interpreting estimates.

## Request and observation model

Use a `pharmflux.fit/v0.1` envelope with `problem.kind: "saem"`. Each subject
supplies its own sensitivity run, regimen, fixed model inputs, and observation
rows. The population request declares fixed-effect bounds, random effects,
eta covariance and variance bounds, seed, and work limits. Use the generated
fit-request schema together with runtime validation.

Fixed residual error uses `additive_sd` in the output's units and a
dimensionless `proportional_sd`. Their prediction-dependent variance is
`additive_sd² + (proportional_sd * prediction)²`; it must be positive.
Subject effects satisfy `individual_parameter = fixed_effect * exp(eta)`.
Lognormal subject effects do not imply lognormal observation errors.

`additive_error_fit` optionally estimates one common additive Gaussian SD.
Every included observation must use the selected output, zero proportional
error, and the declared starting SD. Supply positive, finite bounds with units.
The estimator supports this contract for one or two subject effects; it does
not estimate a general residual covariance matrix.

## Status, results, and limits

`burn_in_iterations` sets exploration length independently of
`max_iterations`. The request needs at least eight iterations and four
objective evaluations per subject; realistic work limits depend on the
cohort and model. The quickstart's eight iterations illustrate request and
result handling and are not a recommended estimation configuration.

A converged result reports `saem_marginal_objective_function_value`. A limited
run may report `saem_complete_data_surrogate`, which is not interchangeable
with the marginal objective. Inspect the returned subject effects, covariance,
fit status, stochastic diagnostics, and fitted observations when present.
Record the seed, examine stability across seeds, and verify recovery for the
intended equation class, sampling design, and cohort.

This public version does not implement SAEM uncertainty: `uncertainty: true`
returns an unsupported-operation error. It also does not accept censored
observations, lognormal observation-error fields, or a fitted-covariate request
field. Model-authored covariates may be supplied in the explicit subject run;
the row adapter does not infer them. Subject priors are unsupported. Successful
computation does not establish identifiability or biological validity.
