# Gaussian likelihood for narrow fitting

`CompiledSensitivityDocument::gaussian_objective` evaluates a negative log
likelihood and its gradient in the selected parameters' declared numeric units.
This objective supports the native individual and pooled fitters described below.

The request contains a sensitivity run, an observation table and optional
independent Gaussian priors. Each observation names an output and a zero-based
row in the run's ordered result grid, plus a measured quantity and error model.
Duplicate rows represent independent replicate measurements. Omit missing
measurements. Correlated residuals, censored observations, and lognormal
observation errors are outside this public contract.

The error coefficients are known constants for this operation:

- `additive_sd` is a quantity in the observable's units.
- `proportional_sd` is dimensionless.
- Both are nonnegative and cannot both be zero.

For prediction μ, variance is `additive_sd² + (proportional_sd · μ)²`.
The score includes the normal density normalization and the variance derivative,
not only a weighted residual term. Pure proportional error at zero prediction
is a domain error. Data and additive standard deviations convert to the model's
declared output unit before evaluation. Density values use that unit basis.
Optional priors name selected parameters and provide mean and positive standard
deviation quantities. Priors are normal densities on natural parameter values
in the model's declared parameter units. No parameter-coordinate Jacobian is
added. No priors means a maximum-likelihood objective. Residual coefficients
are not estimated by this request.

The result separates negative log likelihood, negative log prior and their
sum, and returns named gradient values with inverse-parameter units. Execution
identity binds the underlying sensitivity request, data, priors and budgets.
Invalid units, nonfinite/zero variance, unknown outputs, out-of-range rows,
duplicate prior names, overflow and exhausted storage budgets fail without a
partial objective. Observation storage is reserved before simulation.

The executable guides check synthetic parameter recovery with additive
Gaussian error. A successful synthetic example does not establish recovery
for every equation class, dataset, or sampling design.

## Native bounded individual fitting

`CompiledSensitivityDocument::fit_gaussian` minimizes the objective with
projected inverse BFGS and an Armijo line search. Each selected parameter needs
an explicit finite lower/upper quantity pair, in selection order. The search box
must contain the initial value and remain inside declared model bounds. Fixed
parameters cannot be selected for fitting.

This API explicitly uses affine coordinates that map each natural parameter
box to [0, 1]. It does not use the model's log/logit/probit transform metadata.
Priors remain densities on natural parameters, with no coordinate Jacobian.
The gradient tolerance applies to the infinity norm of the projected gradient
in these affine coordinates and therefore depends on search-box widths.

The result distinguishes convergence, evaluation limit, iteration limit and
line-search failure, preserving the best accepted parameters and their score.
A nonconverged result is not a successful fit. Domain, solver and invariant
failures during trial steps cause backtracking; malformed requests and exhausted
simulation budgets return errors. Failed trials count toward the evaluation
limit. Each evaluation uses the requested simulation budgets; total work is
bounded by those budgets times the evaluation limit (at most 10,000). Iterations
are separately bounded by 10,000. No global-optimum or identifiability claim is
made. Integration tolerances must support the requested gradient accuracy;
line-search failure near an optimum can indicate integration noise.

Recovery checks use independent analytic observations for clearance and volume,
including repeated boluses and a reset. Further checks cover a boundary optimum,
a normal-prior MAP fit, evaluation exhaustion, invalid-unit rejection and
recovery. These use tightened integration tolerances (rtol 1e-12). The fitter has native Rust and CLI interfaces; transform-based search and uncertainty estimates remain open.

## Pooled fitting

`fit_pooled_gaussian` accepts a primary `GaussianFitRequest` and up to 1,023
additional subject objectives. All subjects use the same compiled model and
share the selected parameters. Each subject retains its own regimen, covariates,
initial-state overrides, unselected parameter values, observations, residual
coefficients and simulation budgets. Selected initial parameter values must
agree after unit conversion; conflicting values are rejected. The primary
request supplies the shared search box and optimizer controls.

Declare priors only in the primary objective. Additional-subject priors are
rejected so the shared prior contributes exactly once. Likelihoods and gradients
sum over subjects in request order. Each aggregate evaluation counts once toward
the optimizer evaluation limit, and its maximum solver work is the sum of the
subject budgets. Subjects execute sequentially to keep the accumulation order
and peak simulation storage predictable. Independent fits/tests may run in
parallel.

The result retains the total likelihood, one prior, combined gradient and total
observation count. Aggregate execution hashes bind every subject's execution
identity, and the fit hash binds the full pooled request. A subject failure
cannot return a partial pooled score. Pooled MLE/MAP tests recover clearance and
volume from two different dose regimens, check score/gradient sums and prior
counting, and reject conflicting initial values, additional priors and invalid
subject units. This is shared-parameter pooled estimation; random effects and
population NLME estimation remain outside this API.

## Versioned file interface

The CLI accepts `pharmflux fit MODEL REQUEST.json`, where MODEL is model JSON
or `.pfx` text. JSON models permit 4,000,000 bytes; text and requests permit
1,000,000 bytes. The request schema
is `pharmflux.fit/v0.1`. Its `problem` contains `kind` (`individual` or `pooled`)
and `request` (the corresponding typed request above). A complete executable
example is `conformance/requests/synthetic-fit.json`.

Results use `pharmflux.fit-result/v0.1` with a `fit` object. Inspect `fit.status`:
only `converged` meets the requested optimizer criterion. Exit code zero means
the operation returned a result, including a nonconverged result; it does not
assert successful estimation. Invalid input and scientific execution failures
retain the CLI's structured stderr and nonzero exit behavior. Files are never
modified. Fit hashes bind the envelope version and problem kind in addition to
the optimizer request. Generated request/result schemas are in `spec/schema`.

The Python `CompiledSensitivities` class exposes `fit(request)` (dict result) and
`fit_json(request)` (JSON text). Both accept this same versioned envelope. Compile
with the request's selected parameter order, then reuse the immutable instance
across independent fits or threads. Native parsing, fitting and serialization
release the GIL; requests retain the one-million-byte bound. `PharmfluxError`
preserves structured failures. Nonconverged fits return an explicit status.
