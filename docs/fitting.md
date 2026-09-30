# Parameter fitting

PharmFlux 0.1.0 includes the following bounded fitting contracts. They use
unit-bearing observations and explicit simulation settings.

| Request kind | Supported scope |
| --- | --- |
| `individual` | Gaussian observation likelihood and optional independent normal priors; fitted parameters have finite bounds. |
| `pooled` | Shared fitted parameters across subjects with separate regimens and observations. |
| Scalar fit | Derivative-free search for one bounded parameter through `CompiledModel.fit_scalar()` or its scalar request. |
| `laplace` | Population fit with one lognormal random effect and a Laplace approximation. |
| `focei` | Population fit with one or two lognormal random effects; the two-effect correlated covariance is fixed. |
| `saem` | Population fit with one or two lognormal random effects within the supported request shapes. |

Population parameters follow `individual_parameter = fixed_effect * exp(eta)`.
The covariance matrix `omega` is in the order of `random_effects`. A bounded
common additive residual standard deviation can be estimated only for the
specified one-effect FOCEI and SAEM contracts. Consult the versioned types in
[`fit.rs`](../crates/pharmflux-core/src/fit.rs) for exact fields and validation.

## Run an individual fit

The included synthetic request estimates clearance and volume:

```sh
cargo run --release --locked -p pharmflux-cli -- fit \
  conformance/models/synthetic-one-compartment.json \
  conformance/requests/synthetic-fit.json > fit.json
```

In Python, after [installing from source](installation.md):

```python
import json
from pathlib import Path
import pharmflux

source = Path("conformance/models/synthetic-one-compartment.json").read_text()
request = json.loads(Path("conformance/requests/synthetic-fit.json").read_text())
result = pharmflux.CompiledSensitivities(source, ["cl", "v"]).fit(request)
print(result["fit"]["status"], result["fit"]["parameters"])
```

The request nests a `pharmflux.sensitivity/v0.1` request and a
`pharmflux.run/v0.1` simulation request. Observation `row` indexes that
simulation's ordered output grid; each row names an output, measured value,
and additive/proportional Gaussian error. The result's `fit` member contains
status, estimates, objective, evaluation counts, and a request hash.

## Population results and data

Population fit results use the `population` member. Inspect `status`,
`fixed_effects`, `omega`, subject effects, and `objective_kind` together;
objectives from different methods have different meanings and must not be
compared as interchangeable numbers. When available,
`fitted_observations` records row-aligned observed values and population and
individual predictions. Optional FOCEI uncertainty is local curvature from
the observed outer objective, not an interval calibration. SAEM's
`saem_diagnostics` reports proposal and stochastic-drift information; it is
not a convergence certificate.

Python's `population_fit_dataset()` converts an explicitly mapped tabular
record subset into a population request and returns source-row mappings.
Pass that object to `CompiledSensitivities.fit_population_data()`. Its API is
declared in the [Python type stub](../bindings/python/python/pharmflux/__init__.pyi).
The synthetic population tests under
[`crates/pharmflux/tests/`](../crates/pharmflux/tests/) show supported request
shapes. Check bounds, residual-error assumptions, identifiability, and
scientific plausibility before using any fitted result.
