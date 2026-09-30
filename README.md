# PharmFlux

PharmFlux is an open-source Rust library for pharmacology models. It compiles
unit-aware model definitions, simulates dosing regimens, and fits supported
individual, pooled, and population models. The same engine is available to
Rust, command-line, Python, and browser applications.

This is an early **0.1.0 source release**. The instructions below build it
from this repository.

## Documentation

- [Installation and first run](docs/installation.md): Rust, Python, and
  WebAssembly builds.
- [Models and simulation](docs/modeling.md): model format, units, dosing,
  solver selection, and result shape.
- [Fitting](docs/fitting.md): supported fit contracts, example requests, and
  result interpretation.
- [Browser integration and privacy](docs/browser.md): Web Worker setup and
  what local execution does and does not guarantee.

## What is included

- **Modeling and simulation:** text or JSON models, explicit units, boluses,
  infusions, resets, covariates, observations on either side of an event,
  analytic linear PK paths, and numerical ODE solvers. Available solver choices
  include BDF, Tsit45, ESDIRK34, TR-BDF2, Rosenbrock23, and Rodas5P. Solver
  capabilities differ; unsupported combinations return an error.
- **Fitting:** individual and pooled Gaussian fits, a bounded single-parameter
  fit, and scoped FOCEI, Laplace, and SAEM population fits. Population methods
  support specific random-effect and residual-error forms; see the request
  types and tests before using a new model shape. A converged fit still needs
  scientific assessment of the data and parameter identifiability.
- **Analysis:** forward sensitivities where supported, bounded parameter scans,
  and Morris screening.
- **Reproducibility:** versioned requests and results, input validation,
  execution identity, explicit solver settings, and structured errors.

## Try it from the source checkout

Install Rust and Cargo. From the repository root, run the
included synthetic examples:

```sh
cargo run --release --locked -p pharmflux-cli -- run \
  conformance/models/synthetic-pbpk-24.pfx \
  conformance/requests/synthetic-pbpk-24.json > simulation.json

cargo run --release --locked -p pharmflux-cli -- fit \
  conformance/models/synthetic-one-compartment.json \
  conformance/requests/synthetic-fit.json > fit.json
```

The fit example estimates clearance and volume. Its status and estimates are
under `fit` in `fit.json`. The examples use illustrative parameters and make no
clinical prediction. [Synthetic fixture terms](conformance/FIXTURE-LICENSE.md)
apply to the included model and request files.

For Python 3.10 or newer, install from this checkout with
`python -m pip install .`. The package uses the same native Rust engine:

```python
import json
from pathlib import Path
import pharmflux

model = Path("conformance/models/synthetic-one-compartment.json").read_text()
request = json.loads(Path("conformance/requests/synthetic-fit.json").read_text())
result = pharmflux.CompiledSensitivities(model, ["cl", "v"]).fit(request)
print(result["fit"]["status"], result["fit"]["parameters"])
```

Python also exposes `CompiledModel.run()` and `CompiledModel.fit_scalar()`.
`population_fit_dataset()` converts an explicitly mapped tabular record subset
into a population fit request while retaining source-row mappings. The
command-line `fit` operation accepts the same versioned fit request used by
the Rust and Python sensitivity APIs.

## Repository structure

| Path | Purpose |
| --- | --- |
| `crates/pharmflux-core/` | Model, regimen, unit, run, fit, and result types. |
| `crates/pharmflux/` | Model parser, compiler, simulation runtime, sensitivities, and fit algorithms. |
| `bindings/cli/` | `pharmflux` command-line interface. |
| `bindings/python/` | Native Python package and tabular population-fit helper. |
| `bindings/wasm/` | WebAssembly and JavaScript bindings for local browser execution. Generated WASM files are not committed. |
| `bindings/r/src/rust/` | Low-level R bridge source. A complete installable R package is not part of this release. |
| `conformance/` | Synthetic models, requests, and scientific regression cases. |
| `vendor/diffsol/` | Pinned solver source used by the Rust runtime. See [provenance](vendor/DIFFSOL-PROVENANCE.md). |

The WebAssembly source includes simulation, analysis, and fitting bindings.
The supplied browser worker currently exposes simulation; see the
[browser guide](docs/browser.md) for its contract and build steps. Model text,
requests, and results can remain on the user's device when an application uses
that worker without uploading them. The application controls its own network
requests, telemetry, and persistence, so it must preserve that boundary.
Diffsol's README and crate histories under `vendor/` describe that upstream
solver project; the PharmFlux interfaces are listed above.

## Checks and contributions

Run `cargo fmt --all -- --check`, `cargo test --workspace --locked`, and
`node --test bindings/wasm/js/client.test.mjs` from the repository root. The
public workflow runs these checks. See [CONTRIBUTING.md](CONTRIBUTING.md) for
scientific test expectations, issue reports, and contribution steps.

## License and attribution

PharmFlux code and synthetic fixtures are licensed under
[Apache 2.0](LICENSE). Vendored Diffsol retains its upstream MIT license;
its source, citation, and local changes are described in the
[provenance record](vendor/DIFFSOL-PROVENANCE.md).

Developed by [UniBio Intelligence](https://unibiointelligence.com).

## Acknowledgment

PharmFlux development received support from Anthropic's rare disease grant.
