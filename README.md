# PharmFlux

**An open-source engine for building and simulating pharmacology models.**

PharmFlux helps researchers describe how drugs move through the body and how
biological systems respond. It is being developed for pharmacokinetic (PK),
pharmacodynamic (PD), and mechanistic modeling, from compartment models to
physiologically based pharmacokinetic (PBPK) and quantitative systems
pharmacology (QSP) models.

## What it does

Write a model as equations with explicit units, define a dosing regimen, and
simulate how concentrations and other model outputs change over time.
Use it to explore dose schedules, compare parameter choices, and build
reproducible simulation workflows.

- **Readable models:** describe parameters, states, equations, and outputs in
  a text format, with structured JSON for programmatic use.
- **A shared Rust engine:** bring the same model into command-line, Python,
  and browser workflows through native and WebAssembly bindings.
- **Explicit scientific inputs and outputs:** keep units, dosing events,
  solver settings, and execution provenance alongside the simulation.

The goal is to make pharmacology models easier to inspect, reuse, and integrate
into research software.

## Try it locally

With Rust installed, clone this repository and run the included synthetic
examples:

```sh
cargo run --release --locked -p pharmflux-cli -- run \
  conformance/models/synthetic-pbpk-24.pfx \
  conformance/requests/synthetic-pbpk-24.json > simulation.json

cargo run --release --locked -p pharmflux-cli -- fit \
  conformance/models/synthetic-one-compartment.json \
  conformance/requests/synthetic-fit.json > fit.json
```

The simulation result contains time points, observation sides, named outputs,
and execution identity. The fitting result reports a status and fitted
parameters. These examples use synthetic inputs; estimation methods have
explicitly supported model, error, and random-effect scopes.

The Python package can be built with `maturin` from the repository root.
`CompiledModel.fit_scalar` fits one bounded parameter from independent,
additive-Gaussian observations using native simulations when compiling
sensitivity equations is impractical. It reports convergence and the number
of simulation evaluations.
See `pyproject.toml` for its locked build configuration. Run
`cargo test --workspace --locked` to check the Rust workspace.

## Local execution and privacy

PharmFlux can run simulations on your own machine or directly in a browser
through WebAssembly (WASM). Browser applications can compute results locally
without uploading model inputs to a simulation server, and use Web Workers to
keep the interface responsive while a simulation runs.

This makes it possible to build interactive modeling tools that keep models
and data on the user's device. Any application embedding PharmFlux controls
its own data collection, storage, and network behavior.

## Current scope

This is an early source release. The simulation engine covers unit-aware
models, dosing and events, analytic linear PK, and explicit and stiff ODE
solvers. Individual and pooled Gaussian fitting, plus bounded forms of FOCEI,
Laplace, and SAEM, are available through the Rust API. These population
methods have limited supported model and random-effect shapes; check the
returned diagnostics and validate results for each scientific use case.

The ODE solver source is vendored from [Diffsol](vendor/DIFFSOL-PROVENANCE.md)
under the MIT license, with PharmFlux Rosenbrock extensions.

## Contribute

Try PharmFlux with your own models and tell us what works and what could be
better. Open an issue with a modeling use case, a question, or a feature
suggestion. Contributions of model examples, correctness tests, documentation,
bug fixes, and integrations are welcome. For larger changes, start with an
issue so we can discuss the approach.

## License

PharmFlux is licensed under either the
[MIT License](LICENSE-MIT) or the
[Apache License 2.0](LICENSE-APACHE), at your choice.

Developed by [UniBio Intelligence](https://unibiointelligence.com).
