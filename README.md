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

## Local execution and privacy

PharmFlux can run simulations on your own machine or directly in a browser
through WebAssembly (WASM). Browser applications can compute results locally
without uploading model inputs to a simulation server, and use Web Workers to
keep the interface responsive while a simulation runs.

This makes it possible to build interactive modeling tools that keep models
and data on the user's device. Any application embedding PharmFlux controls
its own data collection, storage, and network behavior.

## Contribute

Try PharmFlux with your own models and tell us what works and what could be
better. Open an issue with a modeling use case, a question, or a feature
suggestion. Contributions of model examples, correctness tests, documentation,
bug fixes, and integrations are welcome. For larger changes, start with an
issue so we can discuss the approach.

## License

PharmFlux is licensed under the
[Apache License 2.0](https://www.apache.org/licenses/LICENSE-2.0).

Developed by [UniBio Intelligence](https://unibiointelligence.com).
