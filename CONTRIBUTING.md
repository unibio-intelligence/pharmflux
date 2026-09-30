# Contributing to PharmFlux

PharmFlux welcomes focused bug fixes, documentation, synthetic model examples,
scientific regression tests, and new capabilities within a clearly stated
supported scope.

## Start with a reproducible case

For a bug or scientific discrepancy,
[open an issue](https://github.com/unibio-intelligence/pharmflux/issues) with a
minimal model and request, the expected and observed behavior, units, dosing event order, solver
and tolerances, PharmFlux version, and operating system. Use synthetic or
shareable data. If a proposed feature changes a model format, solver contract,
or fit method, discuss its intended scope in an issue first.

## Develop and check a change

Install Rust and Cargo. Python changes need Python 3.10 or newer; building
the native Python package also needs the build dependencies declared in
`pyproject.toml`. From the repository root:

```sh
cargo fmt --all -- --check
cargo test --workspace --locked
node --test bindings/wasm/js/client.test.mjs
```

For a Python binding change, build and install the package in an isolated
Python environment, then run the relevant tests in `bindings/python/tests/`.
The Rust workspace tests cover the library and command-line interface.

A scientific change should include an independent check where practical: an
analytic solution, conservation or unit invariant, or a small synthetic
reference case. Test event boundaries and failure behavior when they are
relevant. State the method's supported equation, error, and random-effect
scope in the code and documentation. Keep solver selection and scientific
errors explicit rather than silently changing methods or returning partial
results.

## Submit a change

Fork the public repository, work on a focused branch, and open a pull request
against `main`. Explain the behavior the change adds or corrects. Include the
model and request needed to reproduce a scientific result, the checks you ran,
and any limits that remain. Update the README or examples when a public API or
supported scope changes. Do not include patient data, credentials,
confidential research data, or third-party assets without redistribution rights.

The `vendor/diffsol/` tree is pinned upstream solver source. For solver changes,
identify the affected upstream version, add numerical regression coverage,
and update [its provenance record](vendor/DIFFSOL-PROVENANCE.md). Solver changes
may also be suitable for contribution to the
[Diffsol project](https://github.com/martinjrobins/diffsol).

PharmFlux code and synthetic fixtures are licensed under
[Apache 2.0](LICENSE). Vendored Diffsol retains its upstream MIT license.
