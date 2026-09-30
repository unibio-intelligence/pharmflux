# Installation and first run

PharmFlux 0.1.0 is a source release. Start from a checkout of this repository.
Rust builds use the pinned `Cargo.lock`; the first build downloads dependencies.
Rust and Cargo are required for the command-line and Python native builds.

## Command-line program

From the repository root:

```sh
cargo build --release --locked -p pharmflux-cli
./target/release/pharmflux --version
./target/release/pharmflux validate conformance/models/synthetic-one-compartment.pfx
./target/release/pharmflux run \
  conformance/models/synthetic-pbpk-24.pfx \
  conformance/requests/synthetic-pbpk-24.json > simulation.json
```

The executable writes results to standard output and structured errors to
standard error. `./target/release/pharmflux --help` lists operations. The
[model guide](modeling.md) explains the request format.

## Rust library

This release builds from a workspace; its crates are not published to a Rust
registry. In a Rust project, use path dependencies into a checkout:

```toml
[dependencies]
pharmflux = { path = "/path/to/pharmflux/crates/pharmflux" }
pharmflux-core = { path = "/path/to/pharmflux/crates/pharmflux-core" }
```

Replace `/path/to/pharmflux` with the checkout location. `pharmflux-core`
contains versioned requests and results; `pharmflux` contains parsing,
compilation, simulation, and fitting. See their public Rust APIs and the
[synthetic conformance cases](../conformance/) for examples.

## Python package

Python 3.10 or newer is required. From the repository root, in an isolated
Python environment:

```sh
python -m pip install .
python -c "import pharmflux; print(pharmflux.__version__)"
```

This builds the native extension from the checked-out Rust source; no wheel is
published by these instructions. `CompiledModel` runs simulations and scalar
fits. `CompiledSensitivities` runs sensitivity and individual or population
fits. See [fitting](fitting.md) for a complete example.

## WebAssembly and JavaScript

Install the Rust `wasm32-unknown-unknown` target and
[`wasm-pack`](https://rustwasm.github.io/wasm-pack/installer/), then build from
the repository root:

```sh
rustup target add wasm32-unknown-unknown
cd bindings/wasm
wasm-pack build --target web --out-dir pkg/engine
```

The generated `bindings/wasm/pkg/engine/` directory is consumed by the
supplied JavaScript worker; generated assets are not committed. Serve the
JavaScript module and `.wasm` file from your application. The package in
`bindings/wasm/` is marked private and is not published to a JavaScript
registry. See [browser integration and privacy](browser.md) before embedding
it. The direct WASM bindings contain additional operations; the supplied
browser worker currently exposes simulation only.

## R bridge status

`bindings/r/src/rust/` contains a low-level bridge for R integration. This
release does not contain a complete installable R package; use the Rust,
command-line, Python, or browser interfaces above.
