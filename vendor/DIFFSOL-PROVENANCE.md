# Vendored Diffsol source

`diffsol/` starts from [Diffsol 0.17.1](https://github.com/martinjrobins/diffsol)
at upstream commit `bc6256a49e16779e77b58b786f0a993ad4242b63` (MIT license).
The public source export keeps the solver and linear-algebra crates used by
PharmFlux. It omits upstream examples, the C API crate, and benchmark drivers.
The C inputs still needed by Diffsol's optional Sundials feature are retained.

The solver extension changes these upstream files:

- `diffsol/crates/diffsol/src/lib.rs`
- `diffsol/crates/diffsol/src/error.rs`
- `diffsol/crates/diffsol/src/ode_solver/mod.rs`
- `diffsol/crates/diffsol/src/ode_solver/problem.rs`
- `diffsol/crates/diffsol/src/ode_solver/runge_kutta.rs`

It adds `diffsol/crates/diffsol/src/ode_solver/rosenbrock23.rs` and
`diffsol/crates/diffsol/src/ode_solver/rodas5p.rs`. The changes add two
identity-mass ODE methods, their storage and solver constructors, an
error-order option in the step controller, and a typed unsupported-output
error. Rosenbrock23 uses the ode23s quadratic dense output; Rodas5P uses its
fourth-order continuous extension.
Both enforce the configured minimum step size before attempting a step.

The public export also changes `diffsol/Cargo.toml`,
`diffsol/crates/diffsol/Cargo.toml`, and `diffsol/crates/diffsol/build.rs` to
remove the C API workspace member and unused benchmark targets. The vendored
README is shortened for this export. PharmFlux supplies the workspace
`Cargo.lock`; upstream examples and documentation remain at the pinned
commit above.

PharmFlux hashes the vendored source into its compiler identity and reports
that hash in `backend_version`. See `diffsol/LICENSE.txt` for the upstream
license.
