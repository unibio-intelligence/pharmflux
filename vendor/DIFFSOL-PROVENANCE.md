# Vendored Diffsol source

`diffsol/` contains the Diffsol 0.17.1 source from
[martinjrobins/diffsol](https://github.com/martinjrobins/diffsol), based on
upstream commit `bc6256a49e16779e77b58b786f0a993ad4242b63` and extended
through commit `076d21ccc3e227e5893aaf38af9f4a31856dff2f`. The solver
implementation files match that commit byte for byte. The public
export omits the upstream workspace `Cargo.lock` and examples, so its root
`Cargo.toml` removes the examples workspace member; PharmFlux supplies its
own lockfile. The vendored README and benchmark note have small documentation
updates for this trimmed source export.

The PharmFlux public source export includes the Diffsol solver crates, root
workspace metadata, README, citation, and MIT license needed to build
PharmFlux. Upstream examples, documentation, and CI files remain available
at the commit above.

The extension adds Rosenbrock23 and Rodas5P methods for identity-mass ODEs,
refreshes their Jacobians for each step, and corrects the Rosenbrock23
third-stage time-derivative scale to `gamma*h`. The method implementations,
coefficients, limits, and source attributions are documented in
`diffsol/crates/diffsol/src/ode_solver/rosenbrock23.rs` and
`diffsol/crates/diffsol/src/ode_solver/rodas5p.rs`.

Diffsol is MIT-licensed; see `diffsol/LICENSE.txt`. The PharmFlux build hashes
the vendored files into its compiler identity and reports their source hash
in `backend_version`.
