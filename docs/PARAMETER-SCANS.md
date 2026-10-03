# Bounded parameter scans

`CompiledDocument::scan` and `pharmflux scan MODEL REQUEST.json` evaluate an
explicit ordered list of 1–25 parameter points. The shared run defines regimen,
solver, initial-state overrides, covariates and base parameter values. Each
point has a unique 1–128 byte ID and a parameter quantity map that overrides the
base map. Defaults still supply unspecified model parameters. There is no
implicit grid expansion. The fixture `conformance/requests/synthetic-scan.json`
shows three clearance values for the owned one-compartment model.

Requests use `pharmflux.scan/v0.1`; results use `pharmflux.scan-result/v0.1`.
Every point returns its ID, merged explicit parameter map and full simulation
result with execution identity. The parameter map preserves supplied units;
normalized bound values are represented by the simulation identity. The scan
hash binds the full ordered request and scan algorithm. Points with equal
scientific inputs may have equal simulation identities while retaining distinct
point IDs. Point order is preserved.

All point parameter/initial-state bindings are validated before any simulation.
Run-level protocol and execution validation still occurs in the ordinary run
path. Each point receives the smaller of the run budget and its equal share of
the total scan budget, independently for callbacks, events and numeric output
values. This conservative reservation bounds aggregate work and retained numeric
results. Unused budget from one point is not transferred to another. Execution
is sequential within one scan; independent scans can run concurrently.

Any point failure aborts the operation without returning partial results. The
error expression identifies its point. CLI files are never modified and retain
the one-million-byte request limit (JSON model documents permit four million
bytes). Generated request/result schemas describe
structure; runtime checks enforce point limits, units, domains and budgets.

Qualification covers independent analytic trajectories, mixed-unit overrides,
ordered results, distinct bound identities, malformed points, excess/empty point
lists, aggregate-budget exhaustion, recovery, CLI parity and independent schema
validation. This scan operation does not itself estimate global sensitivity indices.

## Python execution

`CompiledModel.scan(request)` returns a dict; `scan_json(request)` returns
JSON. Both accept the versioned envelope, release the GIL during native
execution, and permit sharing a compiled model across independent scans.
See the Python and R user guides for executable synthetic examples.
