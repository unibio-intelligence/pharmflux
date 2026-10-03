# Scientific run and result envelopes

Draft versions: `pharmflux.run/v0.1` and `pharmflux.result/v0.1`. The authoritative
Rust records are in `pharmflux-core::run`; they have no solver dependency.
Generated structural JSON Schemas are in `schema/`; see `schema/README.md`
for validation limits and reproduction. The executable request example is
`conformance/requests/synthetic-mixed-scales.json`.

## Request

A request contains a schema, optional correlation `request_id`, explicit solver,
explicit budgets, parameter and initial-state overrides as quantities, and the unit-bearing
regimen. The regimen carries covariates, fixed changes, sampling sides,
administrations, resets and tolerance settings. Unknown fields are rejected.
The ID, when provided, is 1–128 bytes. Native and WASM execution share unsigned
32-bit budget limits for callback counts, output values and events.

The implemented ODE solvers are `diffsol_bdf`, `diffsol_tsit45`,
`diffsol_esdirk34`, `diffsol_tr_bdf2`, `diffsol_rosenbrock23`, and
`diffsol_rodas5p`. The last two use locally vendored Diffsol source and accept
identity-mass ODEs. Fitting and sensitivity operations currently require BDF.
Unknown names fail decoding and do not select a fallback. The optional seed field reserves the stochastic contract,
but every non-null seed is currently rejected as unsupported. No stochastic
behavior or implicit random seed is implemented.

`CompiledDocument.execute` returns a result only after the complete scientific
operation succeeds. Failure is the existing typed error with no partial result.
The native `execute_request` example accepts a model text/JSON file and a run
JSON file and emits the result JSON. It is an offline execution entry point,
also available through the implemented `pharmflux run MODEL REQUEST.json`
command. See ../docs/CLI.md for the CLI binding and packaging limits.

## Result

The native result contains schema, correlation ID, scientific identity, model
ID/description, time unit, explicit pre/post sides, times, and named output
columns with units. Display fields are separate from the identity record.

The browser worker accepts `{kind: "run", source, request}`. It transfers typed
numeric arrays and returns the same metadata and identity, with `names`,
`units` and `columns` as the columnar view. Browser sides are bytes (0 pre,
1 post); native JSON uses the names `pre`/`post`. The scientific path rejects
additional legacy protocol/regimen/override/callback fields outside the request
envelope to avoid conflicting settings.

Legacy conformance `kind: "text"` payloads remain available; they do not acquire
a scientific execution identity merely because the new envelope exists.

## Canonical hashes

Core canonicalization uses RFC 8785 JCS (`serde_jcs` 0.2.0), followed by SHA-256
(`sha2` 0.10.9). Hash strings are `sha256:` followed by 64 lowercase hex digits.
Tests cover canonical key ordering including UTF-16 order, numeric negative
zero normalization and a known SHA-256 vector. Integer JSON values outside the
exact interoperable range ±(2^53−1) are rejected; large identifiers should be
strings. Scientific binary64 values use the JCS number representation.

- `model_content_hash`: canonical serialized, validated model document. This
  includes its descriptions and declared units. Whitespace and formatter/JSON
  round trips do not change it; changing document content does. It is not a
  proof that different equations or unit spellings are mathematically equivalent.
- `bound_value_hash`: the complete numeric binding buffer, including defaults,
  converted covariates and derived individual parameters, plus the resolved initial
  state vector in declared state units and compiled order.
  Its interpretation is bound to the model content hash.
- `run_request_hash`: schema, regimen, solver, budgets and seed. Correlation IDs
  are excluded. Parameter and initial-state override spelling is excluded because the
  complete resolved bindings and initial vector are recorded by `bound_value_hash`. Covariate changes and their
  times remain in the regimen hash. Semantically equivalent but differently
  written regimens can have different hashes; no semantic equivalence is assumed.

## Compiler and backend provenance

Identity also records specification version, compiler package version, compiler
source SHA-256, rustc version, target, profile, encoded-rustflags SHA-256,
backend/version and algorithm. Native and WASM can share the scientific content
hashes while correctly recording different build targets.

The current workspace build script hashes ordered relative paths, lengths and
bytes for Rust core/engine/WASM sources, build script, workspace/crate manifests,
lockfile, Cargo configuration, and the complete vendored Diffsol source. This distinguishes development source changes
under the same package version. It is a source/build fingerprint, not a binary
signature, OCI attestation or proof of installed dependencies. Standalone crate
packaging and release attestation remain unfinished. The backend pins Diffsol
0.17.1 plus the local Rosenbrock methods; `backend_version` includes the exact
vendored-source SHA-256. See `vendor/DIFFSOL-PROVENANCE.md` for the upstream
base, local patches, license, and review record.

References: [RFC 8785](https://www.rfc-editor.org/rfc/rfc8785),
[serde_jcs](https://docs.rs/serde_jcs/0.2.0/serde_jcs/),
[SHA-2 implementation](https://docs.rs/sha2/0.10.9/sha2/).

## Initial-state overrides

Optional `initial_states` maps declared state names to quantities, for example
`{"central": {"value": 6, "unit": "mg"}}`. It is a run-level initial condition,
not a time-zero reset. Compatible units convert once into state storage units;
unknown names, incompatible units and non-finite values fail.

An explicit override replaces only that state's initial expression. Other initial
expressions evaluate from the current parameter/covariate binding. The overridden
expression still undergoes structural and unit compilation but is not evaluated
at binding. This allows a valid explicit initial condition even when that
expression would be undefined for the run's parameter values. Parameter bounds,
derived parameters and dose scales still undergo their ordinary validation.

Initial values precede time-zero pre-observations, covariate changes, resets and
doses. Covariate changes carry the trajectory forward without reinitialising it.
Declared invariants are checked before any event can hide an invalid initial
value. Native `bind_with_initial_states(parameters, covariates, initial_states)`
and subsequent parameter `rebind` preserve explicit initial overrides; creating a
new binding/run with an empty map restores default initial expressions.

The bound-value hash now includes resolved initial values. Initial quantities that resolve to the same numeric vector, including explicit
values equal to defaults, share the same scientific identity for otherwise
identical runs. This intentionally changes historical
bound hashes from earlier draft compilers; their recorded build fingerprints
remain distinct. No stable hash compatibility was promised for those drafts.

## Absolute simulation start

`regimen.start` is an optional quantity. Omission means zero in the model time
unit. A raw protocol's numeric `start` similarly defaults to zero. The state
initial condition applies at start, before start-time events. All samples,
observations, reset/covariate changes and expanded administrations use absolute
model time and must lie within the inclusive window. Infusions must end within
it. No pre-start dose, active infusion, or state history is inferred.

Start and end must be finite, with end strictly greater than start. Negative
windows are supported by the engine. The public raw `regimen::expand` utility
retains its zero-start meaning; `expand_in_window` accepts an explicit start.
Time expressions are not shifted, so nonautonomous models keep their meaning.

## Fixed solver checkpoints

`regimen.checkpoints` is an optional array of time quantities (default empty).
Each must be finite, dimensionally a time, and inside the declared start/end
window. Checkpoints force integration boundaries and solver restarts without
changing state, parameters, or active inputs. Coincident checkpoints and events
share a boundary; covariate binding selection still uses covariate change times.
Each supplied checkpoint consumes the event budget, including duplicates.
Ordinary observation times remain dense-output requests and do not force restarts.
Checkpoint quantities are included in the scientific request identity.

## Piecewise smooth state Jacobians

`abs(x)` uses derivative -1 for x < 0, 0 at x = 0, and 1 for x > 0,
multiplied by the derivative of its argument. `min` and `max` use the
derivative of the selected argument; exact ties select the first argument
in source order. These are local Jacobian conventions at nonsmooth points,
not a claim of classical differentiability there. All primal arguments retain
their domain checks. Arbitrary state-dependent conditionals remain unsupported.

Native tests check multiple arguments, ties, invalid unselected primal domains,
and BDF integration across a branch boundary against analytic integrals.
Browser qualification and trajectory comparisons are recorded separately.

## Hill state Jacobian

For positive x, `hill(x,n,k)` is differentiated through all three arguments.
The implementation evaluates both H(x,n,k) and H(k,n,x), avoiding cancellation
from subtracting a rounded-to-one Hill value. At negative x all partials vanish.
At zero, the n and k partials vanish; the x partial is zero for n>1, 1/k for
n=1, and the positive-x slope at x=1e-8*k for n<1. This last convention
applies to the BDF state Jacobian and does not modify the evaluated RHS.

The expression-only parameter derivative API rejects zero-boundary x
derivatives for n<=1 because it cannot validate a model's nonnegative
invariant. End-to-end parameter sensitivity remains a separate M3 capability.
