# Morris elementary-effect screening

`CompiledDocument::morris` evaluates explicit Morris trajectories through the
bounded parameter-scan engine. Inputs identify parameter ranges as quantities,
ordered trajectories of point IDs, an even grid level count, and one output
column/row. Each trajectory changes each selected parameter once, with one
parameter changed per step. Every point must belong to exactly one trajectory,
and every point explicitly supplies exactly the selected parameters. Other
model parameters and the regimen remain fixed through the shared scan request.

The normalized grid has p even levels (4–100), spacing 1/(p−1), and step magnitude
Δ=p/[2(p−1)]. Steps may go in either direction. Parameters are normalized using
user-supplied natural-unit ranges after quantity conversion. Grid checks allow
1e-8 error in grid-index coordinates and 1e-12 at normalized range endpoints to
accommodate arithmetic/unit-conversion roundoff. Grid coordinates are snapped
for design validation and step denominators; simulations use the supplied
parameter quantities without adjustment.

For each step, EE=(output_after−output_before)/signed_normalized_step. For each
parameter, the result retains all EEs, their signed mean `mu`, the mean of their
absolute values `mu_star`, and sample standard deviation `sigma` (denominator
n−1). Because inputs are normalized, these quantities have the selected output's
units. They depend on the specified parameter ranges and output selection.
The definitions follow the unscaled, ungrouped Morris metrics documented by
[SALib](https://salib.readthedocs.io/en/main/api.html). No bootstrap confidence
intervals or standard-deviation scaling are computed.

The engine requires at least two trajectories, 1–11 selected parameters and no
more than 25 scan points. The scan's conservative aggregate budget policy still
applies. Designs are checked before simulation. All runs must expose identical
time/side grids so an output row has the same meaning throughout the analysis;
parameter-dependent event grids that differ are rejected. Nonfinite effects or
summary overflow return structured errors. The result retains the complete
scan evidence and a hash of the analysis request and algorithm.

Tests compare effects and all three summaries with independently evaluated
one-compartment PK formulas, including forward/reverse steps and mixed units.
Malformed grid levels, incomplete/reused trajectories, multiple simultaneous
changes, wrong step sizes, incompatible units and invalid outputs/rows fail,
and a subsequent valid request succeeds.

This is a native Rust API accepting explicit or generated trajectories.
Optimizing trajectories, bootstrap intervals and broad model qualification
remain open. Two
trajectories establish a computable statistic, not adequate coverage of an
arbitrary nonlinear model. Interpretation requires a design that covers the
scientific range of interest; sigma can reflect nonlinearity and interactions
and does not separate them into variance contributions.

## Deterministic trajectory generation

`pharmflux_core::morris::generate_morris_design` accepts the shared run, quantity
ranges, trajectory count, even grid level count, a u32 seed, aggregate budgets
and output selection. It returns a complete `MorrisRequest` plus the seed,
generator name and design-request hash. It uses no solver, clock, global RNG,
threads or filesystem. Model-specific names, domains and execution validity
remain the analysis compiler's responsibility.

The generator implements ungrouped, unoptimized Morris sampling: independent
base grid indices and direction choices, with a Fisher–Yates parameter
permutation per trajectory. Its grid geometry follows the vanilla design
shown in [SALib's sampling implementation](https://salib.readthedocs.io/en/latest/_modules/SALib/sample/morris/morris.html).
It does not claim seed-for-seed agreement with NumPy/SALib.

Generator `morris_sha256_counter_v1` hashes the UTF-8 domain string
`pharmflux.morris.sha256-counter.v1`, followed by the seed as four little-endian
bytes and a zero-based counter as eight little-endian bytes. Each digest's
first eight bytes form a little-endian u64. Bounded integer draws reject values
below `2^64 mod n` before reducing modulo n, with at most 256 attempts per draw.
This makes the stream explicit and reproducible without a new RNG dependency.
For each trajectory, permutation draws precede per-parameter base/direction
pairs. Direction draw zero means increasing. Every generated point has a unique
trajectory/step ID; coordinates may repeat across independently drawn paths.

Requests exceeding 25 runs are rejected before allocation/generation. Ranges
must be finite, increasing and unit-compatible. Quantities use the lower bound's
unit; exact grid endpoints use the converted bounds directly. The design hash
binds all generator inputs and its algorithm name. The analysis hash binds the
resulting explicit request, so designs with identical execution points can be
scientifically equivalent while retaining separate generation provenance.

Qualification includes an independently calculated fixed-seed grid, exact
repeatability, seed-dependent identity, 16 seed/grid combinations through actual
simulation and analysis, and invalid range/grid/work-count rejection. This does
not establish that a small random design adequately covers any particular model.

## Versioned CLI execution

`pharmflux morris MODEL REQUEST.json` accepts model JSON or `.pfx` text and a
`pharmflux.morris/v0.1` envelope. `problem.kind` is `generated` or `explicit`;
`problem.request` is respectively a design-generation request or a complete
explicit analysis request. The owned example is
`conformance/requests/synthetic-morris.json`.

`CompiledDocument::execute_morris` runs the same envelope in Rust. Both paths
return `pharmflux.morris-result/v0.1`, containing the analysis result and an
optional `design` record. Generated designs preserve generator name, seed and
design-request hash; explicit inputs have null generation provenance. The
analysis identity binds the envelope version and full problem request. An
explicit design and its generated equivalent can have identical effects but
different analysis hashes, reflecting different request provenance.

Generated request/result schemas live in `spec/schema`. All nested fields are
closed; unknown versions or fields fail. Semantic limits, dimensional checks and
execution validity remain runtime checks. The CLI preserves its file-size
limit and structured stderr behavior. Tests establish exact CLI/library parity
for both modes, equal scientific effects for equivalent designs, and independent
schema validation of an actual emitted six-point result.

## Python execution

`CompiledModel.morris(request)` returns a dict; `morris_json(request)` returns
JSON. Both accept the versioned envelope and detach native execution from the
GIL. Compiled models may be shared across independent analyses. Request
limits and structured failures remain enforced. See the Python user guide
for an executable synthetic example.
