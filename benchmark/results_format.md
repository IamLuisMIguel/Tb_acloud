# Machine-readable submission results, version 2

Submit `results.json` alongside `solution.py` and `report.md`. Generate it from the same implementation used by the evaluator, with this terminal interface:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 solution.py \
  --inputs benchmark/inputs.json --output results.json
```

The command must execute the required calculations and write the JSON. Caching within a run is allowed. Do not copy author-side reference data or hardcode answers. Importing `solution.py` must still have no workload side effects.

`results.schema.json` is the JSON Schema Draft 2020-12 contract. Every listed field is required and unlisted fields are rejected. Array order is arbitrary; record identity is determined by the keys below. Emit actual finite JSON numbers, not strings, NaNs, or infinities. The only numerical field allowed to be null is the explicitly undefined Chern number of a singular transition loop. Duplicate object keys or duplicate record identities are invalid. Use normal full-precision floating-point serialization, such as `json.dump(..., allow_nan=False)`, rather than rounded report values.

## Top-level fields

| Field | Contents |
| --- | --- |
| `schema_version` | The integer `2` |
| `input_sha256` | Lowercase SHA-256 hex digest of the exact bytes of the input JSON file |
| `task_sha256` | Lowercase SHA-256 hex digest of the exact bytes of the supplied `task.md` |
| `provenance` | Object containing nonempty `python`, `numpy`, `scipy`, and `command` strings, plus `solution_sha256`, the lowercase SHA-256 digest of the submitted `solution.py` bytes |
| `static` | 16 static records |
| `curvature` | 4 pointwise curvature records |
| `adiabatic_lines` | 4 fixed-twist adiabatic integral records |
| `dynamics` | 71 named trajectory records |
| `averages` | 6 twist-average records |
| `frame_static` | 4 flux-frame records, for `pump_L6` and `trivial_L6` at both meshes |
| `frame_dynamics` | 3 moving-frame records identified by baseline trajectory ID |
| `transition_diagnostics` | 3 records for the prescribed offsets of `free_L6` |

The provenance strings record the versions actually used and the command that generated the results. They are not numerical matching targets. If you edit `solution.py` after generation, regenerate the results to update its fingerprint. The task digest binds the extended specification as well as the unchanged baseline input data. The fingerprints identify supplied files; they do not prove that the calculations were executed.

## Static records

Each record contains `case`, `mesh`, `rephase`, and `values`. Its unique identity is `(case, mesh, rephase)`. Include every named case at each prescribed mesh, once with `rephase: false` and once with `rephase: true`. Compute both variants rather than simply duplicating the values.

`values` has exactly the five scalar keys returned by `static_observables`: `chern`, `sampled_min_gap`, `min_link_overlap`, `max_abs_plaquette_phase`, and `kubo_integral`.

## Curvature and adiabatic line records

Each `curvature` record contains `case`, `theta`, `phi`, and `value`. Include the pointwise curvature at each of the four public `operator_probes` coordinates. Each case occurs once. Coordinates must agree with the inputs to absolute tolerance `1e-12`.

Each `adiabatic_lines` record contains `case`, `phi`, and `value`. Include every named case at `dynamics.fixed_phi`. The value is the forward integral of the pointwise curvature over `theta` from zero to `2*pi`. Use the specified 128-node Gauss–Legendre rule on that interval for the recorded value, so it matches the evaluator's discrete definition. You may use additional quadrature refinements in your convergence analysis.

## Trajectory records and identifiers

Each `dynamics` record contains `id`, `case`, `phi`, `T`, `direction`, and `values`. The public validator prints the complete parameter-to-ID mapping without any computed answers:

```bash
python3 benchmark/validate_results.py --inputs benchmark/inputs.json --list-plan
```

The IDs are `trajectory_000` through `trajectory_070`. Their order follows the workload: the three forward fixed-twist pump periods; the reverse fixed-twist pump; the free, eight-site, and trivial fixed-twist controls; the forward twist calculations ordered by twist count, period, and midpoint index; and finally the reverse-pump and trivial-loop calculations interleaved at each of the eight midpoint indices. Use the printed plan to avoid ambiguities. The array itself need not be sorted. Each ID must occur exactly once and agree with its prescribed parameters.

`values` has exactly the six keys returned by `evolve`: `q_total`, `q_equilibrium`, `q_geometric`, `final_norm_squared`, `excitation_probability`, and `final_state`. The first five are finite real numbers. Serialize each complex amplitude as the two-element array `[real, imaginary]`, preserving the canonical basis order. Thus `final_state` is a list of 20 pairs for a six-site ring or 70 pairs for an eight-site ring. A complex number encoded as a string or as separate unlabelled real and imaginary arrays is not accepted.

The public checker verifies the charge subtraction, the state-vector length, its norm, the reported norm, and the probability range. The author-side evaluator additionally checks the excitation probability against the returned vector and compares the vector modulo global phase with the independent reference.

## Twist-average records

Each `averages` record contains `case`, `twists`, `T`, `direction`, `q_total`, `q_equilibrium`, and `q_geometric`. Its unique identity is `(case, twists, T, direction)`.

Include the four forward pump averages for `twists = 8, 16` and `T = 30, 90`, plus the reverse pump average and the forward trivial-loop average, both at eight twists and `T = 90`. Derive all three charges from the corresponding individual midpoint trajectories recorded in `dynamics`. The checker recomputes each arithmetic mean and the charge subtraction, requiring agreement within absolute tolerance `2e-10`.

Transfer probes used by the author-side evaluator do not belong in this file. Do not include unpublished trajectories or unknown extra records to fill the required counts.

## Validation and acceptance

Check your file with:

```bash
python3 benchmark/validate_results.py results.json \
  --inputs benchmark/inputs.json --solution solution.py
```

The checker reads the schema bundled beside itself, validates the hashes and complete record coverage, and checks internal numerical consistency. It returns exit status zero on success and one on failure, with a JSON diagnostic. It contains no numerical oracle; passing it does not establish correct physics.

Full author-side evaluation requires both the executable APIs and this artifact. The artifact must pass the public checks and agree with the stored numerical reference using the tolerances in `task.md`. This is an additional acceptance condition; the existing 100-point physics rubric remains unchanged. A missing, malformed, incomplete, inconsistent, or numerically incorrect artifact makes the full submission fail even if its API score is high.

The report's tables should come from `results.json`. A database or spreadsheet is optional and cannot replace this file. No plots, binary array formats, or database engine are needed for acceptance.

## Flux-frame records

Each `frame_static` record contains `case`, `mesh`, `ground_frames`, and `values`. Its unique identity is `(case, mesh)`. The cases are `pump_L6` and `trivial_L6`, each at meshes 12 and 24. These artifact records use `rephase=False`; the executable API is additionally checked with both values of the flag.

`ground_frames` is a two-dimensional array of node states: the element at `[a][b]` is a list of `binomial(L,L/2)` complex amplitudes in canonical coefficient order. Each amplitude is a `[real, imaginary]` pair. Thus its JSON shape is `(mesh, mesh, binomial(L,L/2), 2)`.

`values` contains exactly `sewn_chern`, `unsewn_chern`, `min_link_overlap`, and `max_abs_plaquette_phase`. The physical boundary identification, the diagnostic unsewn sum, and which links belong to the last two fields are defined in Section 9.1 of `task.md`. The checker verifies the arithmetic from the supplied states and, with the implementation available, checks their ground-state residuals and reevaluates the frame APIs.

## Moving-frame records

Each `frame_dynamics` record has only `trajectory_id` and `values`. The three IDs are `trajectory_000`, `trajectory_003`, and `trajectory_005`. Case, twist, period, and direction are inherited from the corresponding baseline record, avoiding a second parameter convention.

`values` contains exactly:

- `canonical_initial_state`, `canonical_final_state`, `frame_initial_state`, and `frame_final_state`: complex vectors represented as lists of `[real, imaginary]` pairs.
- `q_boundary`, `q_parameter`, and `polarization_change`: finite real numbers with the definitions in Section 9.2.
- `ramsey_raw` and `ramsey_covariant`: complex scalars, each encoded as one `[real, imaginary]` pair.

The canonical endpoint must represent the corresponding baseline final state modulo a common phase. That freedom does not permit independent phase alignment between canonical and frame states. Covariance and the complex readout values are checked with their relative phases intact. The frame charge is a total integrated current, not the baseline equilibrium-subtracted charge.

## Transition records

Each `transition_diagnostics` record contains `delta_offset`, `status`, `minimum_gap`, `witness`, `ground_rank`, `ground_projector`, `chern`, and `mesh_values`. The offsets `0.49`, `0.50`, and `0.51` each occur once. All other case parameters are those of `free_L6`.

`status` is `gapped` or `singular`. `witness` has finite real `theta` and `phi` coordinates. `ground_rank` is the integer rank of the lowest-energy eigenspace at that witness. `ground_projector` has shape `(20,20,2)` in JSON: every matrix element is a complex pair. It represents the full eigenspace, independently of the eigenvectors chosen within it.

For a gapped loop, `chern` is a finite real number and `mesh_values` has two records, each with `mesh`, `chern`, `kubo_integral`, and `sampled_min_gap`, at meshes 12 and 24. These are the unrounded baseline observables of the modified loop. For a singular rank-one bundle, `chern` is null and `mesh_values` is an empty list. Do not fill undefined quantities with zero or with finite values obtained by perturbing the specified loop.

## Version-2 checking and compatibility

The unchanged baseline input file still has its original data-schema version and SHA-256 fingerprint. Output version 2 adds requirements; it is not a retrospective change to version-1 acceptance. The previous output artifact is intentionally incomplete for this new task.

Without `--solution`, the public validator checks the schema, coverage, fingerprints, and array-level identities. With `--solution solution.py`, it also imports the candidate and checks the new APIs, ground-state residuals, moving-frame generator and current definitions, and transition projectors. It uses only the allowed dependencies and the submitted canonical model interface. The existing private evaluator still checks the original canonical numerical results against its independent oracle, so passing public identities alone is not full acceptance.

The original 100-point canonical score is unchanged. Version 2 also requires all extended artifact and API checks to pass. A high baseline score with a failed extension is a failed version-2 submission. The validator does not grade the correctness of a mathematical proof in the report; the continuous-gap and invariant-validity outputs have deterministic numerical and categorical checks.
