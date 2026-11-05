# Machine-readable submission results

Submit `results.json` alongside `solution.py` and `report.md`. Generate it from the same implementation used by the evaluator, with this terminal interface:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 solution.py \
  --inputs benchmark/inputs.json --output results.json
```

The command must execute the required calculations and write the JSON. Caching within a run is allowed. Do not copy author-side reference data or hardcode answers. Importing `solution.py` must still have no workload side effects.

`results.schema.json` is the JSON Schema Draft 2020-12 contract. Every listed field is required and unlisted fields are rejected. Array order is arbitrary; record identity is determined by the keys below. Emit actual finite JSON numbers, not strings, nulls, NaNs, or infinities. Duplicate object keys or duplicate record identities are invalid. Use normal full-precision floating-point serialization, such as `json.dump(..., allow_nan=False)`, rather than rounded report values.

## Top-level fields

| Field | Contents |
| --- | --- |
| `schema_version` | The integer `1` |
| `input_sha256` | Lowercase SHA-256 hex digest of the exact bytes of the input JSON file |
| `provenance` | Object containing nonempty `python`, `numpy`, `scipy`, and `command` strings, plus `solution_sha256`, the lowercase SHA-256 digest of the submitted `solution.py` bytes |
| `static` | 16 static records |
| `curvature` | 4 pointwise curvature records |
| `adiabatic_lines` | 4 fixed-twist adiabatic integral records |
| `dynamics` | 71 named trajectory records |
| `averages` | 6 twist-average records |

The provenance strings record the versions actually used and the command that generated the results. They are not numerical matching targets. If you edit `solution.py` after generation, regenerate the results to update its fingerprint. The fingerprints identify supplied files; they do not prove that the calculations were executed.

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
