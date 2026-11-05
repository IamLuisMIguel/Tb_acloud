# Author-side specification and verification

## What has been established

The selected idea is now a fixed, numerically graded benchmark. The English natural-language prompt to give the candidate is `instructions.md`. Its binding mathematical specification is in `task.md`; literal model inputs, workloads, and tolerances are in `benchmark/inputs.json`. The candidate must also generate `results.json` under the public schema and format contract. No disorder is generated during evaluation.

The benchmark evaluates a finite, precisely specified computational problem. It does not ask a model to establish an infinite-system phase diagram, prove continuous-torus gap bounds, or obtain exact finite-time quantization. Those would introduce claims that the proposed finite workload cannot certify.

All mandatory calculations have been executed in the current environment. The stored reference covers four cases, two static meshes, four pointwise reference coordinates, four fixed-twist adiabatic line integrals, and 71 named trajectories. Additional deterministic `L = 4` transfer checks are calculated at grading time.

## Distribution: keep the answers out of the task

The candidate distribution consists of exactly these files:

```text
instructions.md
task.md
benchmark/inputs.json
benchmark/requirements.txt
benchmark/results.schema.json
benchmark/results_format.md
benchmark/validate_results.py
```

Do not give the candidate this directory, `idea.md`, a full repository checkout containing the author files, calibration logs, or the stored reference results. The directory name `author_only` is organizational; it is not an access-control mechanism.

Create an archive with an explicit file allowlist:

```bash
python3 author_only/package_candidate.py --output /tmp/tb_acloud_candidate.tar.gz
```

The archive has fixed metadata, file order, and gzip timestamp. Repeated packaging of unchanged public files produces the same SHA-256 hash. The archive excludes the solutions, evaluation source, and reference observables. This provides a concrete handoff artifact; it does not itself create an isolated evaluation service.

## Author files

| File | Purpose |
| --- | --- |
| `reference.py` | Executable operator, spectral, topological, and dynamical oracle |
| `reference_results.json` | Tight-tolerance finite-time reference values and state vectors |
| `grade.py` | Functional evaluator, consistency checks, and weighted numerical score |
| `validate_reference.py` | Independent identities, convergence checks, and synthetic fault rejection |
| `calibration.json` | Recorded verification outcomes and numerical error margins |
| `package_candidate.py` | Deterministic public-only archive generator |
| `export_reference_artifact.py` | Converts the stored oracle to the candidate JSON format for author-side testing |
| `test_results.py` | Artifact regression tests, including missing data and numerically incorrect but internally consistent results |

## Reproduce the reference and calibration

Use the Python environment described in `benchmark/requirements.txt`. From the repository root:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 author_only/reference.py \
  --output author_only/reference_results.json --tight

OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 author_only/reference.py \
  --output /tmp/pump-reference-baseline.json

OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 author_only/validate_reference.py \
  --baseline /tmp/pump-reference-baseline.json \
  --output author_only/calibration.json
```

The reference fixture includes the SHA-256 digest of the exact bytes of `inputs.json`. The grader rejects a fixture generated for different inputs. Changing whitespace in that JSON also changes the digest; regenerate the reference after an intentional input change.

The baseline uses DOP853 with `rtol = 1e-10`, `atol = 1e-12`, and 128-point Gauss–Legendre equilibrium-current integration. The tight run uses `rtol = 2e-12`, `atol = 2e-14`, and 256 points. The candidate is free to use a different converged algorithm.

Observed reference-generation times for the initial static-and-dynamical calibration were approximately 109 seconds for the tight run and 40 seconds for the baseline on this cloud machine, with one BLAS thread. Subsequent software revisions added the fixed-twist curvature line integrals. These timings are feasibility evidence, not a portable runtime guarantee or a grading cutoff. Set an execution budget after timing the actual evaluation hardware; do not score wall-clock speed as physics correctness.

## Independent checks already run

The calibration script verifies:

- Full-Fock Jordan–Wigner tensor-product construction against the restricted-basis operator implementation, including the `L = 8` boundary sign.
- Hamiltonian derivatives against centered finite differences.
- Pointwise spectral curvature against finite-difference derivatives of the ground-state projector, including its sign.
- Free many-body ground-state energy against occupied single-particle levels.
- The specified deterministic eigenvector rephasing and reversed-loop Chern orientation.
- All 71 trajectories under tighter ODE tolerances and doubled equilibrium-current quadrature.
- One independent exactly unitary midpoint propagator at 2,000 and 4,000 steps, with Richardson refinement of its current integral.
- Fixed-twist adiabatic line integrals at 128 and 256 quadrature nodes.
- Rejection of synthetic incorrect results: omitted interaction, reversed current derivative, omitted persistent-current subtraction, reversed geometric-charge sign, and transport replaced by an integer.

The largest charge difference between baseline and tight trajectories was approximately `1.24e-9`, compared with a grading tolerance of `3e-6`. The independent propagator's extrapolated charge differed by approximately `1.03e-11` in the tested case. These are empirical checks, not mathematically certified error bounds. The complete magnitudes and environment versions are recorded in `calibration.json`.

## Grade a candidate

The numerical evaluator imports the candidate's `solution.py` and calls the four documented APIs. From the repository root:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 author_only/grade.py \
  --submission /path/to/submission/solution.py \
  --results /path/to/submission/results.json \
  --output /tmp/grade.json
```

The JSON output records individual pass/fail results, group scores, artifact validation, and the overall status. Exit status is zero only if all selected checks and the required artifact checks pass. In a full run, `--results` defaults to `results.json` beside the submission file. Optional `--group operators`, `--group static`, `--group curvature`, or `--group dynamics` runs the corresponding subset; its score is out of that subset's available points, not out of 100. Subset runs check the artifact only when `--results` is explicitly supplied.

An end-to-end smoke test uses the reference itself as the submission:

```bash
python3 author_only/export_reference_artifact.py --output /tmp/reference-artifact.json

OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 author_only/grade.py \
  --submission author_only/reference.py \
  --results /tmp/reference-artifact.json \
  --output /tmp/reference-grade.json
```

This tests evaluator orchestration and data handling. It is not an independent physical validation; the tensor-product, projector, free-limit, and unitary-propagation checks provide that additional evidence. Synthetic fault probes confirm specific rejection paths, but do not substitute for trials against real model-generated submissions.

The executed end-to-end smoke test passed all 106 checks: 7 operator probes, 16 static/gauge checks, 11 curvature checks, and 72 trajectories, including the transfer trajectory. Its numerical score was 100/100.

The JSON extension preserves these API checks and adds a mandatory artifact gate. It checks all 95 recorded numerical results (16 static variants, 4 point curvatures, 4 adiabatic lines, and 71 trajectories) against the oracle, plus complete coverage, source/input fingerprints, and six self-consistent twist averages. A failed artifact does not alter the diagnostic 100-point API rubric but makes the overall submission fail.

The score is 25 points for operator probes, 25 for static/gauge checks, 10 for pointwise and integrated curvature, and 40 for trajectories. Each evaluator check is all-or-nothing across its required fields; checks have equal weight within their group. Scores are diagnostic, and passing requires every check. Explanatory prose is not part of the automated numerical score. This evaluator does not enforce report existence, report language, or correctness of derivations through text matching.

The reference and candidate state vectors are compared modulo global phase, and normalization is checked separately. Signed charges, excitation probabilities, and the consistency of returned fields are also tested. No returned integer is accepted as a substitute for finite-time evolution.

## Machine-readable results and regression checks

`benchmark/results_format.md` defines the candidate's `results.json`, and `benchmark/results.schema.json` is its strict Draft 2020-12 schema. The public validator checks parsing, finite values, record counts and identities, coordinates, source/input fingerprints, final-state dimensions and norms, gauge consistency, charge subtraction, and arithmetic twist averages. It contains no oracle values. The author-side gate then compares the recorded numbers and phase-invariant final states with the existing reference.

The candidate must implement `python3 solution.py --inputs benchmark/inputs.json --output results.json` and generate the report tables from that file. Metadata is provenance information, not proof of execution. API reevaluation and reference comparisons remain necessary. The grader does not run that CLI itself or verify report prose; it validates the submitted artifact and independently calls the APIs.

An author-side fixture can be exported without rerunning the simulation, using the command above. This conversion reuses the stored expected static results for both rephasing variants; it is an adapter for grader tests, not evidence of a new rephasing calculation. It is excluded from the candidate archive, as are the answers it exports.

Run the artifact regression suite with:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 author_only/test_results.py
```

The tests cover successful format and numerical validation; record-order and global-phase invariance; agreement between the public trajectory plan and the reference schedule; missing artifacts; 17 malformed or inconsistent data variants; strict JSON parsing; and an internally consistent but physically incorrect charge that the public checker accepts and the numerical oracle rejects. The validator dependency `jsonschema==4.26.0` is included in the public requirements. Core input bytes and physics tolerances are unchanged, so the original reference fixture remains applicable.

## Calibrated scientific distinctions

For the interacting six-site case at the prescribed fixed twist:

| Quantity | Value |
| --- | --- |
| Static Chern number | `+1` |
| Fixed-twist adiabatic curvature integral | approximately `0.9165221642` |
| Forward geometric charge at `T = 10` | approximately `1.0274864226` |
| Forward geometric charge at `T = 30` | approximately `0.9513765777` |
| Forward geometric charge at `T = 90` | approximately `0.9266761866` |
| Total charge at `T = 90` | approximately `-0.4436405368` |
| Equilibrium-current integral at `T = 90` | approximately `-1.3703167234` |
| Forward geometric charge averaged over 16 twists, `T = 90` | approximately `1.0000004359` |

These values make several tempting shortcuts fail numerically. At fixed twist, even the adiabatic charge need not equal the Chern number. Persistent current can change the sign of total charge. Very small final excitation probability does not certify the integrated charge to the same accuracy. Finite-time reversal need not give exact antisymmetry, and a trivial loop may have nonzero subtracted charge at an individual twist.

The refined Kubo mesh integral is also not forced to equal an integer: for `trivial_L6` at `M = 24` it is approximately `-1.07627e-7`, whereas the overlap-based Chern number is zero within roundoff. Keeping both quantities prevents a solver from disguising a quadrature result by rounding.

## Scope of determinism and remaining limitations

The workload, probes, input bytes, sign conventions, quadrature definitions, tolerances, score, and pass/fail rules are fixed. Evaluation uses no random sampling or subjective interpretation of plots. Different BLAS or eigensolver implementations are accommodated by tolerances and phase-invariant state comparisons.

This is numerical verification over a defined test set, not a formal proof of arbitrary-parameter correctness or of how a candidate obtained its answers. Functional transfer probes reduce the value of tabulating public examples, but cannot prove that an implementation contains no special cases. Independent checks reduce reference risk; they do not eliminate all correlated implementation errors.

The current calibration uses one software stack and one machine. Test a genuinely independent candidate implementation and another supported numerical stack before calling the benchmark portable or finalized. No real frontier-model evaluation has yet established its difficulty or discrimination rate. The rubric and public task are a calibrated version 1, ready for that pilot stage.
