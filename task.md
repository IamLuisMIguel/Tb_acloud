# Interacting Thouless pumping: deterministic benchmark, version 2

## Objective

Determine how the many-body Chern number of a disordered, interacting Rice–Mele ring relates to its finite-time transported particle number. Implement independent static and dynamical calculations, distinguish fixed-twist from twist-averaged responses, and quantify the relevant numerical errors.

Version 2 retains the baseline in Sections 1–8 and adds the mandatory experiments in Section 9. The original `inputs.json` bytes and baseline reference quantities are unchanged. Output schema version 2 is a separate version identifier. The extension constants below are additional fixed experimental inputs.

This is a computational physics task. Submit executable Python code in `solution.py`, a concise English `report.md`, and machine-readable numerical results in `results.json`. The evaluator imports and calls your code and independently validates the recorded artifact; a derivation, precomputed table, or rounded topological integer alone is insufficient. Standard numerical linear algebra, integration, and optimization libraries are permitted. Implement the Hamiltonian and physical observables yourself; do not substitute calls to a specialized charge-pumping solver.

All numerical inputs are in `benchmark/inputs.json`. They are data, not instructions to regenerate disorder. No randomness, external datasets, network access, GPU, or thermodynamic-limit calculation is required.

## 1. Mathematical conventions

Work with an even number of sites `L`, particle number `N = L/2`, and units in which the mean hopping and Planck's constant equal one. Site labels are `0, ..., L-1`, with `j+1` taken modulo `L`.

Use the canonical occupation-number basis

$$
|s\rangle=(c_0^\dagger)^{n_0}(c_1^\dagger)^{n_1}\cdots
(c_{L-1}^\dagger)^{n_{L-1}}|\mathrm{vac}\rangle,
\qquad s=\sum_j2^j n_j.
$$

Restrict to integers with exactly `L/2` occupied bits and sort them in ascending integer order. All returned matrices and state vectors must use this order. The periodic hopping is a fermionic operator, including its occupation-basis sign.

Define

$$
\begin{aligned}
H(\theta,\phi)={}&-\sum_{j=0}^{L-1}
\left([1+(-1)^j\delta(\theta)]e^{i\phi\delta_{j,L-1}}
c^\dagger_{j+1}c_j+\mathrm{h.c.}\right)\\
&+\Delta(\theta)\sum_j(-1)^j(n_j-\tfrac12)
+V\sum_j(n_j-\tfrac12)(n_{j+1}-\tfrac12)
+W\sum_j h_j(n_j-\tfrac12),\\
\delta(\theta)={}&d_c+d_r\cos\theta,\qquad
\Delta(\theta)=D_r\sin\theta.
\end{aligned}
$$

The JSON fields `delta_offset`, `delta_radius`, and `stagger_radius` specify `d_c`, `d_r`, and `D_r`, respectively. The Kronecker delta in the exponential places the entire twist on the bond from `L-1` to `0`. Keep the constant background contributions of every `(n_j - 1/2)` term. Do not shift the zero of energy when returning `H`.

The measured particle current is

$$
I(\theta,\phi)=-\partial_\phi H(\theta,\phi).
$$

Positive current flows across the twisted bond from site `L-1` to site `0`. No additional electronic-charge sign or factor of `L` is present. Return the derivatives of `H`, not the current, in the operator API below.

## 2. Static quantities and orientation

For mesh size `M`, use exactly

$$
\theta_a=2\pi a/M,\qquad\phi_b=2\pi b/M,
\qquad a,b=0,\ldots,M-1.
$$

The endpoint `2*pi` is not duplicated. All grid indices wrap periodically. Let `u_ab` be a normalized lowest-energy eigenstate. Define normalized overlaps

$$
U_{ab}=\frac{\langle u_{ab}|u_{a+1,b}\rangle}{|\langle u_{ab}|u_{a+1,b}\rangle|},\qquad
V_{ab}=\frac{\langle u_{ab}|u_{a,b+1}\rangle}{|\langle u_{ab}|u_{a,b+1}\rangle|}.
$$

Use plaquette phases and the lattice Chern number

$$
f_{ab}=\operatorname{Arg}(U_{ab}V_{a+1,b}U_{a,b+1}^{-1}V_{ab}^{-1}),
\qquad C_M=\frac{1}{2\pi}\sum_{ab}f_{ab}.
$$

Use the principal phase in `(-pi, pi]`; none of the prescribed plaquettes lies near its branch cut. Return the unrounded sum. Record the minimum magnitude of the unnormalized links, including links that cross the grid boundaries.

Also compute the pointwise spectral/projector Berry curvature in the same orientation,

$$
F_{\theta\phi}=-i\operatorname{Tr}\left(P[\partial_\theta P,\partial_\phi P]\right),
\qquad P=|u\rangle\langle u|.
$$

Derive and implement its excited-state spectral representation using Hamiltonian derivatives. This must be a pointwise linear-response calculation, rather than a reinterpretation of a finite plaquette. Excited-state degeneracies are allowed; the ground state is isolated in the prescribed cases.

Return the discrete quadrature

$$
C^{\mathrm{Kubo}}_M=\frac{2\pi}{M^2}\sum_{ab}F_{\theta\phi}(\theta_a,\phi_b).
$$

This quadrature is generally not an integer at finite mesh size. Do not round it or replace it by `C_M`.

For `rephase=True`, before constructing overlaps replace each eigenvector by

$$
u_{ab}\longmapsto
\exp\!\left(2\pi i\frac{(37a+17b+13ab)\bmod101}{101}\right)u_{ab}.
$$

All returned static observables must remain invariant within numerical tolerance. Report the minimum of `E_1 - E_0` over the specified grid as a **sampled minimum gap**, not as a certified gap bound over the continuous torus.

## 3. Finite-time protocol and observables

For a duration `T > 0` and `direction` equal to `+1` or `-1`, set

$$
s=t/T,\qquad p(s)=10s^3-15s^4+6s^5,\qquad
\theta(t)=\mathrm{direction}\,2\pi p(s),\qquad0\leq t\leq T.
$$

The twist is fixed during each trajectory. Start in the normalized instantaneous ground state at `theta = 0`, at that trajectory's own twist. Solve

$$
i\partial_t|\psi(t)\rangle=H(\theta(t),\phi)|\psi(t)\rangle.
$$

Use the actual evolving state, without normalizing it inside the current expectation or concealing accumulated norm error. Integrate

$$
\begin{aligned}
Q_{\mathrm{total}}&=\int_0^T\langle\psi(t)|I(\theta(t),\phi)|\psi(t)\rangle\,dt,\\
Q_{\mathrm{equilibrium}}&=\int_0^T\langle u_0(\theta(t),\phi)|I(\theta(t),\phi)|u_0(\theta(t),\phi)\rangle\,dt,\\
Q_{\mathrm{geometric}}&=Q_{\mathrm{total}}-Q_{\mathrm{equilibrium}}.
\end{aligned}
$$

Here `u_0` is a normalized instantaneous ground state. The label `geometric` denotes this subtraction at finite `T`; it does not assert that all nonadiabatic corrections have vanished. Persistent current is generally nonzero at the prescribed fixed twist.

Return the final state in the canonical basis, its norm squared, and

$$
P_{\mathrm{exc}}=1-
\frac{|\langle u_0(0,\phi)|\psi(T)\rangle|^2}{\langle\psi(T)|\psi(T)\rangle}.
$$

The initial Hamiltonian equals the final Hamiltonian for both directions. Overall state phase is irrelevant and will be removed by the evaluator's fidelity comparison.

For twist averages, use exactly the midpoint nodes

$$
\phi_b=2\pi(b+\tfrac12)/K,\qquad
\overline Q_{K,T}=K^{-1}\sum_{b=0}^{K-1}Q_{\mathrm{geometric}}(\phi_b,T).
$$

Do not substitute a current averaged over physical bonds for an average over independently initialized twists.

## 4. Required calculations

Use the exact cases and workload settings in `benchmark/inputs.json`:

- Static calculations at `M = 12, 24` for `free_L6`, `pump_L6`, `pump_L8`, and `trivial_L6`, both with and without the specified rephasing.
- Pointwise curvature and operator evaluations at the listed probe coordinates. The evaluator also uses deterministic transfer probes in the same model, including `L = 4`.
- Forward trajectories for `pump_L6` at the fixed twist and `T = 10, 30, 90`, and its reverse trajectory at `T = 90`.
- The fixed-twist free control at `T = 30`, the interacting `L = 8` trajectory at `T = 30`, and the trivial-loop control at `T = 90`.
- Forward `pump_L6` twist averages with `K = 8, 16` at `T = 30, 90`.
- Reverse-pump and forward trivial-loop twist averages with `K = 8` at `T = 90`.

These comprise 71 distinct named trajectories. The evaluator may additionally request a small transfer trajectory with `L = 4` and duration between `10` and `30`. Your functions must operate on the supplied parameters, rather than select tabulated answers by case name. Transfer probes use the same definitions with `0 <= V <= 1`, `0 <= W <= 0.3`, disorder entries in `[-1,1]`, and loop parameters within the ranges spanned by the supplied cases and `transfer_probe`.

In the report:

1. Explain the fermionic boundary sign, the current convention, and the relation of the chosen Berry-curvature orientation to adiabatic transport.
2. Compare the two static meshes, including their unrounded Kubo quadratures and overlap diagnostics.
3. Evaluate the fixed-twist adiabatic response `integral_0^(2*pi) F(theta, phi) dtheta` and compare it with the finite-time fixed-twist results. The evaluator independently integrates your pointwise curvature using a fixed 128-node Gauss–Legendre rule.
4. Tabulate total, equilibrium, and subtracted charges. Compare twist-quadrature refinement and longer periods separately. Discuss the reversal and trivial-loop controls.
5. Distinguish finite-size, finite-period, time-integration, and twist-quadrature effects. Do not claim exact integer transport at finite `T`, exact antisymmetry of every finite-time trajectory, or monotonic error reduction without evidence.

No extrapolation to an infinite chain is required. Comparing `L = 6` and `8` is a finite-size diagnostic, not proof of a thermodynamic limit. You may use any converged evolution or quadrature method; the evaluation does not require the reference implementation's algorithm.

## 5. Python interface

Implement the following functions. `case` is a dictionary from `inputs.json` containing only `L`, `V`, `W`, `disorder`, and `loop`; it does not contain its case name. Do not mutate caller-owned inputs.

```python
def build_operators(case, theta, phi):
    return {
        "H": H,
        "dtheta": dH_dtheta,
        "dphi": dH_dphi,
    }

def kubo_curvature(case, theta, phi):
    return real_curvature

def static_observables(case, mesh, rephase=False):
    return {
        "chern": unrounded_chern,
        "sampled_min_gap": minimum_gap_on_grid,
        "min_link_overlap": minimum_unnormalized_link_magnitude,
        "max_abs_plaquette_phase": largest_absolute_principal_phase,
        "kubo_integral": unrounded_kubo_quadrature,
    }

def evolve(case, phi, T, direction, rtol=1e-10, atol=1e-12):
    return {
        "q_total": total_charge,
        "q_equilibrium": equilibrium_charge,
        "q_geometric": total_charge - equilibrium_charge,
        "final_norm_squared": final_norm_squared,
        "excitation_probability": excitation_probability,
        "final_state": complex_final_state,
    }
```

Matrices must be dense two-dimensional NumPy-compatible arrays of dimension `binomial(L, L/2)`. The state must be a one-dimensional complex NumPy-compatible array of that length. Scalar fields must be finite real numbers. The functions must support arbitrary real angles; the operator and curvature calls must respect their periodicity.

The `rtol` and `atol` arguments specify requested integration accuracy. Methods that do not use adaptive ODE tolerances may implement equivalent internal refinement. In all cases the returned observables must satisfy the acceptance tolerances below. Additional Python modules are allowed, provided `solution.py` can import them from its own directory. Importing your module must not launch the complete workload.

## 6. Deterministic numerical evaluation

The evaluator uses fixed inputs and fixed transfer probes, stored author-side reference observables, and independently checked operator formulas. It never generates a fresh random seed during grading. Results need not be bitwise identical across BLAS implementations.

| Quantity | Acceptance criterion |
| --- | --- |
| Every matrix entry | `abs(error) <= 2e-11 + 1e-11 * abs(reference)` |
| Every returned static scalar | Absolute error at most `2e-8` |
| Pointwise curvature | Absolute error at most `2e-8` |
| Fixed-twist adiabatic curvature integral | Absolute error at most `2e-7` under the evaluator's 128-node rule |
| Each signed charge | Absolute error at most `3e-6` against the finite-time reference |
| `q_geometric = q_total - q_equilibrium` | Absolute consistency error at most `2e-10` |
| Excitation probability | Absolute error at most `3e-7` |
| Final norm squared | Absolute deviation from one at most `1e-7` |
| Reported norm versus returned state | Absolute consistency error at most `2e-10` |
| Final state | Phase-independent infidelity at most `1e-8` |

The probability is also reconstructed from the returned final state. NaNs, infinities, wrong shapes, and missing fields fail the corresponding check. Final-state infidelity uses normalized overlap, while norm conservation is checked separately.

The score is divided into operator construction (25%), static topology and gauge invariance (25%), pointwise and integrated linear response (10%), and finite-time evolution (40%). Within each group, equally weighted evaluator checks receive credit only when all quantities in that check pass. A successful numerical submission must pass every check; partial scores help diagnose failures.

The English report is required as an explanation, but prose quality is not included in the deterministic numerical score. Mathematical correctness beyond the specified numerical tests is not mechanically certified. Do not replace numerical evidence with claims about convergence.

## 7. Execution and reproduction

The calibrated environment uses Python 3.12, NumPy 2.3.5, and SciPy 1.17.0. The artifact validator also uses jsonschema 4.26.0. These versions are listed in `benchmark/requirements.txt`. The workload is CPU-only. Use `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1` when reproducing reference timings or comparing refinements.

In `report.md`, give terminal commands that import your implementation, execute the required calculations, and reproduce your tables. Generate numerical results from your functions rather than embedding reference tables. A plot is optional and receives no numerical score.

No reference answers, author implementation, or grading code are part of the candidate inputs.

## 8. Machine-readable artifact

Implement the terminal command `python3 solution.py --inputs benchmark/inputs.json --output results.json`. It must run the workload using your implementation and write the output contract specified in `benchmark/results_format.md` and `benchmark/results.schema.json`. Keep module imports free of workload side effects.

Record all 16 static/rephasing calculations, four pointwise curvatures, four fixed-twist adiabatic integrals, all 71 named trajectories, and six twist averages, together with version metadata and SHA-256 fingerprints of the input file and submitted `solution.py`. Serialize complex state amplitudes as `[real, imaginary]` pairs and preserve full numerical precision. Do not add matrices or private transfer probes to the artifact.

The public command `python3 benchmark/validate_results.py --list-plan` prints the required trajectory IDs and coordinates. Validate your completed artifact with `python3 benchmark/validate_results.py results.json --inputs benchmark/inputs.json --solution solution.py`. This checks structure, coverage, fingerprints, and internal consistency without revealing reference values. The author-side evaluator also compares recorded results to the numerical oracle, in addition to exercising your APIs. Both must pass; the existing numerical score remains out of 100. Use the JSON as the source of the report's numerical tables. A database may be supplied as an optional convenience but does not replace `results.json`.

## 9. Version-2 scientific extension

All three experiments below are mandatory. Numerical methods are your choice. The additional arrays are evidence to be evaluated, rather than a prescribed implementation strategy. The original 100-point API score remains diagnostic; complete acceptance also requires all version-2 checks to pass. Passing only the baseline is insufficient for this revised task.

### 9.1 Flux frame and the parameter-space boundary

Define the dimensionless many-body operators

$$
X=\sum_{j=0}^{L-1}\frac{j}{L}n_j,
\qquad R(\phi)=e^{i\phi X},
\qquad H_R(\theta,\phi)=R(\phi)H(\theta,\phi)R(\phi)^\dagger.
$$

The coefficients are represented in the same canonical occupation basis as the baseline. This is a different state frame for the same physical family. Its fiber identification at the twist boundary is

$$
v(\theta,\phi+2\pi)=S\,v(\theta,\phi),\qquad S=e^{2\pi iX}.
$$

For `pump_L6` and `trivial_L6`, at both baseline meshes, the required evidence is a normalized ground-state frame at every mesh node and two unrounded lattice Chern sums. `sewn_chern` uses the physical boundary identification above and the plaquette orientation of Section 2. `unsewn_chern` is the diagnostic sum obtained by treating the coefficient arrays as ordinarily periodic at the twist boundary. The reported minimum link magnitude and largest plaquette phase refer to the physically identified sum. State which result is an invariant of the specified family and explain any discrepancy. Include the explicitly required rephasing behavior in the API.

### 9.2 Moving frame, current response, and phase-sensitive readout

For each of baseline trajectories `trajectory_000`, `trajectory_003`, and `trajectory_005`, let `s=t/T` and define

$$
a_j=\frac{(-1)^j(j+1)}{5L},\quad
A=\sum_j a_jn_j,\quad\widehat N=\sum_jn_j,\quad
G(s,\phi)=\exp\!\left[i\left(\phi X+\sin(2\pi s)A+\frac{\phi s}{4}\widehat N\right)\right].
$$

The canonical state follows the original Schrödinger equation and protocol. The frame state is defined by `psi_frame(t) = G(s, phi) psi_canonical(t)`, including its relative phase and its transformed initial condition. The frame Hamiltonian is the generator for which this state satisfies the Schrödinger equation. Determine that generator, its total partial derivative with respect to `phi` at fixed `theta`, `s`, and `T`, and the physical boundary-current operator expressed in this frame. The latter measures particles crossing the same physical bond as the baseline.

The required integrated observables are

$$
Q_{\rm boundary}=\int_0^T\langle\psi_f|I_{{\rm boundary},f}|\psi_f\rangle\,dt,
\qquad
Q_{\rm parameter}=-\int_0^T\langle\psi_f|\partial_\phi H_f|\psi_f\rangle\,dt.
$$

These are total currents, without equilibrium-current subtraction. Also report the endpoint change in the canonical expectation of

$$
K(s,\phi)=-iG(s,\phi)^\dagger\partial_\phi G(s,\phi).
$$

The two complex readout amplitudes are defined in a fixed coefficient basis and with a transformed detector, respectively:

$$
\mathcal A_{\rm raw}=\langle\psi_c(T)|\psi_f(T)\rangle,
\qquad
\mathcal A_{\rm cov}=\langle\psi_c(T)|G(1,\phi)^\dagger|\psi_f(T)\rangle.
$$

They are Ramsey overlap diagnostics with the stated readout conventions. The raw amplitude is deliberately frame dependent; it is not claimed to be a gauge-invariant charge. Return canonical and frame states at the beginning and end, the two charges, the change in `K`, and both complex amplitudes. An arbitrary common initial phase is permitted, but the frame and canonical branches must use their defined relative phase. Independent endpoint phase alignment, or dropping a scalar part of the generator, is incompatible with these readouts.

Explain the relationship between the two current responses and the endpoint quantity, and which current corresponds to the baseline particle detector. Do not identify them merely because both involve a flux parameter.

### 9.3 Validity of the ground-state bundle at a transition

Use the `free_L6` case, retaining its original `delta_radius = 0.5` and `stagger_radius = 1`, but set `delta_offset` in turn to `0.49`, `0.50`, and `0.51`. All other parameters are unchanged.

For each loop report whether the rank-one ground-state bundle is gapped or singular over the continuous `(theta, phi)` torus. Return the continuous minimum of `E1-E0`, a minimizing coordinate pair, the ground-space rank and full ground-space projector at that pair, and the rank-one Chern number if it is defined. A singular bundle must have `chern: null`, rather than a regularized or arbitrarily assigned integer. The projector at a degeneracy means the entire lowest-energy eigenspace, not a selected eigenvector from it.

For each gapped loop also return the unrounded baseline static quantities `chern`, `kubo_integral`, and `sampled_min_gap` at meshes 12 and 24. For a singular loop the `mesh_values` array is empty. The report must justify the continuous-gap and invariant-validity decisions and explain what the two finite-mesh response quadratures do or do not establish. Do not infer a continuous gap solely from a positive value on an unrelated sample grid.

### 9.4 Additional executable interfaces

The following names and return contracts are evaluation interfaces. They do not constrain the internal numerical algorithms.

- `frame_operators(case, theta, phi, s, T)` returns `H`, `dphi`, and `boundary_current` as dense complex matrices in the canonical coefficient basis. The frame is `G` from Section 9.2, not just the static `R` from Section 9.1.
- `frame_static(case, mesh, rephase=False)` returns `ground_frames`, an array of shape `(mesh, mesh, binomial(L, L/2))`, and `values` with keys `sewn_chern`, `unsewn_chern`, `min_link_overlap`, and `max_abs_plaquette_phase`. This experiment uses the static `H_R`. When requested, `rephase=True` applies the same node phases as Section 2 before the observables are reported.
- `frame_evolve(case, phi, T, direction, rtol=1e-10, atol=1e-12)` returns `canonical_initial_state`, `canonical_final_state`, `frame_initial_state`, `frame_final_state`, `q_boundary`, `q_parameter`, `polarization_change`, `ramsey_raw`, and `ramsey_covariant`. The four states are complex vectors, both Ramsey values are complex scalars, and the remaining values are real. The canonical endpoint must agree with the corresponding baseline physical trajectory modulo a common phase; the frame-covariance and readout comparisons do not allow an independent frame phase.
- `diagnose_transition(case, delta_offset)` returns `status`, `minimum_gap`, `witness`, `ground_rank`, `ground_projector`, `chern`, and `mesh_values` as defined in Section 9.3 and the JSON contract. `case` is the unmodified `free_L6` dictionary; the argument supplies its replacement offset. The API also supports offsets in `[0.48,0.495]`, exactly `0.5`, and `[0.505,0.52]` for deterministic transfer checks.

The frame-operator API also supports the supplied four-site transfer model and arbitrary finite angles in the same moderate range as the baseline. No function may mutate the caller's case dictionary.

### 9.5 Additional acceptance criteria

Frame Hamiltonian covariance, Hermiticity, and the boundary-current matrix use absolute tolerance `2e-10`. The total flux derivative is checked to `2e-8`. Frame-grid ground-state residuals and componentwise covariance of endpoint vectors use `2e-7`; vector normalization uses `1e-7`. Relative phases in these covariance checks are retained. Complex Ramsey amplitudes use absolute tolerance `3e-7`, and integrated charges and the endpoint-current relation use `3e-6`.

New Chern sums and finite-mesh scalars use `2e-8`. Transition minimum gaps use `2e-8`; minimizing angles are compared modulo `2*pi` within `1e-6`. Ground projectors use componentwise tolerance `2e-7`, including Hermiticity, idempotence, rank, and agreement with the entire lowest-energy subspace. The singular/gapped classification and null-invariant decision are exact categorical outputs.

These are new requirements for version 2, not grounds for reclassifying a successful version-1 submission. All new experiments retain the small CPU-feasible Hilbert spaces; no larger-system or unreasonable precision requirement is being introduced.
