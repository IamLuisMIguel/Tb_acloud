# PhD-level computational physics tasks for a frontier LLM

**Selected direction:** Idea 1 has been developed into [a deterministic task specification](task.md), with fixed [inputs](benchmark/inputs.json). The [author-side verification notes](author_only/README.md) describe the oracle, calibration, and candidate packaging. Author-side material must not be included in the candidate workspace.

## Purpose and recommendation

This document proposes six coherent research-style tasks that require actual computation through a terminal. Each combines a central physics question with numerical implementation, verification, and interpretation. The objective is to expose scientifically consequential mistakes, rather than reward long answers, obscure terminology, or expensive calculations.

**Recommended main task: interacting topological charge pumping.** It connects a many-body topological invariant to real-time transport, with several independent checks and a particularly subtle distinction between fixed-flux and flux-averaged observables.

**Recommended lower-cost alternative: the non-Hermitian skin effect and spectral instability.** An exact similarity transformation provides a strong grading oracle, while a conventional eigensolver can produce a convincing but physically incorrect answer. A small numerical experiment below already demonstrates this failure mode.

These are task-design proposals, not completed benchmark problems. Their ability to discriminate between frontier models must be measured in pilot evaluations; difficulty cannot be established from topic names alone. Only the exploratory calculations explicitly reported below have been executed.

## Comparison

| Idea | Central question | Main source of difficulty | Strongest verification route | Relative computational burden |
| --- | --- | --- | --- | --- |
| 1. Interacting charge pump | When does quantized topology predict measured transport? | Gauge invariance, fermionic signs, adiabaticity, finite-size effects | Independent Chern and dynamical calculations | Medium |
| 2. Non-Hermitian chain | Which apparent spectral features are physical? | Nonnormality, boundary sensitivity, unstable eigenvalues | Exact similarity to a Hermitian matrix | Low to medium |
| 3. Efimov renormalization | Can cutoff-independent three-body predictions emerge from a running coupling? | Integral equations, root tracking, regulator dependence | Discrete scale invariance and cross-regulator agreement | Medium |
| 4. Near-extremal Kerr modes | Can the same physical mode be tracked toward extremality? | Complex roots, coupled angular/radial problems, singular limits | Independent truncation and precision checks | Medium to high |
| 5. BKT transition | Can a transition be located with credible uncertainty? | Autocorrelation and logarithmic finite-size corrections | Multiple observables and statistical calibration | Medium to high |
| 6. Analytic continuation | What real-frequency information is actually identifiable? | Ill-posed inversion and correlated noise | Forward checks and bounds on spectral observables | Low to medium |

Burden is a qualitative comparison, not a measured runtime guarantee. All six admit CPU-only versions; a GPU or cluster should not be a hidden prerequisite. Before fixing the final task, benchmark its required instances on the actual evaluation hardware.

## 1. Interacting Thouless pumping: topology versus finite-time transport

### Physical setup

Use a half-filled ring of spinless fermions with an interacting Rice–Mele Hamiltonian. Set the mean hopping and Planck's constant to one. For even length `L`, define

$$
H(\theta,\phi)=
-\sum_{j=0}^{L-1}\left[t_j(\theta)e^{i\phi\delta_{j,L-1}}
c^\dagger_{j+1}c_j+\mathrm{h.c.}\right]
+\Delta(\theta)\sum_j(-1)^j(n_j-\tfrac12)
+V\sum_j(n_j-\tfrac12)(n_{j+1}-\tfrac12)
+W\sum_j h_j(n_j-\tfrac12),
$$

with site indices taken modulo `L`,

$$
t_j(\theta)=1+(-1)^j\delta_0\cos\theta,
\qquad \Delta(\theta)=\Delta_0\sin\theta.
$$

Here `theta` parametrizes a pumping cycle and `phi` is a boundary twist. A useful starting point is `delta_0 = 0.5`, `Delta_0 = 1`, `V = 1`, and `W = 0.3`. Generate the fixed disorder vector with `numpy.random.default_rng(42).uniform(-1, 1, L)` and store the actual vector in the benchmark inputs. The noninteracting clean case is a required control.

### Coherent sequence of requests

1. Construct the number-conserving many-body Hamiltonian, including the fermionic sign of the boundary hopping. Verify Hermiticity and compare the noninteracting many-body energy against occupied single-particle levels.
2. Compute the lowest two energies over the `(theta, phi)` torus. Report the sampled minimum gap and explain why a sampled minimum is not a rigorous bound over the continuous torus.
3. Calculate the ground-state Chern number using normalized overlaps around mesh plaquettes. Repeat with finer meshes and independent random phases attached to every eigenvector. Report the smallest link overlap and the largest plaquette phase, rather than only a rounded integer.
4. Independently evolve the initial ground state through `theta(t) = 2*pi*t/T`. Derive the current from the implemented twist convention and integrate it. Separate transported charge from the instantaneous equilibrium persistent-current contribution.
5. Compare transport at a fixed twist with transport averaged over twists. Investigate time-step, cycle-duration, twist-quadrature, and system-size errors separately. Explain which comparison is protected by the many-body Chern number and which is only approximately quantized.
6. Reverse the pumping cycle and show the corresponding reversal of the geometric response. Test a contractible loop that remains in a gapped region as a trivial control.

For example, define overlap links

$$
U_\mu(x)=\frac{\langle\psi(x)|\psi(x+\hat\mu)\rangle}
{|\langle\psi(x)|\psi(x+\hat\mu)\rangle|}
$$

and sum the principal phases of

$$
U_\theta(x)U_\phi(x+\hat\theta)
U_\theta(x+\hat\phi)^{-1}U_\phi(x)^{-1}.
$$

The solver must state the orientation and derive the sign relating this definition to its current convention. A sign mismatch should be diagnosed, not silently removed by taking absolute values.

### Why this can challenge a frontier model

A correct-looking implementation can still use the wrong boundary sign, differentiate eigenvectors in an arbitrary gauge, identify a coarse-mesh integer with convergence, or equate finite-time charge at one twist with a topological invariant. Interactions also invalidate shortcuts that replace the many-body state by an occupied-band projector.

The difficult part is reconciling independent observables and controlled limits. A memorized statement that “the pumped charge equals the Chern number” is insufficient.

### Verification and feasible scope

Use `L = 6, 8` as mandatory small instances; consider `L = 10, 12` only after calibration. The half-filled Hilbert spaces have dimensions `20`, `70`, `252`, and `924`, respectively. Start static calculations with `12 x 12` and `24 x 24` torus meshes. The final dynamics workload and tolerances must be chosen after timing a reference implementation.

Useful checks include energy agreement in the free limit, gauge invariance, mesh convergence, norm conservation, an independently differentiated current, loop reversal, and agreement between flux-averaged adiabatic transport and the appropriately signed Chern number. A near-closing gap must trigger a qualified conclusion, rather than an unconditional topological claim.

### Exploratory calculation performed for this proposal

A direct fixed-particle-number implementation using NumPy and SciPy produced the following static results. The reported sign follows the overlap-plaquette orientation written above.

| `L` | `V` | `W` | Mesh | Chern number | Sampled minimum gap | Smallest link-overlap magnitude |
| --- | --- | --- | --- | --- | --- | --- |
| 6 | 0 | 0 | `12 x 12` | 1.000000 | 2.000000 | 0.947730 |
| 6 | 1 | 0.3 | `12 x 12` | 1.000000 | 2.002854 | 0.919921 |
| 6 | 1 | 0.3 | `24 x 24` | 1.000000 | 2.002854 | 0.977601 |
| 8 | 1 | 0.3 | `12 x 12` | 1.000000 | 2.067033 | 0.899430 |
| 8 | 1 | 0.3 | `24 x 24` | 1.000000 | 2.067033 | 0.971342 |

These checks support the proposed small-instance starting parameters. They do not establish continuous-torus gap bounds, transport convergence, large-system behavior, or performance on an LLM evaluation.

### Suggested English task opening

> Determine whether the geometric charge transported by an interacting, weakly disordered Rice–Mele ring agrees with its many-body Chern number. Implement independent static and time-dependent calculations. Establish the applicable averaging procedure and order of limits, identify any discrepancy, and provide a reproducible numerical error budget. An integer from a single Berry-curvature mesh or a single finite-time trajectory is not sufficient evidence.

## 2. Non-Hermitian skin effect: distinguish physical boundary sensitivity from numerical artifacts

### Physical setup

Consider the disordered Hatano–Nelson chain with open boundaries,

$$
H_g=-t\sum_{j=0}^{N-2}
\left(e^g|j+1\rangle\langle j|+e^{-g}|j\rangle\langle j+1|\right)
+\sum_j v_j|j\rangle\langle j|,
$$

where `t = 1` and the onsite potential is real. Close the chain with a separately controlled boundary amplitude and twist to investigate boundary sensitivity and point-gap topology.

For the open chain,

$$
H_g=S H_0 S^{-1},\qquad S_{jj}=e^{gj},
\qquad \kappa_2(S)=e^{|g|(N-1)}.
$$

Consequently its exact eigenvalues are real even when the similarity transformation is extremely ill-conditioned.

### Coherent sequence of requests

1. Derive the similarity transformation and validate the spectrum against the symmetric problem. Quantify how the numerical discrepancy changes with `N` and `g`.
2. Calculate left and right eigenvectors, suitable localization diagnostics, and eigenvalue sensitivity. Distinguish ordinary Anderson localization from boundary accumulation.
3. Gradually turn on a physical boundary link and compare its spectral effect with the apparent complex eigenvalues produced by finite precision in the nominally open chain.
4. Compute selected pseudospectral contours through the smallest singular value of `z*I - H`, using a stated operator norm and perturbation scale. Explain why a small eigenpair residual need not imply a small eigenvalue error.
5. For the closed, twisted chain, evaluate the point-gap winding using stable determinant phases or factorization. Check that the reference energy remains outside the spectrum throughout the twist cycle.
6. Compare disorder realizations and finite lengths. As an optional extension, relate spectral complexification to the localization length estimated from stabilized transfer matrices.

### Why this can challenge a frontier model

An eigensolver can return smooth, plausible complex spectral curves that contradict an exact structural fact. Merely reducing its convergence tolerance does not resolve a badly conditioned eigenvalue problem. Similarly, treating a pseudospectrum as the actual spectrum or applying a translation-invariant Bloch formula to a disordered chain can yield persuasive but incorrect physics.

### Verification and feasible scope

Use open-chain similarity as the strongest oracle, and test `g = 0`, clean-chain limits, twist refinement, and independent singular-value calculations. Start with matrices of size `40–200`; use selected pseudospectral slices instead of an unnecessarily dense two-dimensional image. A restricted boundary-link perturbation must not be equated with the full set of arbitrary perturbations defining a pseudospectrum.

### Exploratory calculation performed for this proposal

With `v_j = 0.5 * rng.uniform(-1, 1, N)`, a fresh `default_rng(42)` for each length, and ordinary double-precision `numpy.linalg.eigvals`, the open chain gave:

| `N` | `g` | Similarity condition number | Largest computed absolute imaginary part |
| --- | --- | --- | --- |
| 40 | 0.2 | `2.44e3` | `0` |
| 80 | 0.6 | `3.85e20` | approximately `0.23974` |

The imaginary parts in the second case are numerical artifacts: the exact similarity relation requires a real spectrum. The precise artifact size is solver- and platform-dependent and must not become a hardcoded grading target. Recognizing and correcting the contradiction is the target.

### Suggested English task opening

> A standard eigensolver reports complex eigenvalues for an open Hatano–Nelson chain with real onsite disorder. Determine whether these values represent a physical transition or numerical spectral instability. Reconcile exact transformations, eigenpair residuals, pseudospectra, and controlled boundary perturbations, then identify which conclusions survive disorder and finite-size checks.

## 3. Efimov physics: renormalize a three-body integral equation

### Physical setup

Study three identical bosons with short-range interactions in leading-order zero-range effective field theory. Regulate a bound-state Skorniakov–Ter-Martirosian equation with momentum cutoff `Lambda` and include a running three-body coupling `H(Lambda)`.

The final task must supply the exact integral equation, integration measure, regulator definitions, and coupling normalization. Otherwise two correct implementations could disagree only because they use different conventions.

### Coherent sequence of requests

1. Discretize the bound-state integral equation on an appropriate momentum mesh. Locate trimer energies through the kernel's eigenvalue condition and verify quadrature convergence.
2. Demonstrate cutoff dependence without the three-body counterterm. Fix one trimer binding momentum as the renormalization datum, tune `H(Lambda)`, and predict other shallow trimers.
3. Track physical states across cutoff changes without relabeling them when new deep states enter the regulated spectrum.
4. Recover the log-periodic running of the coupling and the Efimov scaling factor. Distinguish poles of a coupling parametrization from a breakdown of the physical calculation.
5. Repeat with a second regulator, refitting the same datum. Show agreement of low-energy observables as cutoff artifacts shrink, even though the running bare couplings need not coincide.

At unitarity, the positive Efimov exponent satisfies

$$
s_0\cosh(\pi s_0/2)=\frac{8}{\sqrt{3}}\sinh(\pi s_0/6),
\qquad s_0\approx1.00624.
$$

The momentum scaling factor is `exp(pi/s_0)`, approximately `22.7`; the adjacent binding-energy ratio is its square, approximately `515`. Confusing these two ratios is an especially clear physics error.

### Why this can challenge a frontier model

Knowing the universal ratio does not produce a renormalized calculation. A solver may inadvertently refit every predicted energy, follow the wrong eigenvalue branch, mistake regulator states for universal trimers, or claim convergence while keeping a momentum mesh that does not resolve the enlarged cutoff interval.

### Verification and feasible scope

Separate momentum-grid convergence from cutoff extrapolation. Compare at least two regulators and more than one predicted observable using only the allowed renormalization datum. Restrict claims to states with binding momenta well below the cutoff and within the assumed zero-range regime.

The workload can use integral kernels with a few hundred quadrature nodes, subject to calibration. Do not require above-breakup scattering or contour deformations in the base version; those would substantially enlarge the numerical problem. The final reference solution must demonstrate that its cutoff span contains enough universal states to make the requested scaling test meaningful.

### Suggested English task opening

> Construct a renormalized numerical solution of the three-boson bound-state equation. Fix one trimer energy, predict the remaining shallow spectrum, and establish discrete scale invariance and regulator independence. Explain why agreement with a memorized Efimov ratio is not sufficient evidence of renormalization.

## 4. Near-extremal Kerr quasinormal modes: follow the physical branch

### Physical setup

Compute gravitational Kerr quasinormal modes with spin weight `s = -2` and time dependence `exp(-i*omega*t)`, in units `G = c = M = 1`. Follow selected `l = m = 2` modes from moderate spin toward `a = 1`, coupling the spin-weighted spheroidal angular problem to the radial outgoing/ingoing boundary-value problem.

Use a bounded near-extremal sequence such as `a = 0, 0.5, 0.9, 0.99, 0.999`. More extreme values should be optional until the reference solver and precision requirements have been established.

### Coherent sequence of requests

1. Implement an angular spectral solver and a radial continued-fraction or independently justified spectral solver. State the angular separation-constant convention explicitly.
2. Track a specified branch by continuation rather than independently sorting roots at each spin. Record the strategy for detecting duplicated roots or a branch switch.
3. Establish separate convergence in angular basis size, radial truncation, root tolerance, and arithmetic precision. A root-finder success flag is not a physical error estimate.
4. Recover the Schwarzschild limit and cross-check moderate-spin values against a pinned reference implementation or supplied reference table.
5. Test the expected approach of a selected zero-damped branch to the horizon frequency and examine damping relative to the surface gravity,

$$
r_+=1+\sqrt{1-a^2},\qquad
\Omega_H=\frac{a}{2r_+},\qquad
\kappa=\frac{\sqrt{1-a^2}}{2r_+}.
$$

For that branch, investigate `Re(omega) -> m*Omega_H` and damping of order `kappa`, with a stated fit window and finite-extremality corrections. Do not assume that every mode family shares this behavior.

### Why this can challenge a frontier model

A solver can converge to the wrong overtone, use an inconsistent separation constant, or interpret interpolation of cached values as an independent numerical solution. Near extremality, apparently stable roots can retain substantial truncation error, and mode labels cannot be assigned reliably by sorting only the imaginary parts.

### Verification and feasible scope

Use one or two carefully specified branches, rather than a full mode atlas. Require independent refinement and checks of the defining equations. The `qnm` package is useful as a moderate-spin comparator, but calling it as the sole solver would bypass the intended task. If the implementation reproduces its algorithm, describe the comparison as a regression check rather than an independent-method validation.

Supply any required reference tables locally. Do not let an optional online cache download become a hidden runtime dependency. No full numerical-relativity simulation is needed.

### Suggested English task opening

> Track a specified gravitational Kerr quasinormal-mode branch toward extremality and determine its damping scaling. Implement the coupled angular and radial calculation, quantify independent sources of numerical error, and demonstrate that apparent convergence is not caused by cached interpolation or a change of mode identity.

## 5. BKT transition: infer a critical temperature without false precision

### Physical setup

Simulate the classical two-dimensional XY model on periodic square lattices,

$$
H=-J\sum_{\langle ij\rangle}\cos(\theta_i-\theta_j),
\qquad J=k_B=1.
$$

Estimate the Berezinskii–Kosterlitz–Thouless transition from the helicity modulus and logarithmic finite-size behavior.

### Coherent sequence of requests

1. Implement a cluster method, with any complementary updates justified, and compare a small system against a separate sampler or high-accuracy reference observable.
2. Derive the helicity modulus from a boundary twist, retaining the connected current-fluctuation term. Verify its normalization and low-temperature behavior.
3. Measure energy, stiffness, and suitable correlation observables over several sizes. Estimate integrated autocorrelation times and effective sample sizes for the actual measured observables.
4. Estimate the transition using BKT finite-size formulas. Compare against a naive universal-jump crossing and explain its finite-size bias.
5. Propagate uncertainty using blocked resampling or independent replicas. Repeat fits after excluding small lattices and varying the correction model. Separate statistical uncertainty from extrapolation uncertainty.
6. Check whether the supported low-temperature correlation behavior is consistent with the transition inferred from the stiffness, rather than treating every observable as an unrelated fit.

### Why this can challenge a frontier model

Thousands of correlated configurations can produce extremely small but invalid error bars. A model may apply an ordinary power-law critical fit to a BKT transition, treat all cluster updates as identical units of Monte Carlo time, or quote a literature transition temperature without numerically establishing it.

### Verification and feasible scope

Begin with sizes such as `L = 16, 32, 64`; include `128` only if the budget permits adequate sampling. Require independent-seed diagnostics and a size-removal stability analysis. Avoid demanding a very precise transition temperature from three small lattices.

The task must distinguish the renormalized stiffness in the universal-jump relation from the finite-torus response used by the estimator. Winding-sector and aspect-ratio conventions matter if the requested precision reaches their correction scale. Either include the appropriate correction or choose a tolerance that does not implicitly ignore it.

Grade statistical consistency and the honesty of the error budget, rather than closeness to a known decimal alone. A single confidence interval missing the reference value is not, by itself, proof of incorrect uncertainty estimation; coverage is an ensemble property.

### Suggested English task opening

> Estimate the BKT transition of the periodic two-dimensional XY model with an uncertainty supported by your simulation. Derive the measured stiffness, diagnose autocorrelation, and distinguish logarithmic finite-size bias from statistical error. Explain whether your available data justify the number of digits you report.

## 6. Analytic continuation: establish what the data cannot determine

### Physical setup

Recover information about a nonnegative fermionic spectral density from noisy imaginary-time data,

$$
G(\tau)=-\int_{-\Omega}^{\Omega}
\frac{e^{-\tau\omega}}{1+e^{-\beta\omega}}A(\omega)\,d\omega,
\qquad A(\omega)\geq0,\qquad \int A(\omega)d\omega=1.
$$

Provide the data, covariance matrix, support assumption, and any exact moment constraints. A possible initial synthetic instance has `beta = 20`, `Omega = 6`, about `80` imaginary-time points, and correlated noise with standard deviation near `1e-4`. These choices need calibration to produce informative but genuinely ambiguous cases.

### Coherent sequence of requests

1. Implement a numerically stable forward kernel, including quadrature weights, and whiten the residual with the supplied covariance.
2. Compare two explicitly regularized reconstruction methods, such as maximum entropy and nonnegative Tikhonov regularization. State how regularization parameters are selected.
3. Examine the singular spectrum of the whitened forward operator and relate it to resolution limits. Test sensitivity to frequency-grid refinement and the assumed support.
4. Construct two materially different nonnegative spectra that both satisfy the same data-consistency criterion and supplied constraints. Demonstrate the ambiguity through forward predictions, not just visually different curves.
5. Bound a specified integrated spectral weight, for example the weight in `[-0.5, 0.5]`, over the allowed confidence region. Distinguish a bound on that observable from a claim to resolve every peak.
6. Repeat on additional noise realizations and examine the calibration of the reported conclusions. Identify which apparent spectral features depend mainly on the regularization prior.

### Why this can challenge a frontier model

A reconstruction can fit the data extremely well and still invent a sharp resonance. Treating a regularized optimum as the unique physical answer, ignoring off-diagonal covariance, or judging success only against one concealed spectrum all misrepresent the inverse problem.

### Verification and feasible scope

The synthetic truth is useful for controlled experiments, but should not be the only grading oracle: multiple spectra may be observationally equivalent. Grade positivity, normalization, forward residuals, sensitivity, and justified bounds.

Define the confidence region precisely. For known Gaussian covariance, one valid construction uses a specified chi-square quantile for the full data vector. Do not casually replace that threshold by an unjustified count of fitted degrees of freedom. Infeasible constraints should be reported and diagnosed rather than repaired by silently relaxing the data.

This task is suitable for dense arrays with a few hundred frequency bins and a standard constrained optimizer. The base version does not require a complete dynamical mean-field theory or quantum Monte Carlo pipeline: those would obscure the central identifiability question.

### Suggested English task opening

> Determine which spectral statements are identifiable from the supplied imaginary-time Green's function and correlated uncertainties. Produce regularized reconstructions, explicitly exhibit a nonunique alternative consistent with the data, and compute bounds on a physically meaningful integrated observable. A visually plausible spectrum is not sufficient evidence.

## How to turn one proposal into a fair, difficult benchmark

1. **Freeze the mathematical problem.** Supply units, conventions, boundary conditions, input arrays, required cases, and acceptance criteria. Keep optional extensions distinct from the mandatory task.
2. **Require executable evidence.** Ask for a terminal entry point, machine-readable results, plots, and a report linking every major conclusion to generated output. Record seeds, package versions, hardware, and runtime. A derivation alone must not satisfy a computational request.
3. **Use independent oracles.** Include at least one exact identity or solvable limit and a second calculation that tests a different failure mode. Two wrappers around the same routine are not independent methods.
4. **Probe generalization.** Keep some small parameter instances or disorder vectors for evaluation, while publishing their distributions and valid domain. Compare against a reference calculation on those instances; do not rely on a single hardcoded number.
5. **Make failure informative.** Require a model to distinguish a physical effect from discretization, conditioning, sampling, model assumptions, or non-identifiability. Correctly establishing that a requested precision is unsupported should receive credit when that is the scientifically warranted conclusion.
6. **Calibrate resources and tolerances.** Run the complete reference workflow on the evaluation machine. Set required sizes and tolerances from that evidence, leaving room for a reasonably efficient alternative implementation.
7. **Pilot with actual models.** Inspect failure modes before claiming the task is hard for a state-of-the-art LLM. If all failures arise from installation or timeouts, improve the environment rather than treating that as evidence of deep physics difficulty.

A possible scoring rubric for Idea 1 is:

| Category | Weight | Evidence |
| --- | --- | --- |
| Correct physical implementation | 20% | Particle-number sector, boundary signs, Hermiticity, free-limit comparison |
| Static topology | 20% | Gauge invariance, nonvanishing overlaps, mesh checks, gap diagnostics |
| Dynamical transport | 25% | Current derivation, persistent-current treatment, norm and integration checks |
| Controlled limits | 20% | Duration, time-step, twist-quadrature, and size analysis |
| Reproducibility and interpretation | 15% | Rerunnable commands, numerical artifacts, and conclusions supported by evidence |

For this task, an unsupported assertion of quantized transport should not receive full physics credit merely because its Chern number happens to be correct. Conversely, conventions producing an overall sign difference should be accepted when derived consistently and verified by loop reversal.

## Research basis and source-access limits

The sources below informed feasibility and verification choices. Direct requests to arXiv were blocked by the environment's network policy. No blocked route was bypassed, and full-text reading of those papers is not claimed. Public software documentation and selected source files were inspected through permitted GitHub access. The two exploratory numerical experiments above were also run locally.

### Public computational sources actually inspected

- [QuSpin README](https://github.com/QuSpin/QuSpin/blob/master/README.md): confirms support for number-conserving many-body calculations, sparse Hamiltonians, and time-dependent quantum dynamics. This supports the computational route for Idea 1; it does not independently validate the proposed pump parameters.
- [qnm README](https://github.com/duetosymmetry/qnm/blob/master/README.md) and [radial solver source](https://github.com/duetosymmetry/qnm/blob/master/qnm/radial.py): describe Kerr angular/radial calculations, continued-fraction truncations, convergence controls, caches, and reference-data comparisons. These directly inform Idea 4's separation of root convergence, truncation, and cached interpolation.
- [TRIQS/maxent README](https://github.com/TRIQS/maxent/blob/unstable/README.md) and [kernel implementation](https://github.com/TRIQS/maxent/blob/unstable/python/triqs_maxent/kernels.py): provide a research implementation of analytic continuation and sign-dependent stable exponential evaluation. The implemented kernel, rather than a documentation formula alone, was checked when specifying Idea 6.
- [ALPS README](https://github.com/ALPSim/ALPS/blob/master/README.md): documents an established computational framework for spin and many-body simulation. It is a software reference for the Monte Carlo workflow, not evidence that the proposed BKT finite-size fit has already been validated.
- [Public tridiagonal pseudospectrum project](https://github.com/jinshanmu/Hatano-Nelson/blob/main/README.md): an exploratory source concerning nonnormal tridiagonal matrices and pseudospectral calculations. Only its README was inspected; its manuscript and formal-proof claims were not audited. Idea 2's exact similarity oracle and reported numerical experiment do not depend on accepting those claims.

These are living repositories. A final benchmark should pin any software it actually depends on to a tested release or commit.

### Foundational literature for the final task design

The following are bibliographic pointers for further verification; their full texts were not retrieved during this session.

- **Charge pumping:** D. J. Thouless, “Quantization of particle transport,” *Physical Review B* **27**, 6083 (1983), [DOI](https://doi.org/10.1103/PhysRevB.27.6083).
- **Discrete Chern calculation:** T. Fukui, Y. Hatsugai, and H. Suzuki, “Chern Numbers in Discretized Brillouin Zone: Efficient Method of Computing (Spin) Hall Conductances,” *Journal of the Physical Society of Japan* **74**, 1674–1677 (2005), [DOI](https://doi.org/10.1143/JPSJ.74.1674).
- **Non-Hermitian localization:** N. Hatano and D. R. Nelson, “Localization Transitions in Non-Hermitian Quantum Mechanics,” *Physical Review Letters* **77**, 570–573 (1996), [DOI](https://doi.org/10.1103/PhysRevLett.77.570).
- **Three-body renormalization:** P. F. Bedaque, H.-W. Hammer, and U. van Kolck, “Renormalization of the Three-Body System with Short-Range Interactions,” *Physical Review Letters* **82**, 463–467 (1999), [DOI](https://doi.org/10.1103/PhysRevLett.82.463).
- **Efimov universality:** E. Braaten and H.-W. Hammer, “Universality in Few-body Systems with Large Scattering Length,” *Physics Reports* **428**, 259–390 (2006), [arXiv:cond-mat/0410417](https://arxiv.org/abs/cond-mat/0410417).
- **Kerr continued fractions:** E. W. Leaver, “An Analytic Representation for the Quasi-Normal Modes of Kerr Black Holes,” *Proceedings of the Royal Society A* **402**, 285–298 (1985), [DOI](https://doi.org/10.1098/rspa.1985.0119).
- **Kerr reference software:** L. C. Stein, “qnm: A Python package for calculating Kerr quasinormal modes, separation constants, and spherical-spheroidal mixing coefficients,” *Journal of Open Source Software* **4**, 1683 (2019), [DOI](https://doi.org/10.21105/joss.01683). The software README was inspected; the article full text was not.
- **BKT finite-size numerics:** H. Weber and P. Minnhagen, “Monte Carlo determination of the critical temperature for the two-dimensional XY model,” *Physical Review B* **37**, 5986–5989 (1988), [DOI](https://doi.org/10.1103/PhysRevB.37.5986).
- **Analytic continuation:** M. Jarrell and J. E. Gubernatis, “Bayesian inference and the analytic continuation of imaginary-time quantum Monte Carlo data,” *Physics Reports* **269**, 133–195 (1996), [DOI](https://doi.org/10.1016/0370-1573(95)00074-7).
