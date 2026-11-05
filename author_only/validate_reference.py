"""Independent identities, convergence evidence, and grader fault probes."""

import argparse
import copy
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from scipy.linalg import eigh

import grade
import reference as ref


def jordan_wigner_hamiltonian(case, theta, phi):
    """Independent full-Fock-space construction using tensor products."""
    L = case["L"]
    identity = np.eye(2)
    parity = np.diag([1.0, -1.0])
    lower = np.array([[0.0, 1.0], [0.0, 0.0]])
    annihilators = []
    for j in range(L):
        a = np.ones((1, 1))
        for k in reversed(range(L)):
            a = np.kron(a, identity if k > j else lower if k == j else parity)
        annihilators.append(a)
    full_identity = np.eye(1 << L)
    densities = [a.T @ a - 0.5 * full_identity for a in annihilators]
    H = np.zeros_like(full_identity, dtype=complex)
    loop = case["loop"]
    delta = loop["delta_offset"] + loop["delta_radius"] * np.cos(theta)
    stagger = loop["stagger_radius"] * np.sin(theta)
    for j in range(L):
        k = (j + 1) % L
        forward = annihilators[k].T @ annihilators[j]
        z = -(1 + (-1) ** j * delta) * np.exp(1j * phi * (j == L - 1))
        H += z * forward + z.conjugate() * forward.T
        H += ((-1) ** j * stagger + case["W"] * case["disorder"][j]) * densities[j]
        H += case["V"] * (densities[j] @ densities[k])
    selected = ref.basis(L)
    return H[np.ix_(selected, selected)]


def projector(case, theta, phi):
    _, u = eigh(ref.build_operators(case, theta, phi)["H"], subset_by_index=[0, 0])
    return u @ u.conj().T


def expect_rejection(operation):
    try:
        operation()
    except (AssertionError, ValueError, TypeError):
        return
    raise AssertionError("A deliberately incorrect result was accepted")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    inputs = json.loads(ref.INPUTS.read_text())
    tight = json.loads(Path(__file__).with_name("reference_results.json").read_text())
    baseline = json.loads(args.baseline.read_text())
    assert baseline["metadata"]["input_sha256"] == tight["metadata"]["input_sha256"]
    evidence = {"checks": [], "maxima": {}}

    def record(label, operation):
        operation()
        evidence["checks"].append(label)

    maxima = {"jordan_wigner_error": 0.0, "derivative_error": 0.0, "projector_curvature_error": 0.0}
    for probe in inputs["operator_probes"]:
        case = inputs["cases"][probe["case"]]
        theta, phi = probe["theta"], probe["phi"]
        ops = ref.build_operators(case, theta, phi)
        H = jordan_wigner_hamiltonian(case, theta, phi)
        error = float(np.max(np.abs(H - ops["H"])))
        maxima["jordan_wigner_error"] = max(maxima["jordan_wigner_error"], error)
        assert error < 1e-12
        for variable, derivative in [("theta", "dtheta"), ("phi", "dphi")]:
            step = 2e-5
            angles = np.array([theta, phi])
            axis = 0 if variable == "theta" else 1
            angles[axis] += step
            plus = ref.build_operators(case, *angles)["H"]
            angles[axis] -= 2 * step
            minus = ref.build_operators(case, *angles)["H"]
            error = float(np.max(np.abs((plus - minus) / (2 * step) - ops[derivative])))
            maxima["derivative_error"] = max(maxima["derivative_error"], error)
            assert error < 2e-8
        step = 2e-5
        P = projector(case, theta, phi)
        Pt = (projector(case, theta + step, phi) - projector(case, theta - step, phi)) / (2 * step)
        Pp = (projector(case, theta, phi + step) - projector(case, theta, phi - step)) / (2 * step)
        curvature = float((-1j * np.trace(P @ (Pt @ Pp - Pp @ Pt))).real)
        error = abs(curvature - ref.kubo_curvature(case, theta, phi))
        maxima["projector_curvature_error"] = max(maxima["projector_curvature_error"], error)
        assert error < 2e-8
    evidence["checks"] += ["independent full-Fock Jordan-Wigner matrices", "finite-difference operator derivatives", "independent projector curvature and sign"]

    free = inputs["cases"]["free_L6"]
    theta, phi = 0.73, 1.31
    L = free["L"]
    one = np.diag(np.sin(theta) * (-1.0) ** np.arange(L)).astype(complex)
    for j in range(L):
        k = (j + 1) % L
        z = -(1 + (-1) ** j * 0.5 * np.cos(theta)) * np.exp(1j * phi * (j == L - 1))
        one[k, j] += z
        one[j, k] += z.conjugate()
    energies = eigh(one, eigvals_only=True)
    many = eigh(ref.build_operators(free, theta, phi)["H"], eigvals_only=True)[0]
    record("free ground-state energy from occupied single-particle levels", lambda: np.testing.assert_allclose(many, np.sum(energies[:L // 2]), atol=1e-12))

    case = inputs["cases"]["pump_L6"]
    original = ref.static_observables(case, 12)
    rephased = ref.static_observables(case, 12, rephase=True)
    for key in original:
        np.testing.assert_allclose(original[key], rephased[key], atol=2e-12, rtol=0)
    reversed_case = copy.deepcopy(case)
    reversed_case["loop"]["stagger_radius"] *= -1
    close_chern = ref.static_observables(reversed_case, 12)["chern"]
    assert abs(close_chern + 1) < 1e-12
    evidence["checks"] += ["deterministic eigenvector rephasing invariance", "reversed-loop Chern orientation"]

    difference = {key: 0.0 for key in ["q_total", "q_equilibrium", "q_geometric", "excitation_probability", "final_norm_squared", "final_state_infidelity"]}
    assert len(tight["dynamics"]) == len(baseline["dynamics"])
    for a, b in zip(tight["dynamics"], baseline["dynamics"]):
        assert all(a[k] == b[k] for k in ["case", "phi", "T", "direction"])
        for key in difference:
            if key == "final_state_infidelity":
                u = np.array([complex(*x) for x in a["values"]["final_state"]])
                v = np.array([complex(*x) for x in b["values"]["final_state"]])
                delta = max(0.0, float(1 - abs(np.vdot(u, v)) ** 2 / (np.vdot(u, u).real * np.vdot(v, v).real)))
            else:
                delta = abs(a["values"][key] - b["values"][key])
            difference[key] = max(difference[key], delta)
    assert max(difference[k] for k in ["q_total", "q_equilibrium", "q_geometric"]) < 3e-8
    assert difference["excitation_probability"] < 3e-8
    assert difference["final_state_infidelity"] < 1e-10
    evidence["checks"].append("all 71 trajectories: tighter ODE tolerances and doubled equilibrium quadrature")

    # Independent, exactly unitary frozen-midpoint propagation at two resolutions.
    entry = tight["dynamics"][0]
    phi, T = entry["phi"], entry["T"]
    samples = [ref.build_operators(case, angle, phi) for angle in [0, np.pi / 2, np.pi]]
    components = {}
    for key in ["H", "dphi"]:
        c = (samples[0][key] + samples[2][key]) / 2
        components[key] = (c, (samples[0][key] - samples[2][key]) / 2, samples[1][key] - c)
    totals = []
    for steps in [2000, 4000]:
        _, u = eigh(samples[0]["H"], subset_by_index=[0, 0])
        state = u[:, 0]
        dt = T / steps
        charge = 0.0
        for n in range(steps):
            angle = 2 * np.pi * ref.drive((n + 0.5) / steps)
            operators = {key: c + np.cos(angle) * a + np.sin(angle) * b for key, (c, a, b) in components.items()}
            energies, vectors = eigh(operators["H"])
            coefficients = vectors.conj().T @ state
            midpoint = vectors @ (np.exp(-0.5j * energies * dt) * coefficients)
            charge += dt * np.vdot(midpoint, -operators["dphi"] @ midpoint).real
            state = vectors @ (np.exp(-1j * energies * dt) * coefficients)
        assert abs(np.vdot(state, state).real - 1) < 1e-10
        totals.append(charge)
    extrapolated = (4 * totals[1] - totals[0]) / 3
    maxima["independent_dynamics_charge_error"] = abs(extrapolated - entry["values"]["q_total"])
    assert maxima["independent_dynamics_charge_error"] < 3e-7
    evidence["checks"].append("independent unitary midpoint evolution and Richardson-refined current integral")

    for line in tight["adiabatic_lines"]:
        c = inputs["cases"][line["case"]]
        error = abs(ref.adiabatic_charge(c, line["phi"], nodes=256) - line["value"])
        assert error < 1e-10
    evidence["checks"].append("fixed-flux adiabatic curvature quadrature: 128 versus 256 nodes")

    tolerances = inputs["tolerances"]
    probe = inputs["operator_probes"][1]
    for fault in ["drop_interaction", "flip_current"]:
        def bad_operators(c, t, p, fault=fault):
            if fault == "drop_interaction":
                c["V"] = 0
            result = ref.build_operators(c, t, p)
            if fault == "flip_current":
                result["dphi"] *= -1
            return result
        candidate = SimpleNamespace(build_operators=bad_operators)
        expect_rejection(lambda: grade.check_operators(candidate, case, probe["theta"], probe["phi"], tolerances))
        evidence["checks"].append("grader rejects " + fault)
    entry = tight["dynamics"][2]
    for fault in ["omit_persistent_current", "reverse_geometric_sign", "replace_transport_by_integer"]:
        result = copy.deepcopy(entry["values"])
        result["final_state"] = np.array([complex(*z) for z in result["final_state"]])
        if fault == "omit_persistent_current":
            result["q_equilibrium"] = 0.0
            result["q_geometric"] = result["q_total"]
        elif fault == "reverse_geometric_sign":
            result["q_geometric"] *= -1
        else:
            result["q_geometric"] = 1.0
        candidate = SimpleNamespace(evolve=lambda *a, **k: result)
        expect_rejection(lambda: grade.check_evolution(candidate, case, entry, tolerances))
        evidence["checks"].append("grader rejects " + fault)
    evidence["maxima"] = {**maxima, "baseline_vs_tight": difference}
    evidence["reference_runtime_seconds"] = tight["metadata"]["elapsed_seconds"]
    evidence["baseline_runtime_seconds"] = baseline["metadata"]["elapsed_seconds"]
    evidence["environment"] = {k: tight["metadata"][k] for k in ["python", "numpy", "scipy"]}
    evidence["input_sha256"] = tight["metadata"]["input_sha256"]
    args.output.write_text(json.dumps(evidence, indent=2, allow_nan=False) + "\n")
    print(json.dumps(evidence, indent=2))


if __name__ == "__main__":
    main()
