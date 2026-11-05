"""Author-side oracle. Exclude this directory from the candidate workspace."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
import time
from pathlib import Path

import numpy as np
import scipy
from scipy.integrate import solve_ivp
from scipy.linalg import eigh
from numpy.polynomial.legendre import leggauss

ROOT = Path(__file__).resolve().parents[1]
INPUTS = ROOT / "benchmark" / "inputs.json"


def basis(L):
    return [x for x in range(1 << L) if x.bit_count() == L // 2]


def hop_matrix(L, source, target):
    """Apply annihilation then creation, counting lower occupied orbitals."""
    states = basis(L)
    index = {x: i for i, x in enumerate(states)}
    result = np.zeros((len(states), len(states)), dtype=complex)
    for column, state in enumerate(states):
        if not (state >> source) & 1 or (state >> target) & 1:
            continue
        sign = (-1) ** ((state & ((1 << source) - 1)).bit_count())
        intermediate = state ^ (1 << source)
        sign *= (-1) ** ((intermediate & ((1 << target) - 1)).bit_count())
        result[index[intermediate | (1 << target)], column] = sign
    return result


def build_operators(case, theta, phi):
    L = case["L"]
    states = basis(L)
    occupation = np.array([[(x >> j) & 1 for j in range(L)] for x in states])
    centered = occupation - 0.5
    loop = case["loop"]
    delta = loop["delta_offset"] + loop["delta_radius"] * np.cos(theta)
    ddelta = -loop["delta_radius"] * np.sin(theta)
    stagger = loop["stagger_radius"] * np.sin(theta)
    dstagger = loop["stagger_radius"] * np.cos(theta)
    stagger_operator = centered @ ((-1.0) ** np.arange(L))
    diagonal = (
        case["V"] * np.sum(centered * np.roll(centered, -1, axis=1), axis=1)
        + case["W"] * (centered @ np.asarray(case["disorder"]))
        + stagger * stagger_operator
    )
    H = np.diag(diagonal).astype(complex)
    dtheta = np.diag(dstagger * stagger_operator).astype(complex)
    dphi = np.zeros_like(H)
    for j in range(L):
        A = hop_matrix(L, j, (j + 1) % L)
        phase = np.exp(1j * phi) if j == L - 1 else 1.0
        coefficient = -(1 + (-1) ** j * delta) * phase
        derivative = -((-1) ** j) * ddelta * phase
        H += coefficient * A + coefficient.conjugate() * A.conj().T
        dtheta += derivative * A + derivative.conjugate() * A.conj().T
        if j == L - 1:
            dc = 1j * coefficient
            dphi += dc * A + dc.conjugate() * A.conj().T
    return {"H": H, "dtheta": dtheta, "dphi": dphi}


def spectral_quantities(operators):
    energies, vectors = eigh(operators["H"])
    u = vectors[:, 0]
    a = vectors.conj().T @ operators["dtheta"] @ u
    b = vectors.conj().T @ operators["dphi"] @ u
    curvature = 2 * np.imag(np.sum(a[1:].conj() * b[1:] / (energies[1:] - energies[0]) ** 2))
    return energies, u, float(curvature)


def kubo_curvature(case, theta, phi):
    return spectral_quantities(build_operators(case, theta, phi))[2]


def adiabatic_charge(case, phi, nodes=128, curvature_function=kubo_curvature):
    x, w = leggauss(nodes)
    return float(np.pi * sum(weight * curvature_function(case, np.pi * (node + 1), phi) for node, weight in zip(x, w)))


def static_observables(case, mesh, rephase=False):
    states = np.empty((mesh, mesh, math.comb(case["L"], case["L"] // 2)), complex)
    gap = np.inf
    curvature_sum = 0.0
    for a in range(mesh):
        for b in range(mesh):
            theta, phi = 2 * np.pi * a / mesh, 2 * np.pi * b / mesh
            energies, u, curvature = spectral_quantities(build_operators(case, theta, phi))
            if rephase:
                u = u * np.exp(2j * np.pi * ((37 * a + 17 * b + 13 * a * b) % 101) / 101)
            states[a, b] = u
            gap = min(gap, energies[1] - energies[0])
            curvature_sum += curvature
    U = np.sum(states.conj() * np.roll(states, -1, axis=0), axis=2)
    V = np.sum(states.conj() * np.roll(states, -1, axis=1), axis=2)
    min_overlap = min(np.min(np.abs(U)), np.min(np.abs(V)))
    U /= np.abs(U)
    V /= np.abs(V)
    plaquettes = np.angle(U * np.roll(V, -1, axis=0) * np.roll(U, -1, axis=1).conj() * V.conj())
    return {
        "chern": float(plaquettes.sum() / (2 * np.pi)),
        "sampled_min_gap": float(gap),
        "min_link_overlap": float(min_overlap),
        "max_abs_plaquette_phase": float(np.max(np.abs(plaquettes))),
        "kubo_integral": float(2 * np.pi * curvature_sum / mesh ** 2),
    }


def drive(s):
    return 10 * s ** 3 - 15 * s ** 4 + 6 * s ** 5


def evolve(case, phi, T, direction, rtol=1e-10, atol=1e-12, eq_nodes=128):
    """DOP853 plus a separate ground-state quadrature for persistent current."""
    samples = [build_operators(case, theta, phi) for theta in [0, np.pi / 2, np.pi]]
    components = {}
    for key in ["H", "dphi"]:
        constant = (samples[0][key] + samples[2][key]) / 2
        cosine = (samples[0][key] - samples[2][key]) / 2
        sine = samples[1][key] - constant
        components[key] = (constant, cosine, sine)

    def at(s, key):
        theta = direction * 2 * np.pi * drive(s)
        c, a, b = components[key]
        return c + np.cos(theta) * a + np.sin(theta) * b

    _, vectors = eigh(samples[0]["H"], subset_by_index=[0, 0])
    initial = vectors[:, 0]
    dimension = len(initial)

    def rhs(t, y):
        state = y[:dimension]
        H = at(t / T, "H")
        current = -at(t / T, "dphi")
        return np.r_[-1j * (H @ state), float(np.vdot(state, current @ state).real)]

    solution = solve_ivp(rhs, [0.0, T], np.r_[initial, 0j], method="DOP853", rtol=rtol, atol=atol)
    if not solution.success:
        raise RuntimeError(solution.message)
    final = solution.y[:dimension, -1]
    q_total = float(solution.y[dimension, -1].real)
    nodes, weights = leggauss(eq_nodes)
    equilibrium_integral = 0.0
    for s, weight in zip((nodes + 1) / 2, weights / 2):
        _, u = eigh(at(s, "H"), subset_by_index=[0, 0])
        equilibrium_integral += weight * np.vdot(u[:, 0], -at(s, "dphi") @ u[:, 0]).real
    q_equilibrium = float(T * equilibrium_integral)
    norm_squared = float(np.vdot(final, final).real)
    probability = float(1 - abs(np.vdot(initial, final)) ** 2 / norm_squared)
    return {
        "q_total": q_total,
        "q_equilibrium": q_equilibrium,
        "q_geometric": q_total - q_equilibrium,
        "final_norm_squared": norm_squared,
        "excitation_probability": probability,
        "final_state": final,
    }


def trajectory_specs(data):
    """Ordered unique probes shared by reference generation and grading."""
    d = data["dynamics"]
    fixed = d["fixed_phi"]
    specs = []

    def add(case, phi, T, direction):
        item = (case, float(phi), float(T), int(direction))
        if item not in specs:
            specs.append(item)

    for T in d["periods"]:
        add("pump_L6", fixed, T, 1)
    add("pump_L6", fixed, d["reverse_period"], -1)
    add("free_L6", fixed, 30.0, 1)
    add("pump_L8", fixed, d["L8_period"], 1)
    add("trivial_L6", fixed, d["control_period"], 1)
    for M in d["twist_counts"]:
        for T in d["average_periods"]:
            for b in range(M):
                add("pump_L6", 2 * np.pi * (b + 0.5) / M, T, 1)
    for b in range(8):
        phi = 2 * np.pi * (b + 0.5) / 8
        add("pump_L6", phi, d["reverse_period"], -1)
        add("trivial_L6", phi, d["control_period"], 1)
    return specs


def serializable_trajectory(result):
    data = dict(result)
    data["final_state"] = [[float(z.real), float(z.imag)] for z in data["final_state"]]
    return data


def generate(data, tight=False):
    start = time.perf_counter()
    result = {"static": [], "curvature": [], "adiabatic_lines": [], "dynamics": []}
    for name, case in data["cases"].items():
        for mesh in data["static_meshes"]:
            values = static_observables(case, mesh)
            result["static"].append({"case": name, "mesh": mesh, "values": values})
            print("static", name, mesh, values, flush=True)
    for probe in data["operator_probes"]:
        case = data["cases"][probe["case"]]
        value = kubo_curvature(case, probe["theta"], probe["phi"])
        result["curvature"].append({**probe, "value": value})
    for name, case in data["cases"].items():
        phi = data["dynamics"]["fixed_phi"]
        result["adiabatic_lines"].append({"case": name, "phi": phi, "value": adiabatic_charge(case, phi)})
    for number, (name, phi, T, direction) in enumerate(trajectory_specs(data)):
        options = {"rtol": 2e-12, "atol": 2e-14, "eq_nodes": 256} if tight else {}
        values = evolve(data["cases"][name], phi, T, direction, **options)
        result["dynamics"].append({"case": name, "phi": phi, "T": T, "direction": direction, "values": serializable_trajectory(values)})
        if number % 10 == 0:
            print("trajectory", number, name, "T", T, "direction", direction, "q", values["q_geometric"], flush=True)
    result["metadata"] = {
        "input_sha256": hashlib.sha256(INPUTS.read_bytes()).hexdigest(),
        "python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__,
        "tight": tight, "elapsed_seconds": time.perf_counter() - start,
    }
    return result


def summarize(result):
    lines = []
    data = json.loads(INPUTS.read_text())
    for M in data["dynamics"]["twist_counts"]:
        phis = [float(2 * np.pi * (b + 0.5) / M) for b in range(M)]
        for T in data["dynamics"]["average_periods"]:
            group = [x for x in result["dynamics"] if x["case"] == "pump_L6" and x["T"] == T and x["direction"] == 1 and x["phi"] in phis]
            lines.append({"case": "pump_L6", "twists": M, "T": T, "q_average": float(np.mean([x["values"]["q_geometric"] for x in group]))})
    return lines


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--tight", action="store_true")
    args = parser.parse_args()
    data = json.loads(INPUTS.read_text())
    results = generate(data, tight=args.tight)
    results["averages"] = summarize(results)
    args.output.write_text(json.dumps(results, indent=2, allow_nan=False) + "\n")
    print(json.dumps({"metadata": results["metadata"], "averages": results["averages"]}, indent=2))


if __name__ == "__main__":
    main()
