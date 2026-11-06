"""Public format and consistency checks; contains no reference answers."""

import argparse
import copy
import hashlib
import importlib.util
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.linalg import eigh
from jsonschema import Draft202012Validator

HERE = Path(__file__).resolve().parent
CHARGES = ("q_total", "q_equilibrium", "q_geometric")
FRAME_CASES = ("pump_L6", "trivial_L6")
FRAME_TRAJECTORIES = ("trajectory_000", "trajectory_003", "trajectory_005")
TRANSITION_OFFSETS = (0.49, 0.5, 0.51)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def no_duplicates(pairs):
    result = {}
    for key, value in pairs:
        require(key not in result, f"Duplicate JSON key: {key}")
        result[key] = value
    return result


def reject_constant(value):
    raise ValueError(f"Invalid JSON numeric constant: {value}")


def load_json(path):
    return json.loads(Path(path).read_text(), object_pairs_hook=no_duplicates, parse_constant=reject_constant)


def check_finite(value):
    if isinstance(value, float):
        require(math.isfinite(value), "Nonfinite number")
    elif isinstance(value, dict):
        for item in value.values():
            check_finite(item)
    elif isinstance(value, list):
        for item in value:
            check_finite(item)


def trajectory_plan(data):
    d = data["dynamics"]
    specs = []

    def add(case, phi, T, direction):
        item = {"case": case, "phi": float(phi), "T": float(T), "direction": direction}
        if item not in specs:
            specs.append(item)

    for T in d["periods"]:
        add("pump_L6", d["fixed_phi"], T, 1)
    add("pump_L6", d["fixed_phi"], d["reverse_period"], -1)
    add("free_L6", d["fixed_phi"], 30, 1)
    add("pump_L8", d["fixed_phi"], d["L8_period"], 1)
    add("trivial_L6", d["fixed_phi"], d["control_period"], 1)
    for K in d["twist_counts"]:
        for T in d["average_periods"]:
            for b in range(K):
                add("pump_L6", 2 * math.pi * (b + 0.5) / K, T, 1)
    for b in range(8):
        phi = 2 * math.pi * (b + 0.5) / 8
        add("pump_L6", phi, d["reverse_period"], -1)
        add("trivial_L6", phi, d["control_period"], 1)
    return [{"id": f"trajectory_{n:03d}", **spec} for n, spec in enumerate(specs)]


def average_plan(data):
    d = data["dynamics"]
    return [
        {"case": "pump_L6", "twists": K, "T": T, "direction": 1}
        for K in d["twist_counts"] for T in d["average_periods"]
    ] + [
        {"case": "pump_L6", "twists": 8, "T": d["reverse_period"], "direction": -1},
        {"case": "trivial_L6", "twists": 8, "T": d["control_period"], "direction": 1},
    ]


def keyed(rows, fields):
    result = {}
    for row in rows:
        key = tuple(row[field] for field in fields)
        require(key not in result, f"Duplicate record: {key}")
        result[key] = row
    return result


def close(actual, expected, tolerance=1e-12, label="value"):
    require(abs(actual - expected) <= tolerance, f"{label}: got {actual}, expected {expected}, atol {tolerance}")


def check_coordinates(actual, expected, fields):
    for field in fields:
        if isinstance(expected[field], str):
            require(actual[field] == expected[field], f"Incorrect {field}")
        else:
            close(actual[field], expected[field], label=field)


def _complex(value):
    a = np.asarray(value, dtype=float)
    require(a.shape[-1] == 2, "Complex numbers must be [real, imaginary] pairs")
    return a[..., 0] + 1j * a[..., 1]


def _near(actual, expected, tolerance=2e-7, label="array"):
    actual, expected = np.asarray(actual), np.asarray(expected)
    require(actual.shape == expected.shape, f"Wrong shape: {label}")
    require(np.all(np.isfinite(actual)), f"Nonfinite {label}")
    require(np.max(np.abs(actual - expected)) <= tolerance, f"Incorrect {label}")


def _frame_data(L, phi, s=0.0, T=1.0):
    states = [x for x in range(1 << L) if x.bit_count() == L // 2]
    occ = np.array([[(x >> j) & 1 for j in range(L)] for x in states])
    X = occ @ (np.arange(L) / L)
    A = occ @ (((-1.0) ** np.arange(L)) * (np.arange(L) + 1) / (5 * L))
    G = np.exp(1j * (phi * X + A * np.sin(2 * np.pi * s) + phi * s * (L // 2) / 4))
    dot_angle = A * (2 * np.pi / T) * np.cos(2 * np.pi * s) + phi * (L // 2) / (4 * T)
    K = X + s * (L // 2) / 4
    return G, dot_angle, K


def _frame_metrics(states, L):
    M = states.shape[0]
    U = np.sum(states.conj() * np.roll(states, -1, axis=0), axis=2)
    V = np.sum(states.conj() * np.roll(states, -1, axis=1), axis=2)
    unsewn_V = V.copy()
    sewing = _frame_data(L, 2 * np.pi)[0]
    V[:, -1] = np.sum(states[:, -1].conj() * sewing * states[:, 0], axis=1)
    minimum = float(min(np.min(abs(U)), np.min(abs(V))))
    require(minimum > 1e-12 and np.min(abs(unsewn_V)) > 1e-12, "Undefined frame link")
    U /= abs(U)
    V /= abs(V)
    unsewn_V /= abs(unsewn_V)
    def flux(v):
        return np.angle(U * np.roll(v, -1, axis=0) * np.roll(U, -1, axis=1).conj() * v.conj())
    phases = flux(V)
    return {"sewn_chern": float(phases.sum() / (2 * np.pi)),
            "unsewn_chern": float(flux(unsewn_V).sum() / (2 * np.pi)),
            "min_link_overlap": minimum, "max_abs_plaquette_phase": float(np.max(abs(phases)))}


def _check_frame_static(record, case, api=None, encoded=True):
    M, L = record["mesh"], case["L"]
    dim = math.comb(L, L // 2)
    states = _complex(record["ground_frames"]) if encoded else np.asarray(record["ground_frames"], dtype=complex)
    require(states.shape == (M, M, dim), "Wrong frame-grid dimensions")
    require(np.max(abs(np.sum(abs(states) ** 2, axis=2) - 1)) <= 1e-7, "Frame states not normalized")
    values = _frame_metrics(states, L)
    for key in values:
        close(record["values"][key], values[key], 2e-8, "Frame topology: " + key)
    if api is not None:
        for a in range(M):
            for b in range(M):
                theta, phi = 2 * np.pi * a / M, 2 * np.pi * b / M
                H = np.asarray(api.build_operators(copy.deepcopy(case), theta, phi)["H"])
                e0 = eigh(H, eigvals_only=True, subset_by_index=[0, 0])[0]
                G = _frame_data(L, phi)[0]
                bare = G.conj() * states[a, b]
                require(np.max(abs(H @ bare - e0 * bare)) <= 2e-7, "Frame vector is not a ground state")


def _check_frame_dynamic(values, case, original, api=None, encoded=True):
    L = case["L"]
    phi, T = original["phi"], original["T"]
    dim = math.comb(L, L // 2)
    convert = _complex if encoded else lambda value: np.asarray(value, dtype=complex)
    vectors = {key: convert(values[key]) for key in ["canonical_initial_state", "canonical_final_state", "frame_initial_state", "frame_final_state"]}
    for key, vector in vectors.items():
        require(vector.shape == (dim,), "Wrong state dimension: " + key)
        close(float(np.vdot(vector, vector).real), 1, 1e-7, "State norm: " + key)
    ci, cf = vectors["canonical_initial_state"], vectors["canonical_final_state"]
    fi, ff = vectors["frame_initial_state"], vectors["frame_final_state"]
    G0, _, K0 = _frame_data(L, phi, 0, T)
    G1, _, K1 = _frame_data(L, phi, 1, T)
    _near(fi, G0 * ci, label="initial frame covariance")
    _near(ff, G1 * cf, label="final frame covariance (phase-sensitive)")
    old = _complex(original["values"]["final_state"])
    fidelity = abs(np.vdot(old, cf)) ** 2 / (np.vdot(old, old).real * np.vdot(cf, cf).real)
    require(1 - fidelity <= 1e-8, "Canonical endpoint disagrees with baseline trajectory")
    close(values["q_boundary"], original["values"]["q_total"], 3e-6, "Boundary charge")
    change = float(np.vdot(cf, K1 * cf).real - np.vdot(ci, K0 * ci).real)
    close(values["polarization_change"], change, 3e-6, "Frame polarization change")
    close(values["q_parameter"] - values["q_boundary"], change, 3e-6, "Flux-derivative current relation")
    _near(convert(values["ramsey_raw"]), np.vdot(cf, ff), 3e-7, "Raw Ramsey amplitude")
    _near(convert(values["ramsey_covariant"]), np.vdot(cf, G1.conj() * ff), 3e-7, "Covariant Ramsey amplitude")
    if api is not None:
        H = np.asarray(api.build_operators(copy.deepcopy(case), 0.0, phi)["H"])
        e0 = eigh(H, eigvals_only=True, subset_by_index=[0, 0])[0]
        require(np.max(abs(H @ ci - e0 * ci)) <= 2e-7, "Canonical initial state is not a ground state")


def _check_transition(record, data, api=None, encoded=True):
    offset = record["delta_offset"]
    singular = abs(offset - 0.5) < 1e-12
    require(record["status"] == ("singular" if singular else "gapped"), "Incorrect bundle validity classification")
    close(record["minimum_gap"], 4 * abs(offset - 0.5), 2e-8, "Continuous minimum gap")
    for coordinate in ["theta", "phi"]:
        close(math.remainder(record["witness"][coordinate] - math.pi, 2 * math.pi), 0.0, 1e-6, "Gap witness")
    require(record["ground_rank"] == (2 if singular else 1), "Incorrect ground-space rank")
    if singular:
        require(record["chern"] is None and record["mesh_values"] == [], "A singular rank-one bundle has no Chern result")
    else:
        close(record["chern"], 1.0 if offset < 0.5 else 0.0, 2e-8, "Gapped Chern number")
        meshes = keyed(record["mesh_values"], ["mesh"])
        require(set(meshes) == {(M,) for M in data["static_meshes"]}, "Incomplete transition mesh diagnostics")
    P = _complex(record["ground_projector"]) if encoded else np.asarray(record["ground_projector"], dtype=complex)
    dim = math.comb(6, 3)
    require(P.shape == (dim, dim), "Incorrect ground projector dimensions")
    _near(P, P.conj().T, label="Hermitian ground projector")
    _near(P @ P, P, label="Idempotent ground projector")
    close(float(np.trace(P).real), record["ground_rank"], 2e-7, "Projector rank")
    if api is not None:
        case = copy.deepcopy(data["cases"]["free_L6"])
        case["loop"]["delta_offset"] = offset
        H = np.asarray(api.build_operators(case, math.pi, math.pi)["H"])
        e, U = eigh(H)
        rank = record["ground_rank"]
        target = U[:, :rank] @ U[:, :rank].conj().T
        _near(P, target, label="Lowest-energy ground projector")
        close(e[1] - e[0], record["minimum_gap"], 2e-8, "Witness gap")
        if not singular:
            for row in record["mesh_values"]:
                expected = api.static_observables(copy.deepcopy(case), row["mesh"])
                for key in ["chern", "kubo_integral", "sampled_min_gap"]:
                    close(row[key], expected[key], 2e-8, "Unrounded transition mesh: " + key)


def _check_frame_operator(api, case, theta, phi, s, T):
    baseline = api.build_operators(copy.deepcopy(case), theta, phi)
    result = api.frame_operators(copy.deepcopy(case), theta, phi, s, T)
    H, dp, I = (np.asarray(result[key]) for key in ["H", "dphi", "boundary_current"])
    G, dot_angle, _ = _frame_data(case["L"], phi, s, T)
    baseH, baseI = np.asarray(baseline["H"]), -np.asarray(baseline["dphi"])
    _near(H * G[None, :] - G[:, None] * baseH, -np.diag(dot_angle * G), 2e-10, "Schrodinger frame covariance")
    _near(I, G[:, None] * baseI * G.conj()[None, :], 2e-10, "Measured boundary-current operator")
    h = 2e-4
    def shifted(d):
        return np.asarray(api.frame_operators(copy.deepcopy(case), theta, phi + d, s, T)["H"])
    numerical = (shifted(-2*h) - 8*shifted(-h) + 8*shifted(h) - shifted(2*h)) / (12*h)
    _near(dp, numerical, 2e-8, "Flux derivative of the frame Hamiltonian")
    for value in [H, dp, I]:
        _near(value, value.conj().T, 2e-10, "Hermitian frame operator")


def _load_api(path):
    sys.path.insert(0, str(Path(path).resolve().parent))
    spec = importlib.util.spec_from_file_location("frame_submission", Path(path).resolve())
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _validate_extensions(document, data, api=None):
    frames = keyed(document["frame_static"], ["case", "mesh"])
    require(set(frames) == {(name, M) for name in FRAME_CASES for M in data["static_meshes"]}, "Incorrect static-frame coverage")
    baseline_static = keyed(document["static"], ["case", "mesh", "rephase"])
    for (name, M), row in frames.items():
        _check_frame_static(row, data["cases"][name], api)
        close(row["values"]["sewn_chern"], baseline_static[name, M, False]["values"]["chern"], 2e-8, "Frame-independent Chern number")
    dynamics = keyed(document["frame_dynamics"], ["trajectory_id"])
    require(set(dynamics) == {(name,) for name in FRAME_TRAJECTORIES}, "Incorrect dynamical-frame coverage")
    baseline_dynamics = keyed(document["dynamics"], ["id"])
    for (name,), row in dynamics.items():
        original = baseline_dynamics[name,]
        _check_frame_dynamic(row["values"], data["cases"][original["case"]], original, api)
    transitions = keyed(document["transition_diagnostics"], ["delta_offset"])
    require(set(transitions) == {(x,) for x in TRANSITION_OFFSETS}, "Incorrect transition-diagnostic coverage")
    for row in transitions.values():
        _check_transition(row, data, api)
    if api is not None:
        for function in ["frame_operators", "frame_static", "frame_evolve", "diagnose_transition"]:
            require(callable(getattr(api, function, None)), "Missing version-2 API: " + function)
        for name, M in frames:
            case = data["cases"][name]
            for rephase in [False, True]:
                answer = api.frame_static(copy.deepcopy(case), M, rephase=rephase)
                _check_frame_static({"mesh": M, **answer}, case, api, encoded=False)
                close(answer["values"]["sewn_chern"], baseline_static[name, M, False]["values"]["chern"], 2e-8, "Reevaluated frame topology")
        for name in FRAME_TRAJECTORIES:
            original = baseline_dynamics[name,]
            case = data["cases"][original["case"]]
            answer = api.frame_evolve(copy.deepcopy(case), original["phi"], original["T"], original["direction"], rtol=1e-10, atol=1e-12)
            _check_frame_dynamic(answer, case, original, api, encoded=False)
        transfer = {key: copy.deepcopy(data["transfer_probe"][key]) for key in ["L", "V", "W", "disorder", "loop"]}
        for case, theta, phi, s, T in [(data["cases"]["pump_L6"], 1.13, data["dynamics"]["fixed_phi"], 0.37, 10.0), (data["cases"]["pump_L8"], 2.27, -0.73, 0.61, 30.0), (transfer, -0.83, 1.39, 0.19, 17.0)]:
            _check_frame_operator(api, case, theta, phi, s, T)
        for offset in [0.485, 0.5, 0.515]:
            answer = api.diagnose_transition(copy.deepcopy(data["cases"]["free_L6"]), offset)
            _check_transition({"delta_offset": offset, **answer}, data, api, encoded=False)
    return {"frame_static": len(frames), "frame_dynamics": len(dynamics), "transition_diagnostics": len(transitions), "extension_apis_checked": api is not None}


def validate(document, inputs_path, solution_path=None):
    check_finite(document)
    schema = load_json(HERE / "results.schema.json")
    Draft202012Validator(schema).validate(document)
    data = load_json(inputs_path)
    require(document["input_sha256"] == hashlib.sha256(Path(inputs_path).read_bytes()).hexdigest(), "Input SHA-256 mismatch")
    require(document["task_sha256"] == hashlib.sha256((HERE.parent / "task.md").read_bytes()).hexdigest(), "Task specification SHA-256 mismatch")
    if solution_path is not None:
        require(document["provenance"]["solution_sha256"] == hashlib.sha256(Path(solution_path).read_bytes()).hexdigest(), "Solution SHA-256 mismatch")
    tol = data["tolerances"]
    static = keyed(document["static"], ["case", "mesh", "rephase"])
    expected_static = {(case, M, rephase) for case in data["cases"] for M in data["static_meshes"] for rephase in [False, True]}
    require(set(static) == expected_static, "Missing or unexpected static record")
    for case, M, rephase in expected_static:
        values = static[case, M, rephase]["values"]
        require(values["sampled_min_gap"] > 0, "Sampled gap must be positive")
        require(0 < values["min_link_overlap"] <= 1 + tol["static_absolute"], "Invalid link overlap")
        require(0 <= values["max_abs_plaquette_phase"] <= math.pi, "Invalid principal phase")
        original = static[case, M, False]["values"]
        for name in values:
            close(values[name], original[name], tol["static_absolute"], "Rephasing invariance: " + name)
    probes = keyed(document["curvature"], ["case"])
    require(set(probes) == {(p["case"],) for p in data["operator_probes"]}, "Incorrect curvature cases")
    for probe in data["operator_probes"]:
        check_coordinates(probes[probe["case"],], probe, ["theta", "phi"])
    lines = keyed(document["adiabatic_lines"], ["case"])
    require(set(lines) == {(name,) for name in data["cases"]}, "Incorrect adiabatic line cases")
    for row in lines.values():
        close(row["phi"], data["dynamics"]["fixed_phi"], label="Adiabatic line twist")

    trajectories = keyed(document["dynamics"], ["id"])
    plan = trajectory_plan(data)
    require(set(trajectories) == {(p["id"],) for p in plan}, "Missing or unexpected trajectory ID")
    for spec in plan:
        row = trajectories[spec["id"],]
        check_coordinates(row, spec, ["case", "phi", "T", "direction"])
        values = row["values"]
        close(values["q_geometric"], values["q_total"] - values["q_equilibrium"], 2e-10, "Charge subtraction")
        L = data["cases"][row["case"]]["L"]
        require(len(values["final_state"]) == math.comb(L, L // 2), "Incorrect final-state dimension")
        norm = sum(real * real + imag * imag for real, imag in values["final_state"])
        close(norm, 1.0, tol["norm_absolute"], "Final-state norm")
        close(values["final_norm_squared"], norm, 2e-10, "Reported norm")
        require(-tol["probability_absolute"] <= values["excitation_probability"] <= 1 + tol["probability_absolute"], "Invalid excitation probability")
    averages = keyed(document["averages"], ["case", "twists", "T", "direction"])
    expected_averages = average_plan(data)
    require(set(averages) == {tuple(p[k] for k in ["case", "twists", "T", "direction"]) for p in expected_averages}, "Missing or unexpected twist average")
    for spec in expected_averages:
        row = averages[spec["case"], spec["twists"], spec["T"], spec["direction"]]
        group = []
        for b in range(spec["twists"]):
            phi = 2 * math.pi * (b + 0.5) / spec["twists"]
            matches = [r for r in document["dynamics"] if r["case"] == spec["case"] and r["T"] == spec["T"] and r["direction"] == spec["direction"] and abs(r["phi"] - phi) <= 1e-12]
            require(len(matches) == 1, "Average has missing or ambiguous midpoint trajectory")
            group.append(matches[0]["values"])
        for key in CHARGES:
            close(row[key], sum(v[key] for v in group) / len(group), 2e-10, "Twist average: " + key)
        close(row["q_geometric"], row["q_total"] - row["q_equilibrium"], 2e-10, "Average charge subtraction")
    extra = _validate_extensions(document, data, _load_api(solution_path) if solution_path is not None else None)
    return {"passed": True, "static": len(static), "curvature": len(probes), "adiabatic_lines": len(lines), "dynamics": len(trajectories), "averages": len(averages), **extra}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("results", type=Path, nargs="?")
    parser.add_argument("--inputs", type=Path, default=HERE / "inputs.json")
    parser.add_argument("--solution", type=Path)
    parser.add_argument("--list-plan", action="store_true")
    args = parser.parse_args()
    try:
        if args.list_plan:
            print(json.dumps(trajectory_plan(load_json(args.inputs)), indent=2))
            return
        require(args.results is not None, "Provide a results.json path")
        print(json.dumps(validate(load_json(args.results), args.inputs, args.solution), indent=2))
    except Exception as exc:
        print(json.dumps({"passed": False, "error": str(exc)}))
        raise SystemExit(1)


if __name__ == "__main__":
    main()
