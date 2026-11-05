"""Deterministic numerical grader; keep outside the candidate distribution."""

from __future__ import annotations

import argparse
import copy
import hashlib
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np

import reference


def load_module(path):
    spec = importlib.util.spec_from_file_location("candidate_solution", path.resolve())
    module = importlib.util.module_from_spec(spec)
    sys.path.insert(0, str(path.resolve().parent))
    spec.loader.exec_module(module)
    return module


class Checks:
    def __init__(self):
        self.rows = []

    def run(self, group, name, operation):
        try:
            operation()
            self.rows.append({"group": group, "name": name, "passed": True})
        except Exception as exc:
            self.rows.append({"group": group, "name": name, "passed": False, "error": str(exc)})


def number(value):
    if isinstance(value, (bool, complex, np.complexfloating)):
        raise AssertionError("Expected a finite real scalar")
    result = float(value)
    assert np.isfinite(result), "Nonfinite scalar"
    return result


def close(actual, expected, tolerance):
    actual = number(actual)
    assert abs(actual - expected) <= tolerance, f"got {actual:.12g}; expected {expected:.12g}; atol {tolerance}"


def transfer_cases(data):
    """Fixed rational parameters, not random sampling or wall-clock seeds."""
    probes = []
    for k in range(3):
        case = copy.deepcopy(data["transfer_probe"])
        theta, phi = case.pop("theta"), case.pop("phi")
        case["V"] = (6 + k) / 10
        case["W"] = (21 + k) / 100
        case["disorder"][0] = (3 + k) / 20
        probes.append((case, theta + k / 7, phi - k / 11))
    return probes


def check_operators(candidate, case, theta, phi, tolerance):
    expected = reference.build_operators(case, theta, phi)
    got = candidate.build_operators(copy.deepcopy(case), theta, phi)
    for key in ("H", "dtheta", "dphi"):
        actual = np.asarray(got[key], dtype=complex)
        assert actual.shape == expected[key].shape, f"Wrong shape for {key}"
        assert np.all(np.isfinite(actual)), f"Nonfinite {key}"
        np.testing.assert_allclose(actual, expected[key], atol=tolerance["operators_absolute"], rtol=tolerance["operators_relative"], err_msg=key)


def check_static(candidate, case, entry, tolerance, rephase):
    got = candidate.static_observables(copy.deepcopy(case), entry["mesh"], rephase=rephase)
    for key, value in entry["values"].items():
        close(got[key], value, tolerance["static_absolute"])


def check_evolution(candidate, case, entry, tolerance):
    got = candidate.evolve(copy.deepcopy(case), entry["phi"], entry["T"], entry["direction"], rtol=1e-10, atol=1e-12)
    expected = entry["values"]
    for key in ("q_total", "q_equilibrium", "q_geometric"):
        close(got[key], expected[key], tolerance["charge_absolute"])
    close(got["q_geometric"], number(got["q_total"]) - number(got["q_equilibrium"]), 2e-10)
    close(got["excitation_probability"], expected["excitation_probability"], tolerance["probability_absolute"])
    close(got["final_norm_squared"], 1.0, tolerance["norm_absolute"])
    state = np.asarray(got["final_state"], dtype=complex)
    gold = np.array([complex(*z) for z in expected["final_state"]])
    assert state.shape == gold.shape and np.all(np.isfinite(state)), "Invalid final state"
    norm = float(np.vdot(state, state).real)
    close(norm, 1.0, tolerance["norm_absolute"])
    close(got["final_norm_squared"], norm, 2e-10)
    infidelity = 1 - abs(np.vdot(gold, state)) ** 2 / (np.vdot(gold, gold).real * norm)
    assert infidelity <= tolerance["final_state_infidelity"], f"Final-state infidelity {infidelity}"
    H0 = reference.build_operators(case, 0.0, entry["phi"])["H"]
    _, vectors = reference.eigh(H0, subset_by_index=[0, 0])
    excitation = 1 - abs(np.vdot(vectors[:, 0], state)) ** 2 / norm
    close(got["excitation_probability"], excitation, tolerance["probability_absolute"])


def grade(candidate, data, gold, selected="all"):
    assert gold["metadata"]["input_sha256"] == hashlib.sha256(reference.INPUTS.read_bytes()).hexdigest(), "Reference fixture does not match inputs.json"
    checks = Checks()
    tolerance = data["tolerances"]
    if selected in ("all", "operators"):
        for n, probe in enumerate(data["operator_probes"]):
            case = data["cases"][probe["case"]]
            checks.run("operators", f"public-{n}", lambda c=case, p=probe: check_operators(candidate, c, p["theta"], p["phi"], tolerance))
        for n, (case, theta, phi) in enumerate(transfer_cases(data)):
            checks.run("operators", f"transfer-{n}", lambda c=case, t=theta, p=phi: check_operators(candidate, c, t, p, tolerance))
    if selected in ("all", "static"):
        for entry in gold["static"]:
            for rephase in (False, True):
                case = data["cases"][entry["case"]]
                label = f"{entry['case']}-M{entry['mesh']}-rephase{rephase}"
                checks.run("static", label, lambda c=case, e=entry, r=rephase: check_static(candidate, c, e, tolerance, r))
    if selected in ("all", "curvature"):
        probes = [(data["cases"][p["case"]], p["theta"], p["phi"], p["value"]) for p in gold["curvature"]]
        probes += [(c, t, p, reference.kubo_curvature(c, t, p)) for c, t, p in transfer_cases(data)]
        for n, (case, theta, phi, value) in enumerate(probes):
            checks.run("curvature", f"point-{n}", lambda c=case, t=theta, p=phi, v=value: close(candidate.kubo_curvature(copy.deepcopy(c), t, p), v, tolerance["curvature_absolute"]))
        for entry in gold["adiabatic_lines"]:
            case = data["cases"][entry["case"]]
            checks.run("curvature", f"adiabatic-line-{entry['case']}", lambda c=case, e=entry: close(reference.adiabatic_charge(copy.deepcopy(c), e["phi"], curvature_function=candidate.kubo_curvature), e["value"], 2e-7))
    if selected in ("all", "dynamics"):
        for n, entry in enumerate(gold["dynamics"]):
            case = data["cases"][entry["case"]]
            label = f"trajectory-{n}-{entry['case']}-T{entry['T']}-d{entry['direction']}"
            checks.run("dynamics", label, lambda c=case, e=entry: check_evolution(candidate, c, e, tolerance))
        case, _, _ = transfer_cases(data)[1]
        entry = {"phi": 0.47, "T": 17.5, "direction": -1}
        entry["values"] = reference.serializable_trajectory(reference.evolve(case, **entry, rtol=2e-12, atol=2e-14, eq_nodes=256))
        checks.run("dynamics", "transfer-L4", lambda: check_evolution(candidate, case, entry, tolerance))
    weights = {"operators": 25, "static": 25, "curvature": 10, "dynamics": 40}
    groups = {}
    for group, weight in weights.items():
        rows = [x for x in checks.rows if x["group"] == group]
        if rows:
            passed = sum(x["passed"] for x in rows)
            groups[group] = {"passed": passed, "total": len(rows), "points": weight * passed / len(rows), "available": weight}
    return {"passed": all(x["passed"] for x in checks.rows), "score": sum(x["points"] for x in groups.values()), "available": sum(x["available"] for x in groups.values()), "groups": groups, "checks": checks.rows}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--submission", type=Path, required=True)
    parser.add_argument("--reference", type=Path, default=Path(__file__).with_name("reference_results.json"))
    parser.add_argument("--group", choices=["all", "operators", "static", "curvature", "dynamics"], default="all")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    candidate = load_module(args.submission)
    result = grade(candidate, json.loads(reference.INPUTS.read_text()), json.loads(args.reference.read_text()), args.group)
    rendered = json.dumps(result, indent=2, allow_nan=False)
    if args.output:
        args.output.write_text(rendered + "\n")
    print(rendered)
    raise SystemExit(0 if result["passed"] else 1)


if __name__ == "__main__":
    main()
