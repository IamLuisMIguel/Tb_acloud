"""Public format and consistency checks; contains no reference answers."""

import argparse
import hashlib
import json
import math
from pathlib import Path

from jsonschema import Draft202012Validator

HERE = Path(__file__).resolve().parent
CHARGES = ("q_total", "q_equilibrium", "q_geometric")


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


def validate(document, inputs_path, solution_path=None):
    check_finite(document)
    schema = load_json(HERE / "results.schema.json")
    Draft202012Validator(schema).validate(document)
    data = load_json(inputs_path)
    require(document["input_sha256"] == hashlib.sha256(Path(inputs_path).read_bytes()).hexdigest(), "Input SHA-256 mismatch")
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
    return {"passed": True, "static": len(static), "curvature": len(probes), "adiabatic_lines": len(lines), "dynamics": len(trajectories), "averages": len(averages)}


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
