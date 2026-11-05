"""Convert the author-only oracle fixture to the candidate results format."""

import argparse
import copy
import hashlib
import json
import math
import sys
from pathlib import Path

import reference

sys.path.insert(0, str(reference.ROOT / "benchmark"))
import validate_results as artifact


def convert(gold, inputs_path, solution_path):
    data = artifact.load_json(inputs_path)
    require_hash = hashlib.sha256(Path(inputs_path).read_bytes()).hexdigest()
    if gold["metadata"]["input_sha256"] != require_hash:
        raise ValueError("Oracle inputs do not match the requested artifact")
    plan = artifact.trajectory_plan(data)
    if len(plan) != len(gold["dynamics"]):
        raise ValueError("Oracle trajectory count does not match the public plan")
    document = {
        "schema_version": 1,
        "input_sha256": require_hash,
        "provenance": {
            **{key: gold["metadata"][key] for key in ["python", "numpy", "scipy"]},
            "command": "python3 author_only/export_reference_artifact.py --output /tmp/reference-artifact.json",
            "solution_sha256": hashlib.sha256(Path(solution_path).read_bytes()).hexdigest(),
        },
        "static": [{**copy.deepcopy(row), "rephase": rephase} for row in gold["static"] for rephase in [False, True]],
        "curvature": copy.deepcopy(gold["curvature"]),
        "adiabatic_lines": copy.deepcopy(gold["adiabatic_lines"]),
        "dynamics": [],
        "averages": [],
    }
    for spec, row in zip(plan, gold["dynamics"]):
        artifact.check_coordinates(row, spec, ["case", "phi", "T", "direction"])
        document["dynamics"].append({"id": spec["id"], **copy.deepcopy(row)})
    for spec in artifact.average_plan(data):
        nodes = [2 * math.pi * (b + 0.5) / spec["twists"] for b in range(spec["twists"])]
        group = [row["values"] for row in document["dynamics"] if row["case"] == spec["case"] and row["T"] == spec["T"] and row["direction"] == spec["direction"] and any(abs(row["phi"] - phi) <= 1e-12 for phi in nodes)]
        artifact.require(len(group) == spec["twists"], "Incorrect oracle average membership")
        document["averages"].append({**spec, **{key: sum(row[key] for row in group) / len(group) for key in artifact.CHARGES}})
    artifact.validate(document, inputs_path, solution_path)
    return document


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", type=Path, default=Path(__file__).with_name("reference_results.json"))
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    document = convert(artifact.load_json(args.reference), reference.INPUTS, Path(reference.__file__))
    document["provenance"]["command"] = "python3 " + " ".join(sys.argv)
    args.output.write_text(json.dumps(document, indent=2, allow_nan=False) + "\n")
    print(json.dumps(artifact.validate(document, reference.INPUTS, Path(reference.__file__)), indent=2))


if __name__ == "__main__":
    main()
