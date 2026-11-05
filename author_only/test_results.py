"""Regression checks for the new artifact contract and its numerical gate."""

import copy
import json
import math
import tempfile
import unittest
from pathlib import Path

import grade
import reference
from export_reference_artifact import artifact, convert


class ResultArtifactTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.gold = artifact.load_json(Path(__file__).with_name("reference_results.json"))
        cls.solution = Path(reference.__file__)
        cls.data = artifact.load_json(reference.INPUTS)
        cls.document = convert(cls.gold, reference.INPUTS, cls.solution)

    def validate(self, document):
        return artifact.validate(document, reference.INPUTS, self.solution)

    def author_check(self, document):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "results.json"
            path.write_text(json.dumps(document, allow_nan=False))
            return grade.check_artifact(path, self.solution, self.data, self.gold)

    def test_reference_artifact_passes_both_validators(self):
        public = self.validate(self.document)
        self.assertEqual(public["dynamics"], 71)
        self.assertEqual(public["averages"], 6)
        result = self.author_check(self.document)
        self.assertTrue(result["passed"], result)
        self.assertEqual(result["numerical_checks"], 95)

    def test_plan_matches_independent_reference_schedule(self):
        actual = [(r["case"], r["phi"], r["T"], r["direction"]) for r in artifact.trajectory_plan(self.data)]
        self.assertEqual(actual, reference.trajectory_specs(self.data))

    def test_record_order_is_irrelevant(self):
        document = copy.deepcopy(self.document)
        for key in ["static", "curvature", "adiabatic_lines", "dynamics", "averages"]:
            document[key].reverse()
        self.assertTrue(self.author_check(document)["passed"])

    def test_global_phase_is_irrelevant(self):
        document = copy.deepcopy(self.document)
        for row in document["dynamics"]:
            row["values"]["final_state"] = [[-imag, real] for real, imag in row["values"]["final_state"]]
        self.assertTrue(self.author_check(document)["passed"])

    def test_missing_results_rejected_by_author_gate(self):
        result = grade.check_artifact(Path("/tmp/no-such-pump-artifact/results.json"), self.solution, self.data, self.gold)
        self.assertFalse(result["passed"])

    def test_coverage_types_shapes_hashes_and_consistency(self):
        def mutate(document, fault):
            if fault == "missing_trajectory": document["dynamics"].pop()
            if fault == "duplicate_trajectory": document["dynamics"][1] = copy.deepcopy(document["dynamics"][0])
            if fault == "wrong_coordinate": document["dynamics"][0]["phi"] += 0.1
            if fault == "wrong_curvature_coordinate": document["curvature"][0]["theta"] += 0.1
            if fault == "duplicate_curvature_case": document["curvature"][1] = copy.deepcopy(document["curvature"][0])
            if fault == "wrong_state_dimension": document["dynamics"][0]["values"]["final_state"].pop()
            if fault == "complex_as_string": document["dynamics"][0]["values"]["final_state"][0] = "1+0j"
            if fault == "boolean_charge": document["dynamics"][0]["values"]["q_total"] = True
            if fault == "null_charge": document["dynamics"][0]["values"]["q_total"] = None
            if fault == "nonfinite_charge": document["dynamics"][0]["values"]["q_total"] = math.inf
            if fault == "wrong_subtraction": document["dynamics"][0]["values"]["q_geometric"] += 0.1
            if fault == "wrong_norm": document["dynamics"][0]["values"]["final_norm_squared"] += 0.1
            if fault == "wrong_average": document["averages"][0]["q_total"] += 0.1
            if fault == "wrong_input_hash": document["input_sha256"] = "0" * 64
            if fault == "wrong_source_hash": document["provenance"]["solution_sha256"] = "0" * 64
            if fault == "unknown_field": document["unexpected"] = "extra"
            if fault == "broken_rephasing": document["static"][1]["values"]["chern"] += 0.1

        faults = ["missing_trajectory", "duplicate_trajectory", "wrong_coordinate", "wrong_curvature_coordinate", "duplicate_curvature_case", "wrong_state_dimension", "complex_as_string", "boolean_charge", "null_charge", "nonfinite_charge", "wrong_subtraction", "wrong_norm", "wrong_average", "wrong_input_hash", "wrong_source_hash", "unknown_field", "broken_rephasing"]
        for fault in faults:
            with self.subTest(fault=fault):
                document = copy.deepcopy(self.document)
                mutate(document, fault)
                with self.assertRaises(Exception):
                    self.validate(document)

    def test_consistent_but_incorrect_physics_is_rejected(self):
        document = copy.deepcopy(self.document)
        values = document["dynamics"][0]["values"]  # Fixed twist; no average uses this row.
        values["q_total"] += 0.01
        values["q_geometric"] += 0.01
        self.assertTrue(self.validate(document)["passed"])
        self.assertFalse(self.author_check(document)["passed"])

    def test_strict_json_parser(self):
        for text in ['{"a": 1, "a": 2}', '{"a": NaN}', '{"a": Infinity}']:
            with self.subTest(text=text), tempfile.TemporaryDirectory() as directory:
                path = Path(directory) / "invalid.json"
                path.write_text(text)
                with self.assertRaises(ValueError):
                    artifact.load_json(path)


if __name__ == "__main__":
    unittest.main(verbosity=2)
