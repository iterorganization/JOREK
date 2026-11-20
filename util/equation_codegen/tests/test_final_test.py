"""Checks for the discrepancy benchmark's comparison logic."""

import sys
import unittest
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))
sys.path.insert(0, str(PROJECT_ROOT / "src"))
sys.path.insert(0, str(PROJECT_ROOT / "examples"))

from final_test import REFERENCE, block_name, compare  # noqa: E402


def reference(blocks):
    return {"blocks": blocks}


class CompareTest(unittest.TestCase):
    def test_identical_discrepancies_pass(self):
        block = {"findings": [7], "source_only": ["a*b"], "generated_only": ["c*d"]}
        self.assertEqual(
            compare(reference({"x": block}), {"x": dict(block)}), [],
        )

    def test_a_discrepancy_that_disappeared_is_reported(self):
        stored = {"findings": [7], "source_only": ["a*b"], "generated_only": []}
        current = {"source_only": [], "generated_only": []}
        problems = compare(reference({"x": stored}), {"x": current})
        self.assertTrue(any("- a*b" in line for line in problems))

    def test_a_new_discrepancy_is_reported(self):
        stored = {"findings": [7], "source_only": [], "generated_only": []}
        current = {"source_only": ["new*term"], "generated_only": []}
        problems = compare(reference({"x": stored}), {"x": current})
        self.assertTrue(any("+ new*term" in line for line in problems))

    def test_a_block_that_now_agrees_is_reported(self):
        stored = {"findings": [7], "source_only": ["a*b"], "generated_only": []}
        problems = compare(reference({"x": stored}), {})
        self.assertTrue(any("BLOCK NOW AGREES" in line for line in problems))

    def test_a_block_that_started_to_disagree_is_reported(self):
        current = {"source_only": ["a*b"], "generated_only": []}
        problems = compare(reference({}), {"x": current})
        self.assertTrue(any("NEW BLOCK" in line for line in problems))

    def test_repeated_monomials_are_compared_with_multiplicity(self):
        stored = {"findings": [], "source_only": ["a*b", "a*b"], "generated_only": []}
        current = {"source_only": ["a*b"], "generated_only": []}
        self.assertTrue(compare(reference({"x": stored}), {"x": current}))

    def test_block_name_is_stable(self):
        self.assertEqual(
            block_name(("Two-temperature (Ti/Te) model", "amat(var_u,var_u)", 0)),
            "Two-temperature (Ti/Te) model | amat(var_u,var_u) | #0",
        )


class ReferenceTest(unittest.TestCase):
    """The stored reference must stay well formed and fully annotated."""

    def test_reference_is_complete(self):
        import json

        document = json.loads(REFERENCE.read_text(encoding="utf-8"))
        blocks = document["blocks"]
        self.assertEqual(document["block_count"], len(blocks))
        self.assertEqual(
            document["source_only_lines"],
            sum(len(item["source_only"]) for item in blocks.values()),
        )
        self.assertEqual(
            document["generated_only_lines"],
            sum(len(item["generated_only"]) for item in blocks.values()),
        )
        for name, item in blocks.items():
            with self.subTest(block=name):
                self.assertTrue(
                    item["source_only"] or item["generated_only"],
                    "a recorded block must actually disagree",
                )
                self.assertTrue(
                    item["findings"],
                    "every recorded block must name the finding it belongs to",
                )
                self.assertEqual(item["source_only"], sorted(item["source_only"]))
                self.assertEqual(
                    item["generated_only"], sorted(item["generated_only"]),
                )


if __name__ == "__main__":
    unittest.main()
