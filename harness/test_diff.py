import unittest

import diff


def dump(b_factor: float | None) -> dict:
    atom = ["A", "A", 1, "", "GLY", "N", "", "N", 0.0, 0.0, 0.0, 1.0, b_factor]
    return {
        "entry": "1ABC",
        "models": [{"number": 1, "atoms": [atom], "residues": []}],
        "entities": [],
        "assemblies": [],
    }


class MissingValueTest(unittest.TestCase):
    def test_a_missing_b_factor_differs_from_a_number(self):
        self.assertEqual(
            diff.compare(dump(None), dump(20.0)),
            {"b_iso": ["model 1: A/A GLY 1 N: None vs 20.0"]},
        )

    def test_two_missing_b_factors_agree(self):
        self.assertEqual(diff.compare(dump(None), dump(None)), {})

    def test_numbers_within_the_tolerance_agree(self):
        self.assertEqual(diff.compare(dump(20.0), dump(20.0005)), {})
        self.assertIn("b_iso", diff.compare(dump(20.0), dump(20.5)))


if __name__ == "__main__":
    unittest.main()
