import copy
import unittest

import report

CORPUS = {
    "random_seed": 7,
    "holdings": {"count": 261546, "sha256": "abc123"},
    "baseline": ["1abc", "2abc"],
    "pools": {"microheterogeneity": {"source": "https://pdbj.example/sql", "query": "SELECT pdbid", "ids": ["3abc"]}},
    "class_size": 50,
    "entries": {"1abc": {"altlocs": 2}, "2abc": {}, "3abc": {"microheterogeneity": 1, "altlocs": 1}},
    "classes": {"altlocs": ["9abc", "1abc"], "microheterogeneity": ["3abc"]},
    "pins": {"1abc": {}, "2abc": {}, "3abc": {}, "9abc": {}},
}
SEED = {"9abc": ["altlocs", "baseline"]}
RESULTS = {
    "corpus": "full",
    "subset": False,
    "source": "v1-2-gabc",
    "gemmi": "0.7.5",
    "tolerance": 0.001,
    "timeout_seconds": 300.0,
    "entries": {
        "1abc": {},
        "2abc": {"explained": 3},
        "3abc": {"unexplained": {"b_iso": {"count": 2400, "examples": ["model 1: A/1 N: 0.0 vs 20.0", "second", "third", "fourth"]}}},
        "9abc": {},
    },
    "rules": [
        {"example": "2abc", "category": "assembly_set", "match": "only in ours", "reason": "gemmi drops them", "explained": 3, "entries": ["2abc"]},
        {"entry": "4xyz", "category": "models", "reason": "unused", "explained": 0, "entries": []},
    ],
}


class ClassRowsTest(unittest.TestCase):
    def test_rows_count_outcomes_per_class(self):
        rows = dict(report.class_rows(RESULTS, CORPUS, SEED))
        self.assertEqual(rows["Random baseline"], {"entries": 2, "agree": 1, "explained": 1, "unexplained": 0})
        self.assertEqual(rows["Hand-picked seed"], {"entries": 1, "agree": 1, "explained": 0, "unexplained": 0})
        self.assertEqual(rows["Alternate conformations"], {"entries": 3, "agree": 2, "explained": 0, "unexplained": 1})
        self.assertEqual(rows["Microheterogeneity"], {"entries": 1, "agree": 0, "explained": 0, "unexplained": 1})
        self.assertEqual(rows["Insertion codes"], {"entries": 0, "agree": 0, "explained": 0, "unexplained": 0})
        self.assertEqual(rows["All entries"], {"entries": 4, "agree": 2, "explained": 1, "unexplained": 1})
        self.assertEqual(list(rows)[-1], "All entries")


class RenderTest(unittest.TestCase):
    def test_the_report_states_the_corpus_the_rules_and_every_open_difference(self):
        text = report.render(RESULTS, CORPUS, SEED)
        for expected in (
            "4 entries",
            "seed 7",
            "261,546",
            "`abc123`",
            "gemmi 0.7.5",
            "`v1-2-gabc`",
            "| Alternate conformations | 3 | 2 | 0 | 1 |",
            "gemmi drops them",
            "`only in ours`",
            "3 differences in 1 entry: 2abc",
            "**3abc**: b_iso (2,400), for example `model 1: A/1 N: 0.0 vs 20.0`",
        ):
            self.assertIn(expected, text)
        self.assertNotIn("fourth", text)
        self.assertNotIn("unused", text)
        self.assertTrue(text.endswith("\n"))
        self.assertNotIn("—", text)

    def test_a_clean_run_says_so(self):
        clean = copy.deepcopy(RESULTS)
        clean["entries"]["3abc"] = {}
        self.assertIn("No differences are unexplained.", report.render(clean, CORPUS, SEED))


class CheckTest(unittest.TestCase):
    def test_only_a_complete_full_run_can_be_reported(self):
        report.check_results(RESULTS, CORPUS)
        for change in ({"subset": True}, {"corpus": "seed"}, {"entries": {"1abc": {}}}):
            with self.assertRaises(ValueError, msg=change):
                report.check_results(RESULTS | change, CORPUS)


if __name__ == "__main__":
    unittest.main()
