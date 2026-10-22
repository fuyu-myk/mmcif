import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import run

CAT = [sys.executable, "-c", "import sys; sys.stdout.write(open(sys.argv[1]).read())"]
SLEEP = [sys.executable, "-c", "import time; time.sleep(5)"]


def dump(entry: str = "1ABC", models: int = 1) -> str:
    return json.dumps({
        "entry": entry,
        "models": [{"number": number, "atoms": [], "residues": []} for number in range(1, models + 1)],
        "entities": [],
        "assemblies": [],
    })


def write_rules(text: str) -> Path:
    handle = tempfile.NamedTemporaryFile("w", suffix=".toml", delete=False)
    handle.write(text)
    handle.close()
    return Path(handle.name)


class RulesTest(unittest.TestCase):
    def test_an_entry_rule_and_an_entry_free_rule_load(self):
        rules = run.load_rules(write_rules(
            '[[explained]]\nentry = "1M4X"\ncategory = "assembly_set"\nreason = "r"\n'
            '[[explained]]\nexample = "2BUK"\ncategory = "assembly_set"\nmatch = "only in ours"\nreason = "r"\n'
        ))
        self.assertEqual(rules[0], run.Rule("1m4x", None, "assembly_set", None, "r"))
        self.assertEqual(rules[1], run.Rule(None, "2buk", "assembly_set", "only in ours", "r"))

    def test_a_rule_names_either_an_entry_or_an_example(self):
        for table in (
            'category = "c"\nreason = "r"\n',
            'entry = "1abc"\nexample = "1abc"\ncategory = "c"\nreason = "r"\n',
            'entry = "1abc"\ncategory = "c"\nreason = "r"\nnote = "n"\n',
            'entry = "1abc"\nreason = "r"\n',
            'entry = ""\ncategory = "c"\nreason = "r"\n',
            'entry = "1abc"\ncategory = "c"\nmatch = ""\nreason = "r"\n',
            'example = "1abc"\ncategory = "c"\nreason = "r"\n',
        ):
            with self.assertRaises(ValueError, msg=table):
                run.load_rules(write_rules("[[explained]]\n" + table))

    def test_an_entry_free_rule_covers_every_entry(self):
        rule = run.Rule(None, "2buk", "assembly_set", "only in ours", "r")
        self.assertTrue(run.rule_covers(rule, "9xyz", "assembly_set", "assembly 2 only in ours"))
        self.assertFalse(run.rule_covers(rule, "9xyz", "assembly_set", "assembly 2 only in gemmi"))
        self.assertFalse(run.rule_covers(rule, "9xyz", "assembly_copies", "assembly 2 only in ours"))
        entry_rule = run.Rule("1m4x", None, "assembly_set", None, "r")
        self.assertFalse(run.rule_covers(entry_rule, "9xyz", "assembly_set", "anything"))

    def test_explained_differences_are_counted_per_rule(self):
        rules = [
            run.Rule(None, "2buk", "assembly_set", None, "r"),
            run.Rule("1abc", None, "assembly_set", "2", "r"),
        ]
        result = run.apply_rules("1abc", {"assembly_set": ["assembly 2", "assembly 3"], "coordinate": ["x"]}, rules)
        self.assertEqual(result.unexplained, {"coordinate": ["x"]})
        self.assertEqual(result.explained_count, 2)
        self.assertEqual(result.explained_by_rule, {0: 2, 1: 1})


class BinaryDifferencesTest(unittest.TestCase):
    def test_binary_differences_are_named_for_the_two_files(self):
        renamed = run.as_bcif_differences({
            "atoms_only_ours": ["model 1 A/1/N"],
            "assembly_set": ["assembly 2 only in gemmi"],
            "coordinate": ["N: [1, 2, 3] vs [1, 2, 4]"],
        })
        self.assertEqual(renamed, {
            "bcif_atoms_only_bcif": ["model 1 A/1/N"],
            "bcif_assembly_set": ["assembly 2 only in cif"],
            "bcif_coordinate": ["N: [1, 2, 3] vs [1, 2, 4]"],
        })


class ProcessEntryTest(unittest.TestCase):
    def setUp(self):
        self.root = Path(tempfile.mkdtemp())
        (self.root / "data").mkdir()
        (self.root / "out").mkdir()
        patches = (
            mock.patch.object(run, "REPO_ROOT", self.root),
            mock.patch.object(run, "DATA_DIR", self.root / "data"),
            mock.patch.object(run, "OUT_DIR", self.root / "out"),
            mock.patch.object(run, "GEMMI_DUMP_COMMAND", CAT),
        )
        for patch in patches:
            patch.start()
            self.addCleanup(patch.stop)

    def entry(self, cif: str, bcif: str) -> None:
        (self.root / "data" / "1abc.cif").write_text(cif)
        (self.root / "data" / "1abc.bcif").write_text(bcif)

    def test_a_clean_entry_leaves_no_dumps(self):
        self.entry(dump(), dump())
        result = run.process_entry("1abc", CAT, 1e-3, [], 30)
        self.assertEqual(result.unexplained, {})
        self.assertEqual(list((self.root / "out").iterdir()), [])

    def test_a_binary_difference_is_reported_and_its_dumps_kept(self):
        self.entry(dump(), dump(models=2))
        result = run.process_entry("1abc", CAT, 1e-3, [], 30)
        self.assertEqual(result.unexplained, {"bcif_models": ["model 2 only in bcif"]})
        self.assertEqual(
            sorted(path.name for path in (self.root / "out").iterdir()),
            ["1abc.gemmi.json", "1abc.ours.json", "1abc.ours_bcif.json"],
        )

    def test_a_dumper_that_runs_too_long_is_reported_as_a_timeout(self):
        self.entry(dump(), dump())
        result = run.process_entry("1abc", SLEEP, 1e-3, [], 0.5)
        self.assertEqual(set(result.unexplained), {"ours_timeout", "ours_bcif_timeout"})
        self.assertEqual(result.unexplained["ours_timeout"], ["no dump within 0.5 s"])

    def test_a_timeout_keeps_what_the_dumper_wrote_to_stderr(self):
        self.entry(dump(), dump())
        noisy = [sys.executable, "-c", "import sys, time; print('working', file=sys.stderr, flush=True); time.sleep(5)"]
        result = run.process_entry("1abc", noisy, 1e-3, [], 0.5)
        self.assertEqual(result.unexplained["ours_timeout"], ["no dump within 0.5 s: working"])

    def test_a_missing_dumper_and_an_unreadable_dump_are_reported(self):
        self.entry(dump(), dump())
        result = run.process_entry("1abc", [str(self.root / "no-such-dumper")], 1e-3, [], 30)
        self.assertEqual(set(result.unexplained), {"ours_failed", "ours_bcif_failed"})
        binary = [sys.executable, "-c", "import sys; sys.stdout.buffer.write(b'\\xff\\xfe')"]
        result = run.process_entry("1abc", binary, 1e-3, [], 30)
        self.assertEqual(set(result.unexplained), {"ours_invalid_dump", "ours_bcif_invalid_dump"})

    def test_failures_can_be_explained_by_a_rule(self):
        self.entry(dump(), dump())
        rule = run.Rule("1abc", None, "gemmi_failed", None, "r")
        with mock.patch.object(run, "GEMMI_DUMP_COMMAND", [sys.executable, "-c", "raise SystemExit(3)"]):
            result = run.process_entry("1abc", CAT, 1e-3, [rule], 30)
        self.assertEqual(result.unexplained, {})
        self.assertEqual(result.explained_by_rule, {0: 1})


class SelectionTest(unittest.TestCase):
    def test_the_full_corpus_is_every_pinned_entry(self):
        corpus = write_rules(json.dumps({"pins": {"2abc": {}, "1abc": {}}}))
        manifest = write_rules(json.dumps({"9abc": ["baseline"]}))
        with mock.patch.object(run, "CORPUS_PATH", corpus), mock.patch.object(run, "MANIFEST_PATH", manifest):
            self.assertEqual(run.select_entries("full", None), ["1abc", "2abc"])
            self.assertEqual(run.select_entries("seed", None), ["9abc"])
            self.assertEqual(run.select_entries("full", ["2ABC"]), ["2abc"])
            with self.assertRaises(SystemExit):
                run.select_entries("full", ["9abc"])


class ResultsTest(unittest.TestCase):
    def test_a_subset_run_does_not_replace_the_corpus_results(self):
        self.assertEqual(run.default_results_path("full", None).name, "results-full.json")
        self.assertEqual(run.default_results_path("full", ["1abc"]).name, "results-full-subset.json")

    def test_the_results_document_is_complete_and_bounded(self):
        rules = [run.Rule(None, "2buk", "assembly_set", None, "why")]
        results = [
            run.EntryResult("2abc", {"coordinate": [str(n) for n in range(15)]}, 0, {}),
            run.EntryResult("1abc", {}, 3, {0: 3}),
            run.EntryResult("3abc", {}, 0, {}),
        ]
        document = run.results_document(results, rules, {"corpus": "full"})
        self.assertEqual(document["corpus"], "full")
        self.assertEqual(list(document["entries"]), ["1abc", "2abc", "3abc"])
        self.assertEqual(document["entries"]["3abc"], {})
        self.assertEqual(document["entries"]["1abc"], {"explained": 3})
        coordinate = document["entries"]["2abc"]["unexplained"]["coordinate"]
        self.assertEqual(coordinate["count"], 15)
        self.assertEqual(len(coordinate["examples"]), run.diff.EXAMPLES_PER_CATEGORY)
        self.assertEqual(document["rules"], [{
            "example": "2buk",
            "category": "assembly_set",
            "reason": "why",
            "explained": 3,
            "entries": ["1abc"],
        }])


if __name__ == "__main__":
    unittest.main()
