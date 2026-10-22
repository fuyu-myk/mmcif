import hashlib
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import fetch

CORPUS = {
    "baseline": ["2ABC", "1abc"],
    "pools": {"microheterogeneity": {"ids": ["3abc", "1abc"]}},
    "pins": {"9abc": {}, "2abc": {}},
}


class IdsTest(unittest.TestCase):
    def test_candidates_are_the_baseline_and_pools(self):
        self.assertEqual(fetch.candidate_ids(CORPUS), ["1abc", "2abc", "3abc"])

    def test_a_corpus_fetch_also_covers_pinned_entries(self):
        self.assertEqual(fetch.corpus_ids(CORPUS), ["1abc", "2abc", "3abc", "9abc"])


class ChangedFilesTest(unittest.TestCase):
    def test_only_files_whose_hash_differs_from_the_pin_are_reported(self):
        with tempfile.TemporaryDirectory() as directory, mock.patch.object(fetch, "DATA_DIR", Path(directory)):
            Path(directory, "1abc.cif").write_text("data_1abc\n")
            Path(directory, "1abc.bcif").write_bytes(b"binary")
            Path(directory, "2abc.cif").write_text("data_2abc revised\n")
            pins = {
                "1abc": {"cif": hashlib.sha256(b"data_1abc\n").hexdigest(), "bcif": hashlib.sha256(b"binary").hexdigest()},
                "2abc": {"cif": hashlib.sha256(b"data_2abc\n").hexdigest(), "bcif": "absent"},
            }
            self.assertEqual(fetch.sha256_of(Path(directory, "1abc.bcif")), pins["1abc"]["bcif"])
            self.assertEqual(fetch.changed_files(pins), ["2abc.cif"])


if __name__ == "__main__":
    unittest.main()
