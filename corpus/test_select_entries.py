import tempfile
import unittest
from pathlib import Path
from unittest import mock

from gemmi import cif

import classify
import select_entries

ATOM_SITE_TAGS = (
    "id",
    "label_asym_id",
    "label_seq_id",
    "label_comp_id",
    "label_alt_id",
    "auth_asym_id",
    "auth_seq_id",
    "pdbx_PDB_ins_code",
    "pdbx_PDB_model_num",
)


def atom_site(*rows: tuple[str, ...]) -> str:
    header = "loop_\n" + "".join(f"_atom_site.{tag}\n" for tag in ATOM_SITE_TAGS)
    body = "".join(f"{index + 1} {' '.join(row)}\n" for index, row in enumerate(rows))
    return header + body


def block(*sections: str) -> cif.Block:
    return cif.read_string("data_test\n" + "".join(sections)).sole_block()


def residue(asym: str, label: str, comp: str, auth: str, ins: str = "?", altloc: str = ".", model: str = "1", auth_asym: str = "") -> tuple[str, ...]:
    return (asym, label, comp, altloc, auth_asym or asym, auth, ins, model)


class DrawTest(unittest.TestCase):
    def test_draw_is_deterministic_and_ignores_input_order_and_case(self):
        holdings = [f"{n}abc" for n in range(1, 10)] + ["1XYZ"]
        first = select_entries.draw_baseline(holdings, seed=7, size=4)
        second = select_entries.draw_baseline(list(reversed([h.upper() for h in holdings])), seed=7, size=4)
        self.assertEqual(first, second)
        self.assertEqual(first, sorted(first))
        self.assertEqual(len(first), 4)
        self.assertTrue(set(first) <= {h.lower() for h in holdings})

    def test_a_different_seed_draws_a_different_sample(self):
        holdings = [f"{n}abc" for n in range(1, 100)]
        self.assertNotEqual(
            select_entries.draw_baseline(holdings, seed=1, size=10),
            select_entries.draw_baseline(holdings, seed=2, size=10),
        )

    def test_the_fingerprint_names_the_snapshot_not_its_order(self):
        self.assertEqual(
            select_entries.fingerprint(["1ABC", "2xyz"]),
            select_entries.fingerprint(["2xyz", "1abc"]),
        )
        self.assertNotEqual(select_entries.fingerprint(["1abc"]), select_entries.fingerprint(["1abd"]))

    def test_pdbj_rows_become_sorted_lowercase_ids(self):
        response = {"results": [["1CBN"], ["1al4"], ["1cbn"]], "fields": ["pdbid"]}
        self.assertEqual(select_entries.pdbj_ids(response), ["1al4", "1cbn"])


class SeverityTest(unittest.TestCase):
    def test_a_plain_entry_has_no_classes(self):
        plain = block(atom_site(residue("A", "1", "GLY", "1"), residue("A", "2", "ALA", "2")))
        self.assertEqual(select_entries.severities(plain), {})

    def test_insertion_codes_count_residues(self):
        entry = block(atom_site(
            residue("A", "1", "GLY", "100", ins="A"),
            residue("A", "1", "GLY", "100", ins="A"),
            residue("A", "2", "ALA", "100", ins="B"),
            residue("A", "3", "SER", "101"),
        ))
        self.assertEqual(select_entries.severities(entry)["insertion_codes"], 2)

    def test_altlocs_count_residues(self):
        entry = block(atom_site(
            residue("A", "1", "GLY", "1", altloc="A"),
            residue("A", "1", "GLY", "1", altloc="B"),
            residue("A", "2", "ALA", "2", altloc="A"),
            residue("A", "3", "SER", "3"),
        ))
        self.assertEqual(select_entries.severities(entry)["altlocs"], 2)

    def test_microheterogeneity_takes_the_larger_of_observed_and_declared_positions(self):
        observed = block(atom_site(residue("A", "1", "SER", "1"), residue("A", "1", "PRO", "1")))
        self.assertEqual(select_entries.severities(observed)["microheterogeneity"], 1)
        declared = block(
            atom_site(residue("A", "1", "SER", "1")),
            "loop_\n_entity_poly_seq.entity_id\n_entity_poly_seq.num\n_entity_poly_seq.mon_id\n_entity_poly_seq.hetero\n"
            "1 1 SER y\n1 1 PRO y\n1 2 LEU y\n1 2 ILE y\n1 3 GLY n\n",
        )
        self.assertEqual(select_entries.severities(declared)["microheterogeneity"], 2)

    def test_numbering_divergence_counts_distinct_shifts_per_chain(self):
        entry = block(atom_site(
            residue("A", "1", "GLY", "5"),
            residue("A", "2", "ALA", "6"),
            residue("A", "3", "SER", "19"),
            residue("B", "1", "GLY", "1"),
            residue("C", "1", "GLY", "5"),
            residue("D", "1", "GLY", "X1"),
        ))
        self.assertEqual(select_entries.severities(entry)["numbering_divergence"], 4)

    def test_assembly_operators_take_the_largest_expansion(self):
        entry = block(
            atom_site(residue("A", "1", "GLY", "1")),
            "loop_\n_pdbx_struct_assembly_gen.assembly_id\n_pdbx_struct_assembly_gen.oper_expression\n"
            "_pdbx_struct_assembly_gen.asym_id_list\n1 '(1-60)(61-88)' A\n2 1-3 A\n3 1 A\n",
        )
        self.assertEqual(select_entries.severities(entry)["assembly_operators"], 1680)

    def test_branched_entries_count_branched_residues(self):
        entry = block(
            atom_site(residue("A", "1", "GLY", "1")),
            "_pdbx_entity_branch.entity_id 2\n_pdbx_entity_branch.type oligosaccharide\n"
            "loop_\n_pdbx_branch_scheme.asym_id\n_pdbx_branch_scheme.mon_id\nB NAG\nB NAG\nB BMA\n",
        )
        self.assertEqual(select_entries.severities(entry)["branched"], 3)

    def test_multichar_chain_ids_count_distinct_ids(self):
        entry = block(atom_site(
            residue("A", "1", "GLY", "1", auth_asym="A"),
            residue("AA", "1", "GLY", "1", auth_asym="B1"),
            residue("BA", "1", "GLY", "1", auth_asym="B1"),
        ))
        self.assertEqual(select_entries.severities(entry)["multichar_chain_ids"], 3)

    def test_multiple_models_count_models(self):
        entry = block(atom_site(
            residue("A", "1", "GLY", "1", model="1"),
            residue("A", "1", "GLY", "1", model="2"),
            residue("A", "1", "GLY", "1", model="3"),
        ))
        self.assertEqual(select_entries.severities(entry)["multiple_models"], 3)

    def test_every_class_has_a_severity(self):
        self.assertEqual(set(select_entries.SEVERITY), set(classify.PATHOLOGY_PREDICATES))


class ClassifyFileTest(unittest.TestCase):
    def test_unreadable_files_are_reported_by_id(self):
        with tempfile.TemporaryDirectory() as directory, mock.patch.object(select_entries, "DATA_DIR", Path(directory)):
            Path(directory, "1abc.cif").write_text("")
            Path(directory, "2abc.cif").write_text("data_x\n_cell.length_a")
            self.assertEqual(select_entries.classify_file("3abc"), ("3abc", "not downloaded"))
            for identifier in ("1abc", "2abc"):
                found, result = select_entries.classify_file(identifier)
                self.assertEqual(found, identifier)
                self.assertIsInstance(result, str)


class FillTest(unittest.TestCase):
    def test_seed_members_come_first_then_the_most_severe_candidates(self):
        seed = {"1abc": ["altlocs", "baseline"], "2abc": ["multiple_models"]}
        candidates = {
            "3abc": {"altlocs": 5},
            "4abc": {"altlocs": 9},
            "5abc": {"altlocs": 5},
            "6abc": {"altlocs": 1},
            "1abc": {"altlocs": 99},
        }
        filled = select_entries.fill_classes(seed, candidates, size=3)
        self.assertEqual(filled["altlocs"], ["1abc", "4abc", "3abc"])
        self.assertEqual(filled["multiple_models"], ["2abc"])
        self.assertEqual(filled["branched"], [])
        self.assertEqual(set(filled), set(classify.PATHOLOGY_PREDICATES))

    def test_a_class_with_more_seed_members_than_the_size_keeps_them_all(self):
        seed = {"1abc": ["altlocs"], "2abc": ["altlocs"]}
        filled = select_entries.fill_classes(seed, {"3abc": {"altlocs": 4}}, size=1)
        self.assertEqual(filled["altlocs"], ["1abc", "2abc"])


if __name__ == "__main__":
    unittest.main()
