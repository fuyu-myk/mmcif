use std::path::{Path, PathBuf};

use mmcif_ir::{Corner, GroupId, Level, SegmentId, SlotIndex};
use mmcif_spec::{Checked, check_path, check_text};

fn fixtures() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures")
}

fn check(yaml: &str) -> Checked {
    check_text("design.yaml", yaml, &fixtures())
}

fn target_entity() -> &'static str {
    concat!(
        "  - file:\n",
        "      path: target.cif\n",
        "      include:\n",
        "        - chain:\n",
        "            id: A\n",
    )
}

#[test]
fn a_clean_specification_lowers_without_complaint() {
    let checked = check(&format!(
        concat!(
            "entities:\n",
            "{}",
            "  - protein:\n",
            "      id: X\n",
            "      sequence: 10..14C6C3\n",
            "constraints:\n",
            "  - bond:\n",
            "      atom1: [X, 11, SG]\n",
            "      atom2: [X, 18, SG]\n",
        ),
        target_entity()
    ));
    assert!(checked.report.is_empty(), "{}", checked.render());
    let shortest = checked.shortest.expect("the spec instantiates");
    assert_eq!(
        shortest.segment_ids(),
        [SegmentId::new("A"), SegmentId::new("X")]
    );
    assert_eq!(shortest.slot_count(), 6 + 10 + 1 + 6 + 1 + 3);
    assert_eq!(shortest.constraints().len(), 1);
}

#[test]
fn file_slots_are_context_and_designed_slots_are_de_novo() {
    let checked = check(&format!(
        concat!(
            "entities:\n",
            "{}",
            "  - protein:\n",
            "      id: X\n",
            "      sequence: 8\n",
        ),
        target_entity()
    ));
    assert!(checked.report.is_empty(), "{}", checked.render());
    let shortest = checked.shortest.unwrap();
    let target = shortest.segment_slots(&SegmentId::new("A")).unwrap();
    for slot in target.iter() {
        assert_eq!(
            shortest.determinacy(slot).unwrap().corner(),
            Some(Corner::Context)
        );
    }
    let designed = shortest.segment_slots(&SegmentId::new("X")).unwrap();
    for slot in designed.iter() {
        assert_eq!(
            shortest.determinacy(slot).unwrap().corner(),
            Some(Corner::DeNovo)
        );
    }
}

#[test]
fn a_bonded_designed_cysteine_is_a_partial_state() {
    let checked = check(concat!(
        "entities:\n",
        "  - protein:\n",
        "      id: S\n",
        "      sequence: 10..14C6C3\n",
        "constraints:\n",
        "  - bond:\n",
        "      atom1: [S, 11, SG]\n",
        "      atom2: [S, 18, SG]\n",
    ));
    assert!(checked.report.is_empty(), "{}", checked.render());
    let shortest = checked.shortest.unwrap();
    let bonded = shortest.determinacy(SlotIndex::new(10)).unwrap();
    assert_eq!(
        (bonded.identity, bonded.geometry),
        (Level::Fixed, Level::Restricted)
    );
    assert!(bonded.is_partial());
    assert!(shortest.is_constrained(SlotIndex::new(10)).unwrap());
    assert!(!shortest.is_constrained(SlotIndex::new(9)).unwrap());
}

#[test]
fn a_disulfide_between_fixed_cysteines_is_accepted() {
    let checked = check(concat!(
        "entities:\n",
        "  - file:\n",
        "      path: target.cif\n",
        "      include:\n",
        "        - chain:\n",
        "            id: A\n",
        "constraints:\n",
        "  - bond:\n",
        "      atom1: [A, 2, SG]\n",
        "      atom2: [A, 3, SG]\n",
    ));
    assert!(checked.report.is_empty(), "{}", checked.render());
}

#[test]
fn a_bond_between_distant_fixed_atoms_is_refused() {
    let checked = check(concat!(
        "entities:\n",
        "  - file:\n",
        "      path: target.cif\n",
        "      include:\n",
        "        - chain:\n",
        "            id: A\n",
        "constraints:\n",
        "  - bond:\n",
        "      atom1: [A, 2, SG]\n",
        "      atom2: [A, 6, SG]\n",
    ));
    assert_eq!(checked.report.codes(), ["infeasible-geometry"]);
    assert_eq!(
        checked.render(),
        concat!(
            "error[infeasible-geometry]: A 2 SG and A 6 SG cannot satisfy every constraint\n",
            " --> design.yaml:8:5\n",
            "  |\n",
            "8 |   - bond:\n",
            "  |     ^^^^^ needs at least 15.03 angstrom, allows at most 2.50\n",
            "  = help: the structure file fixes them 15.03 angstrom apart\n",
        )
    );
}

#[test]
fn a_bond_on_a_missing_atom_lists_the_atoms_in_the_file() {
    let checked = check(concat!(
        "entities:\n",
        "  - file:\n",
        "      path: target.cif\n",
        "      include:\n",
        "        - chain:\n",
        "            id: A\n",
        "constraints:\n",
        "  - bond:\n",
        "      atom1: [A, 1, SG]\n",
        "      atom2: [A, 3, SG]\n",
    ));
    assert!(checked.report.codes().contains(&"unknown-atom"));
    let diagnostic = checked
        .report
        .diagnostics()
        .iter()
        .find(|entry| entry.code == "unknown-atom")
        .unwrap();
    assert_eq!(diagnostic.message, "GLY has no atom SG");
    assert_eq!(diagnostic.help.as_deref(), Some("atoms in the file: N, CA"));
}

#[test]
fn an_index_past_the_shortest_length_is_refused() {
    let checked = check(concat!(
        "entities:\n",
        "  - protein:\n",
        "      id: S\n",
        "      sequence: 10..14C6C3\n",
        "constraints:\n",
        "  - bond:\n",
        "      atom1: [S, 11, SG]\n",
        "      atom2: [S, 99, SG]\n",
    ));
    assert_eq!(checked.report.codes(), ["unknown-residue"]);
    assert_eq!(
        checked.report.diagnostics()[0].help.as_deref(),
        Some("chain S holds 21 residues at its shortest length")
    );
}

#[test]
fn a_bond_on_an_undeclared_chain_lists_the_declared_ones() {
    let checked = check(concat!(
        "entities:\n",
        "  - protein:\n",
        "      id: S\n",
        "      sequence: 17\n",
        "constraints:\n",
        "  - bond:\n",
        "      atom1: [S, 1, SG]\n",
        "      atom2: [Q, 1, CK]\n",
    ));
    assert_eq!(checked.report.codes(), ["unknown-chain"]);
    assert_eq!(
        checked.report.diagnostics()[0].help.as_deref(),
        Some("declared chains: S")
    );
}

#[test]
fn a_selection_that_matches_nothing_is_refused() {
    let checked = check(concat!(
        "entities:\n",
        "  - file:\n",
        "      path: target.cif\n",
        "      include:\n",
        "        - chain:\n",
        "            id: A\n",
        "            res_index: 40..50\n",
    ));
    assert_eq!(checked.report.codes(), ["empty-selection"]);
    assert_eq!(
        checked.report.diagnostics()[0].help.as_deref(),
        Some("chain A covers residues 1 to 6")
    );
}

#[test]
fn redesigning_residues_outside_the_selection_is_refused() {
    let checked = check(concat!(
        "entities:\n",
        "  - file:\n",
        "      path: target.cif\n",
        "      include:\n",
        "        - chain:\n",
        "            id: A\n",
        "            res_index: 1..3\n",
        "      design:\n",
        "        - chain:\n",
        "            id: A\n",
        "            res_index: 5..6\n",
    ));
    assert_eq!(checked.report.codes(), ["design-outside-selection"]);
}

#[test]
fn a_redesigned_file_residue_becomes_a_de_novo_slot() {
    let checked = check(concat!(
        "entities:\n",
        "  - file:\n",
        "      path: target.cif\n",
        "      include:\n",
        "        - chain:\n",
        "            id: A\n",
        "      design:\n",
        "        - chain:\n",
        "            id: A\n",
        "            res_index: 4..5\n",
    ));
    assert!(checked.report.is_empty(), "{}", checked.render());
    let shortest = checked.shortest.unwrap();
    let corners: Vec<Option<Corner>> = (0..6)
        .map(|slot| shortest.determinacy(SlotIndex::new(slot)).unwrap().corner())
        .collect();
    assert_eq!(
        corners,
        [
            Some(Corner::Context),
            Some(Corner::Context),
            Some(Corner::Context),
            Some(Corner::DeNovo),
            Some(Corner::DeNovo),
            Some(Corner::Context)
        ]
    );
}

#[test]
fn secondary_structure_groups_the_residues_it_names() {
    let checked = check(concat!(
        "entities:\n",
        "  - file:\n",
        "      path: target.cif\n",
        "      include:\n",
        "        - chain:\n",
        "            id: A\n",
        "      design:\n",
        "        - chain:\n",
        "            id: A\n",
        "            res_index: 4..6\n",
        "      secondary_structure:\n",
        "        - chain:\n",
        "            id: A\n",
        "            helix: 4..6\n",
    ));
    assert!(checked.report.is_empty(), "{}", checked.render());
    let shortest = checked.shortest.unwrap();
    let group = shortest.group(&GroupId::new("ss1:A")).unwrap();
    assert_eq!(
        group,
        [SlotIndex::new(3), SlotIndex::new(4), SlotIndex::new(5)]
    );
    assert_eq!(shortest.constraints().len(), 1);
}

#[test]
fn secondary_structure_on_file_residues_warns_without_failing() {
    let checked = check(concat!(
        "entities:\n",
        "  - file:\n",
        "      path: target.cif\n",
        "      include:\n",
        "        - chain:\n",
        "            id: A\n",
        "      secondary_structure:\n",
        "        - chain:\n",
        "            id: A\n",
        "            helix: 1..3\n",
    ));
    assert_eq!(checked.report.codes(), ["structure-on-fixed-residue"]);
    assert!(!checked.has_errors());
    assert_eq!(
        checked.report.diagnostics()[0].message,
        "3 residues already hold their structure"
    );
}

#[test]
fn a_file_without_an_include_takes_every_chain() {
    let checked = check(concat!(
        "entities:\n",
        "  - file:\n",
        "      path: target.cif\n",
    ));
    assert!(checked.report.is_empty(), "{}", checked.render());
    let shortest = checked.shortest.unwrap();
    assert_eq!(
        shortest.segment_ids(),
        [SegmentId::new("A"), SegmentId::new("B")]
    );
    assert_eq!(shortest.slot_count(), 8);
}

#[test]
fn an_absurd_length_is_refused_before_instantiation() {
    let checked = check(concat!(
        "entities:\n",
        "  - protein:\n",
        "      id: S\n",
        "      sequence: 20000000\n",
    ));
    assert_eq!(checked.report.codes(), ["unsatisfiable-lengths"]);
    assert!(checked.shortest.is_none());
}

#[test]
fn the_same_chain_twice_is_refused_with_both_places() {
    let checked = check(concat!(
        "entities:\n",
        "  - protein:\n",
        "      id: S\n",
        "      sequence: 8\n",
        "  - protein:\n",
        "      id: S\n",
        "      sequence: 9\n",
    ));
    assert_eq!(checked.report.codes(), ["duplicate-chain"]);
    assert_eq!(checked.report.diagnostics()[0].labels.len(), 2);
}

#[test]
fn a_missing_target_file_names_the_path_it_tried() {
    let checked = check(concat!(
        "entities:\n",
        "  - file:\n",
        "      path: nowhere.cif\n",
    ));
    assert_eq!(checked.report.codes(), ["unreadable-target"]);
    assert!(checked.render().contains("nowhere.cif"));
}

#[test]
fn a_rendered_report_points_at_the_offending_line() {
    let checked = check(concat!(
        "entities:\n",
        "  - protein:\n",
        "      id: S\n",
        "      sequence: 17\n",
        "constraints:\n",
        "  - bond:\n",
        "      atom1: [S, 1, SG]\n",
        "      atom2: [Q, 1, CK]\n",
    ));
    assert_eq!(
        checked.render(),
        concat!(
            "error[unknown-chain]: chain Q is not declared\n",
            " --> design.yaml:8:15\n",
            "  |\n",
            "8 |       atom2: [Q, 1, CK]\n",
            "  |               ^ no entity declares this chain\n",
            "  = help: declared chains: S\n",
        )
    );
}

#[test]
fn warnings_alone_leave_the_specification_usable() {
    let checked = check(concat!(
        "entities:\n",
        "  - file:\n",
        "      path: target.cif\n",
        "      binding_types:\n",
        "        - chain:\n",
        "            id: A\n",
        "            binding: 2..3\n",
    ));
    assert_eq!(checked.report.codes(), ["unsupported-key"]);
    assert!(!checked.has_errors());
    assert!(checked.shortest.is_some());
}

#[test]
fn a_specification_is_read_from_its_own_directory() {
    let checked = check_path(fixtures().join("peptide.yaml")).expect("the fixture is readable");
    assert!(checked.report.is_empty(), "{}", checked.render());
    let shortest = checked.shortest.unwrap();
    assert_eq!(
        shortest.segment_ids(),
        [SegmentId::new("A"), SegmentId::new("S")]
    );
    assert_eq!(shortest.constraints().len(), 1);
}

#[test]
fn an_unreadable_specification_is_an_error_not_a_report() {
    let error = check_path(fixtures().join("missing.yaml")).unwrap_err();
    assert!(error.to_string().contains("missing.yaml"));
}

#[test]
fn two_bonds_that_are_fine_alone_can_still_contradict_each_other() {
    let checked = check(concat!(
        "entities:\n",
        "  - file:\n",
        "      path: target.cif\n",
        "      include:\n",
        "        - chain:\n",
        "            id: A\n",
        "  - protein:\n",
        "      id: S\n",
        "      sequence: C4\n",
        "constraints:\n",
        "  - bond:\n",
        "      atom1: [S, 1, SG]\n",
        "      atom2: [A, 2, SG]\n",
        "  - bond:\n",
        "      atom1: [S, 1, SG]\n",
        "      atom2: [A, 6, SG]\n",
    ));
    assert_eq!(checked.report.codes(), ["infeasible-geometry"]);
    let rendered = checked.render();
    assert!(
        rendered.contains("needs at least 15.03 angstrom, allows at most 5.00"),
        "{rendered}"
    );
    assert!(rendered.contains("A 2 SG and A 6 SG"), "{rendered}");
    assert!(rendered.contains("--> design.yaml:11:5"), "{rendered}");
}

#[test]
fn a_length_range_whose_top_end_is_absurd_is_refused_cheaply() {
    let checked = check(concat!(
        "entities:\n",
        "  - protein:\n",
        "      id: S\n",
        "      sequence: 1..20000000\n",
    ));
    assert_eq!(checked.report.codes(), ["unsatisfiable-lengths"]);
    assert!(checked.shortest.is_none());
    assert!(
        checked
            .render()
            .contains("the longest sampled length is 20000000 slots")
    );
}

#[test]
fn a_chain_break_does_not_make_local_bonds_infeasible() {
    let checked = check(concat!(
        "entities:\n",
        "  - file:\n",
        "      path: gapped.cif\n",
        "      include:\n",
        "        - chain:\n",
        "            id: A\n",
        "      design:\n",
        "        - chain:\n",
        "            id: A\n",
        "            res_index: 2,50\n",
        "constraints:\n",
        "  - bond:\n",
        "      atom1: [A, 1, SG]\n",
        "      atom2: [A, 2, CA]\n",
        "  - bond:\n",
        "      atom1: [A, 51, SG]\n",
        "      atom2: [A, 50, CA]\n",
    ));
    assert!(checked.report.is_empty(), "{}", checked.render());
}

#[test]
fn atoms_from_different_files_are_not_measured_against_each_other() {
    let checked = check(concat!(
        "entities:\n",
        "  - file:\n",
        "      path: target.cif\n",
        "      include:\n",
        "        - chain:\n",
        "            id: A\n",
        "  - file:\n",
        "      path: other.cif\n",
        "      include:\n",
        "        - chain:\n",
        "            id: C\n",
        "constraints:\n",
        "  - bond:\n",
        "      atom1: [A, 6, SG]\n",
        "      atom2: [C, 1, SG]\n",
    ));
    assert!(checked.report.is_empty(), "{}", checked.render());
}

#[test]
fn an_atom_cannot_be_bonded_to_itself() {
    let checked = check(concat!(
        "entities:\n",
        "  - file:\n",
        "      path: target.cif\n",
        "      include:\n",
        "        - chain:\n",
        "            id: A\n",
        "constraints:\n",
        "  - bond:\n",
        "      atom1: [A, 2, SG]\n",
        "      atom2: [A, 2, SG]\n",
    ));
    assert_eq!(checked.report.codes(), ["self-reference"]);
    assert_eq!(
        checked.report.diagnostics()[0].message,
        "A 2 SG cannot be constrained against itself"
    );
}

#[test]
fn redesigning_a_chain_that_is_not_included_is_refused() {
    let checked = check(concat!(
        "entities:\n",
        "  - file:\n",
        "      path: target.cif\n",
        "      include:\n",
        "        - chain:\n",
        "            id: A\n",
        "      design:\n",
        "        - chain:\n",
        "            id: Z\n",
        "            res_index: 1..2\n",
    ));
    assert_eq!(checked.report.codes(), ["unknown-chain"]);
    let rendered = checked.render();
    assert!(rendered.contains("--> design.yaml:9:17"), "{rendered}");
    assert!(rendered.contains("9 |             id: Z"), "{rendered}");
}

#[test]
fn an_implicit_chain_points_at_the_file_it_came_from() {
    let checked = check(concat!(
        "entities:\n",
        "  - protein:\n",
        "      id: A\n",
        "      sequence: 8\n",
        "  - file:\n",
        "      path: target.cif\n",
    ));
    assert_eq!(checked.report.codes(), ["duplicate-chain"]);
    let rendered = checked.render();
    assert!(rendered.contains("--> design.yaml:6:13"), "{rendered}");
}

#[test]
fn lengths_that_only_overflow_together_are_refused() {
    let checked = check(concat!(
        "entities:\n",
        "  - protein:\n",
        "      id: S\n",
        "      sequence: 10000000\n",
        "  - protein:\n",
        "      id: T\n",
        "      sequence: 10000000\n",
    ));
    assert_eq!(checked.report.codes(), ["unsatisfiable-lengths"]);
    assert!(
        checked
            .render()
            .contains("the longest sampled length is 20000000 slots")
    );
}

#[test]
fn a_specification_with_too_many_constrained_atoms_says_so() {
    let mut yaml = String::from("entities:\n  - protein:\n      id: S\n      sequence: ");
    yaml.push_str(&"C".repeat(400));
    yaml.push_str("\nconstraints:\n");
    for position in (1..400).step_by(2) {
        yaml.push_str(&format!(
            "  - bond:\n      atom1: [S, {position}, SG]\n      atom2: [S, {}, SG]\n",
            position + 1
        ));
    }
    let checked = check(&yaml);
    assert_eq!(checked.report.codes(), ["geometry-not-checked"]);
    assert!(!checked.has_errors());
}
