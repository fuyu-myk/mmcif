use std::path::{Path, PathBuf};

use mmcif_ir::{Corner, SlotIndex, Spec};
use mmcif_lower::{Backend, lower};
use mmcif_spec::{Checked, check_text};

fn fixtures() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures")
}

fn checked(yaml: &str) -> Checked {
    let checked = check_text("design.yaml", yaml, &fixtures());
    assert!(!checked.has_errors(), "{yaml}\n{}", checked.render());
    checked
}

fn emitted(yaml: &str) -> String {
    let checked = checked(yaml);
    let design = checked.design.as_ref().expect("the spec lowered");
    match lower(Backend::BoltzGen, checked.source, design) {
        Ok(text) => text,
        Err(problems) => {
            let mut report = mmcif_spec::Report::new();
            report.extend(problems);
            panic!("{}", checked.sources.render_all(&report));
        }
    }
}

fn spec_of(yaml: &str) -> Spec {
    checked(yaml).design.expect("the spec lowered").spec
}

fn corners(yaml: &str) -> Vec<Option<Corner>> {
    let checked = checked(yaml);
    let shortest = checked.shortest.expect("the spec instantiates");
    (0..shortest.slot_count() as u32)
        .map(|slot| {
            shortest
                .determinacy(SlotIndex::new(slot))
                .expect("every slot has a determinacy")
                .corner()
        })
        .collect()
}

fn round_trip(yaml: &str) {
    let once = emitted(yaml);
    assert_eq!(spec_of(yaml), spec_of(&once), "specs differ\n{once}");
    assert_eq!(corners(yaml), corners(&once), "determinacy differs\n{once}");
    assert_eq!(once, emitted(&once), "a second pass changed the output");
}

#[test]
fn a_target_with_a_disulfide_binder_survives_a_round_trip() {
    round_trip(
        "entities:
  - file:
      path: target.cif
      include:
        - chain:
            id: A
            res_index: 1..4
  - protein:
      id: S
      sequence: 10..14C6C3
constraints:
  - bond:
      atom1: [S, 11, SG]
      atom2: [S, 18, SG]
",
    );
}

#[test]
fn a_redesigned_region_and_its_secondary_structure_survive() {
    round_trip(
        "entities:
  - file:
      path: target.cif
      include:
        - chain:
            id: A
      design:
        - chain:
            id: A
            res_index: 4..5
      secondary_structure:
        - chain:
            id: A
            helix: 4..5
",
    );
}

#[test]
fn several_chains_and_a_ligand_survive() {
    round_trip(
        "entities:
  - file:
      path: target.cif
      include:
        - chain:
            id: A
        - chain:
            id: B
  - protein:
      id: X
      sequence: AAVTTTTPPP
  - protein:
      id: G
      sequence: 15..20AAAAAAVTTTT18PPP
  - ligand:
      id: Q
      ccd: WHL
",
    );
}

#[test]
fn a_bare_monomer_survives() {
    round_trip("entities:\n  - protein:\n      id: S\n      sequence: 50..80\n");
}

#[test]
fn a_whole_file_include_becomes_explicit_without_changing_meaning() {
    let yaml = "entities:\n  - file:\n      path: target.cif\n";
    let once = emitted(yaml);
    assert!(once.contains("res_index: 1..6"), "{once}");
    assert!(once.contains("res_index: 1..2"), "{once}");
    assert_eq!(spec_of(yaml), spec_of(&once));
}

#[test]
fn an_emitted_specification_is_accepted_by_the_checker() {
    let once = emitted(
        "entities:
  - file:
      path: target.cif
      include:
        - chain:
            id: A
      design:
        - chain:
            id: A
            res_index: 4..6
  - protein:
      id: S
      sequence: 8C4
constraints:
  - bond:
      atom1: [S, 9, SG]
      atom2: [A, 2, SG]
",
    );
    let again = check_text("emitted.yaml", &once, &fixtures());
    assert!(!again.has_errors(), "{once}\n{}", again.render());
    assert!(again.report.is_empty(), "{}", again.render());
}

#[test]
fn a_specification_using_every_supported_feature_survives() {
    round_trip(
        "entities:
  - file:
      path: target.cif
      include:
        - chain:
            id: A
            res_index: 1..4
        - chain:
            id: B
      design:
        - chain:
            id: A
            res_index: 4
      secondary_structure:
        - chain:
            id: A
            helix: 4
  - protein:
      id: X
      sequence: AAVTTTTPPP
  - protein:
      id: G
      sequence: 15..20AAAAAAVTTTT18PPP
  - protein:
      id: R
      sequence: 3..5C6C3
  - ligand:
      id: Q
      ccd: WHL
constraints:
  - bond:
      atom1: [R, 4, SG]
      atom2: [Q, 1, CK]
  - bond:
      atom1: [R, 11, SG]
      atom2: [Q, 1, CH]
",
    );
}
