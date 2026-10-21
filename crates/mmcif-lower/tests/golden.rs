use std::path::{Path, PathBuf};

use mmcif_lower::{Backend, lower};
use mmcif_spec::{Checked, Design, check_text};

fn fixtures() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures")
}

fn checked(yaml: &str) -> Checked {
    let checked = check_text("design.yaml", yaml, &fixtures());
    assert!(!checked.has_errors(), "{}", checked.render());
    checked
}

fn design(checked: &Checked) -> &Design {
    checked.design.as_ref().expect("the spec lowered")
}

fn emit(backend: Backend, yaml: &str) -> String {
    let checked = checked(yaml);
    match lower(backend, checked.source, design(&checked)) {
        Ok(text) => text,
        Err(problems) => {
            let mut report = mmcif_spec::Report::new();
            report.extend(problems);
            panic!("{}", checked.sources.render_all(&report));
        }
    }
}

fn refusal(backend: Backend, yaml: &str) -> String {
    let checked = checked(yaml);
    match lower(backend, checked.source, design(&checked)) {
        Ok(text) => panic!("expected a refusal, got:\n{text}"),
        Err(problems) => {
            let mut report = mmcif_spec::Report::new();
            report.extend(problems);
            checked.sources.render_all(&report)
        }
    }
}

const TARGET_AND_BINDER: &str = "entities:
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
  - protein:
      id: S
      sequence: 10..14C6C3
constraints:
  - bond:
      atom1: [S, 11, SG]
      atom2: [S, 18, SG]
";

#[test]
fn boltzgen_output_matches_its_golden_file() {
    let golden = std::fs::read_to_string(fixtures().join("target_and_binder.boltzgen.yaml"))
        .expect("the golden file is present");
    assert_eq!(emit(Backend::BoltzGen, TARGET_AND_BINDER), golden);
}

#[test]
fn rfd3_output_matches_its_golden_file() {
    let golden = std::fs::read_to_string(fixtures().join("target_and_motif.rfd3.json"))
        .expect("the golden file is present");
    let yaml = "entities:
  - file:
      path: target.cif
      include:
        - chain:
            id: A
      design:
        - chain:
            id: A
            res_index: 4..5
  - protein:
      id: S
      sequence: 10..14
";
    assert_eq!(emit(Backend::Rfd3, yaml), golden);
}

#[test]
fn a_de_novo_monomer_needs_no_input_structure() {
    let yaml = "entities:\n  - protein:\n      id: S\n      sequence: 50..80\n";
    assert_eq!(
        emit(Backend::Rfd3, yaml),
        "{\n  \"design\": {\n    \"contig\": \"50-80\"\n  }\n}\n"
    );
    assert_eq!(
        emit(Backend::BoltzGen, yaml),
        "entities:\n  - protein:\n      id: S\n      sequence: 50..80\n"
    );
}

#[test]
fn chain_breaks_separate_the_contig_pieces() {
    let yaml = "entities:
  - file:
      path: target.cif
      include:
        - chain:
            id: A
        - chain:
            id: B
  - protein:
      id: S
      sequence: 20
";
    let emitted = emit(Backend::Rfd3, yaml);
    assert!(
        emitted.contains("\"contig\": \"A1-6,/0,B1-2,/0,20\""),
        "{emitted}"
    );
}

#[test]
fn rfd3_refuses_a_fixed_residue_in_a_designed_chain() {
    let rendered = refusal(
        Backend::Rfd3,
        "entities:\n  - protein:\n      id: S\n      sequence: 10..14C6C3\n",
    );
    assert!(
        rendered.contains("rfd3 cannot fix the identity of a designed position"),
        "{rendered}"
    );
    assert!(rendered.contains("--> design.yaml:4:23"), "{rendered}");
}

#[test]
fn rfd3_refuses_bonds_and_secondary_structure() {
    let rendered = refusal(Backend::Rfd3, TARGET_AND_BINDER);
    assert!(
        rendered.contains("rfd3 cannot express secondary structure constraint"),
        "{rendered}"
    );
    assert!(
        rendered.contains("rfd3 cannot express a bond constraint"),
        "{rendered}"
    );
}

#[test]
fn rfd3_refuses_a_second_input_structure() {
    let rendered = refusal(
        Backend::Rfd3,
        "entities:
  - file:
      path: target.cif
      include:
        - chain:
            id: A
  - file:
      path: other.cif
      include:
        - chain:
            id: C
",
    );
    assert!(
        rendered.contains("rfd3 takes 1 input structure, not 2"),
        "{rendered}"
    );
    assert!(rendered.contains("merge the chains into one structure file"));
}

#[test]
fn rfd3_refuses_a_ligand_it_cannot_place() {
    let rendered = refusal(
        Backend::Rfd3,
        "entities:\n  - protein:\n      id: S\n      sequence: 20\n  - ligand:\n      id: Q\n      ccd: WHL\n",
    );
    assert!(rendered.contains("rfd3 cannot take a ligand"), "{rendered}");
}

#[test]
fn boltzgen_keeps_a_ligand_and_its_bond() {
    let yaml = "entities:
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
";
    assert_eq!(emit(Backend::BoltzGen, yaml), yaml);
}

#[test]
fn every_backend_is_reachable_by_name() {
    assert_eq!(Backend::named("boltzgen"), Some(Backend::BoltzGen));
    assert_eq!(Backend::named("rfd3"), Some(Backend::Rfd3));
    assert_eq!(Backend::named("nonsense"), None);
    assert_eq!(Backend::BoltzGen.name(), "boltzgen");
    assert_eq!(Backend::Rfd3.name(), "rfd3");
    assert_eq!(Backend::names(), ["boltzgen", "rfd3"]);
    assert!(!Backend::Rfd3.capabilities().bonds);
    assert!(Backend::BoltzGen.capabilities().bonds);
}

#[test]
fn every_backend_names_the_schema_it_was_written_against() {
    for backend in Backend::all() {
        let schema = backend.schema();
        assert!(!schema.is_empty(), "{} has no schema", backend.name());
        assert_ne!(schema, backend.name());
    }
    assert!(Backend::Rfd3.schema().contains("foundry"));
}
