use std::fs::File;
use std::io::BufWriter;
use std::path::{Path, PathBuf};
use std::process::Command;

use mmcif_ccd::{AtomStereo, BondOrder, CcdIndex, Component, build_index};
use mmcif_cif::Document;

const RELEASE: &str = "2026-10-20";
const SORTED_IDS: [&str; 8] = ["A1AAA", "ALA", "ATP", "HEM", "HOH", "MSE", "NAG", "ZN"];

fn fixture_path() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/components-sample.cif")
}

fn sample_index() -> CcdIndex {
    let text = std::fs::read_to_string(fixture_path()).unwrap();
    let mut out = Vec::new();
    let summary = build_index(&text, RELEASE, &mut out).unwrap();
    assert_eq!(summary.components(), 8);
    assert_eq!((summary.atoms(), summary.bonds()), (232, 240));
    assert_eq!(summary.bytes(), out.len() as u64);
    CcdIndex::from_bytes(out).unwrap()
}

fn component(index: &CcdIndex, id: &str) -> Component {
    index
        .get(id)
        .unwrap()
        .unwrap_or_else(|| panic!("{id} is not indexed"))
}

fn count_orders(component: &Component, order: BondOrder) -> usize {
    component
        .bonds()
        .iter()
        .filter(|bond| bond.order() == order)
        .count()
}

#[test]
fn indexes_every_component_in_id_order() {
    let index = sample_index();
    assert_eq!(index.release(), RELEASE);
    assert_eq!(index.len(), 8);
    assert_eq!(index.ids().collect::<Vec<_>>(), SORTED_IDS);
    assert_eq!(index.get("XXX").unwrap(), None);
    assert_eq!(index.get("ala").unwrap(), None);
}

#[test]
fn lookups_match_components_read_directly() {
    let index = sample_index();
    let document = Document::open(fixture_path()).unwrap();
    for block in document.blocks() {
        let expected = Component::from_block(&block).unwrap();
        assert_eq!(index.get(block.name()).unwrap(), Some(expected));
    }
}

#[test]
fn alanine() {
    let ala = component(&sample_index(), "ALA");
    assert_eq!(ala.name(), "ALANINE");
    assert_eq!(ala.kind(), "L-PEPTIDE LINKING");
    assert_eq!(ala.one_letter_code(), Some("A"));
    assert_eq!(ala.parent(), None);
    assert_eq!(ala.formula(), Some("C3 H7 N O2"));
    assert_eq!((ala.formal_charge(), ala.release_status()), (0, "REL"));
    assert_eq!((ala.atoms().len(), ala.bonds().len()), (13, 12));
    let atom = |name: &str| ala.atoms().iter().find(|atom| atom.name() == name).unwrap();
    assert_eq!(atom("CA").stereo(), AtomStereo::S);
    assert!(atom("OXT").is_leaving() && !atom("CA").is_leaving());
    let first = ala.bonds()[0];
    let names = first
        .atom_indices()
        .map(|index| ala.atoms()[usize::from(index)].name());
    assert_eq!(names, ["N", "CA"]);
}

#[test]
fn selenomethionine_names_its_parent() {
    let mse = component(&sample_index(), "MSE");
    assert_eq!(mse.parent(), Some("MET"));
    assert_eq!(mse.one_letter_code(), Some("M"));
    assert!(mse.atoms().iter().any(|atom| atom.element() == "SE"));
}

#[test]
fn heme_coordinates_iron() {
    let hem = component(&sample_index(), "HEM");
    assert_eq!((hem.atoms().len(), hem.bonds().len()), (75, 82));
    let (iron_index, iron) = hem
        .atoms()
        .iter()
        .enumerate()
        .find(|(_, atom)| atom.name() == "FE")
        .unwrap();
    assert_eq!(iron.element(), "FE");
    assert_eq!(iron.charge(), 0);
    let iron_bonds = hem
        .bonds()
        .iter()
        .filter(|bond| bond.atom_indices().contains(&(iron_index as u16)))
        .count();
    assert_eq!(iron_bonds, 4);
    assert_eq!(count_orders(&hem, BondOrder::Double), 15);
    assert_eq!(
        hem.bonds().iter().filter(|bond| bond.is_aromatic()).count(),
        10
    );
}

#[test]
fn n_acetylglucosamine_bond_orders() {
    let nag = component(&sample_index(), "NAG");
    assert_eq!(nag.kind(), "D-saccharide, beta linking");
    assert_eq!(nag.bonds().len(), 30);
    assert_eq!(count_orders(&nag, BondOrder::Single), 29);
    let double: Vec<_> = nag
        .bonds()
        .iter()
        .filter(|bond| bond.order() == BondOrder::Double)
        .map(|bond| {
            bond.atom_indices()
                .map(|i| nag.atoms()[usize::from(i)].name())
        })
        .collect();
    assert_eq!(double, [["C7", "O7"]]);
}

#[test]
fn small_components() {
    let index = sample_index();
    let zinc = component(&index, "ZN");
    assert_eq!((zinc.formal_charge(), zinc.kind()), (2, "NON-POLYMER"));
    assert_eq!((zinc.atoms().len(), zinc.bonds().len()), (1, 0));
    assert_eq!(zinc.atoms()[0].charge(), 2);
    let water = component(&index, "HOH");
    assert_eq!((water.atoms().len(), water.bonds().len()), (3, 2));
    let atp = component(&index, "ATP");
    assert_eq!((atp.atoms().len(), atp.bonds().len()), (47, 49));
    assert_eq!(count_orders(&atp, BondOrder::Double), 7);
}

#[test]
fn five_character_ids_are_supported() {
    let ligand = component(&sample_index(), "A1AAA");
    assert_eq!(ligand.kind(), "non-polymer");
    assert_eq!(ligand.formula(), Some("C17 H18 N8"));
    assert_eq!((ligand.atoms().len(), ligand.bonds().len()), (43, 46));
}

#[test]
fn opens_index_files_from_disk() {
    let path = Path::new(env!("CARGO_TARGET_TMPDIR")).join("mmcif-ccd-open.idx");
    let text = std::fs::read_to_string(fixture_path()).unwrap();
    build_index(&text, RELEASE, BufWriter::new(File::create(&path).unwrap())).unwrap();
    let mapped = CcdIndex::open(&path).unwrap();
    let owned = sample_index();
    assert_eq!(mapped.release(), RELEASE);
    for id in SORTED_IDS {
        assert_eq!(mapped.get(id).unwrap(), owned.get(id).unwrap());
    }
}

#[test]
fn build_binary_writes_an_index() {
    let path = Path::new(env!("CARGO_TARGET_TMPDIR")).join("mmcif-ccd-bin.idx");
    let output = Command::new(env!("CARGO_BIN_EXE_mmcif-ccd-build"))
        .arg(fixture_path())
        .arg(&path)
        .args(["--release", RELEASE])
        .output()
        .unwrap();
    assert!(output.status.success(), "{output:?}");
    let stdout = String::from_utf8(output.stdout).unwrap();
    assert!(
        stdout.starts_with("indexed 8 components (232 atoms, 240 bonds)"),
        "{stdout}"
    );
    let index = CcdIndex::open(&path).unwrap();
    assert_eq!(index.ids().collect::<Vec<_>>(), SORTED_IDS);
}

#[test]
fn build_binary_requires_a_release() {
    let output = Command::new(env!("CARGO_BIN_EXE_mmcif-ccd-build"))
        .arg(fixture_path())
        .arg("unused.idx")
        .output()
        .unwrap();
    assert_eq!(output.status.code(), Some(2));
    assert!(String::from_utf8_lossy(&output.stderr).starts_with("usage:"));
}
