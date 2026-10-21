use mmcif_ccd::{BondOrder, CcdIndex, build_index_from_file};

#[test]
#[ignore = "needs MMCIF_CCD_FILE pointing at components.cif"]
fn indexes_every_component_of_the_full_dictionary() {
    let path = std::env::var("MMCIF_CCD_FILE").expect("MMCIF_CCD_FILE is set");
    let mut bytes = Vec::new();
    let summary = build_index_from_file(&path, "test", &mut bytes).unwrap();
    let text = std::fs::read_to_string(&path).unwrap();
    let blocks = text
        .lines()
        .filter(|line| line.starts_with("data_"))
        .count();
    assert_eq!(summary.components(), blocks);

    let index = CcdIndex::from_bytes(bytes).unwrap();
    assert_eq!(index.len(), blocks);
    let alanine = index.get("ALA").unwrap().unwrap();
    assert_eq!(alanine.one_letter_code(), Some("A"));
    let chromophore = index.get("0YG").unwrap().unwrap();
    assert_eq!(chromophore.one_letter_code(), Some("YG"));
    let heme = index.get("HEM").unwrap().unwrap();
    assert!(
        heme.bonds()
            .iter()
            .any(|bond| bond.order() == BondOrder::Double)
    );
    for id in index.ids().step_by(97) {
        let component = index.get(&id).unwrap().unwrap();
        assert_eq!(component.id(), id);
    }
}
