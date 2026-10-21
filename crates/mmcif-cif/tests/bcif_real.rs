use std::env;

use mmcif_cif::{Category, Document};

fn open_from_env() -> (Document, Document) {
    let path = |name: &str| env::var(name).unwrap_or_else(|_| panic!("{name} is not set"));
    let binary = Document::open_bcif(path("MMCIF_BCIF_FILE")).unwrap();
    let text = Document::open(path("MMCIF_CIF_FILE")).unwrap();
    (binary, text)
}

fn column_text<'a>(category: &Category<'a>, item: &str) -> Vec<&'a str> {
    let column = category
        .column(item)
        .unwrap_or_else(|| panic!("{item} is missing"));
    column.iter().map(|value| value.raw()).collect()
}

#[test]
#[ignore = "needs MMCIF_BCIF_FILE and MMCIF_CIF_FILE"]
fn binary_and_text_atom_sites_agree() {
    let (binary, text) = open_from_env();
    let binary_block = binary.block(0).unwrap();
    let text_block = text.block(0).unwrap();
    let binary_atoms = binary_block.category("atom_site").unwrap();
    let text_atoms = text_block.category("atom_site").unwrap();
    assert_eq!(binary_atoms.row_count(), text_atoms.row_count());

    for item in ["Cartn_x", "Cartn_y", "Cartn_z"] {
        let binary_values = column_text(&binary_atoms, item);
        let text_values = column_text(&text_atoms, item);
        for (row, (b, t)) in binary_values.iter().zip(&text_values).enumerate() {
            let b: f64 = b.parse().unwrap();
            let t: f64 = t.parse().unwrap();
            assert!((b - t).abs() < 1e-3, "{item} row {row}: {b} vs {t}");
        }
    }
    for item in ["label_atom_id", "label_comp_id", "auth_seq_id"] {
        assert_eq!(
            column_text(&binary_atoms, item),
            column_text(&text_atoms, item),
            "{item}"
        );
    }
}
