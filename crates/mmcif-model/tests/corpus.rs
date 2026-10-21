use std::path::PathBuf;

use mmcif_model::{Model, read_bcif, read_mmcif};

const TOLERANCE: f64 = 1e-3;

fn corpus_dir() -> PathBuf {
    std::env::var_os("MMCIF_CORPUS_DIR")
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from(concat!(env!("CARGO_MANIFEST_DIR"), "/../../corpus/data")))
}

type AtomKey = (String, i32, String, String, String, String);

fn atoms(model: &Model) -> Vec<(AtomKey, [f64; 3])> {
    model
        .atoms()
        .map(|atom| {
            let residue = atom.residue();
            let chain = residue.chain();
            let key = (
                chain.label_asym().to_string(),
                residue.auth_seq().get(),
                residue.ins_code().to_string(),
                residue.comp().to_string(),
                atom.name().to_string(),
                atom.altloc().to_string(),
            );
            (key, atom.position())
        })
        .collect()
}

#[test]
#[ignore = "needs the downloaded corpus"]
fn binary_and_text_entries_give_the_same_models() {
    let dir = corpus_dir();
    let mut compared = 0;
    for entry in
        std::fs::read_dir(&dir).unwrap_or_else(|error| panic!("{}: {error}", dir.display()))
    {
        let bcif = entry.unwrap().path();
        if bcif.extension().is_none_or(|e| e != "bcif") {
            continue;
        }
        let cif = bcif.with_extension("cif");
        let text = read_mmcif(&cif).unwrap_or_else(|error| panic!("{}: {error}", cif.display()));
        let binary = read_bcif(&bcif).unwrap_or_else(|error| panic!("{}: {error}", bcif.display()));
        assert_eq!(
            text.models().len(),
            binary.models().len(),
            "{}",
            bcif.display()
        );
        for (text_model, binary_model) in text.models().iter().zip(binary.models()) {
            let expected = atoms(text_model);
            let actual = atoms(binary_model);
            assert_eq!(expected.len(), actual.len(), "{}", bcif.display());
            for ((expected_key, expected_at), (actual_key, actual_at)) in
                expected.iter().zip(&actual)
            {
                assert_eq!(expected_key, actual_key, "{}", bcif.display());
                let close = expected_at
                    .iter()
                    .zip(actual_at)
                    .all(|(a, b)| (a - b).abs() <= TOLERANCE);
                assert!(
                    close,
                    "{} {expected_key:?}: {expected_at:?} vs {actual_at:?}",
                    bcif.display()
                );
            }
        }
        assert_eq!(
            text.entities().len(),
            binary.entities().len(),
            "{}",
            bcif.display()
        );
        assert_eq!(text.assemblies(), binary.assemblies(), "{}", bcif.display());
        compared += 1;
    }
    assert!(compared > 0, "no .bcif files in {}", dir.display());
}

#[test]
#[ignore = "needs the downloaded corpus"]
fn typed_binary_columns_match_the_text_path_exactly() {
    let dir = corpus_dir();
    let mut compared = 0;
    for entry in std::fs::read_dir(&dir).unwrap() {
        let path = entry.unwrap().path();
        if path.extension().is_none_or(|e| e != "bcif") {
            continue;
        }
        let bytes = std::fs::read(&path).unwrap();
        let typed = mmcif_model::parse_bcif(&bytes).unwrap();
        let text = mmcif_cif::Document::from_bcif(&bytes)
            .unwrap()
            .to_cif_string();
        let via_text = mmcif_model::parse_mmcif(&text).unwrap();
        assert!(
            typed == via_text,
            "{} differs between typed and text paths",
            path.display()
        );
        compared += 1;
    }
    assert!(compared > 0);
}
