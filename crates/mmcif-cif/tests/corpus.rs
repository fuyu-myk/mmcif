use std::path::PathBuf;

use mmcif_cif::Document;

fn corpus_files(extension: &str) -> Vec<PathBuf> {
    let dir = std::env::var_os("MMCIF_CORPUS_DIR")
        .map(PathBuf::from)
        .unwrap_or_else(|| {
            PathBuf::from(concat!(env!("CARGO_MANIFEST_DIR"), "/../../corpus/data"))
        });
    let mut files: Vec<PathBuf> = std::fs::read_dir(&dir)
        .unwrap_or_else(|error| panic!("corpus directory {}: {error}", dir.display()))
        .map(|entry| entry.unwrap().path())
        .filter(|path| path.extension().is_some_and(|e| e == extension))
        .collect();
    files.sort();
    assert!(
        !files.is_empty(),
        "no .{extension} files in {}",
        dir.display()
    );
    files
}

#[test]
#[ignore = "needs the downloaded corpus"]
fn every_text_entry_round_trips() {
    for path in corpus_files("cif") {
        let document =
            Document::open(&path).unwrap_or_else(|error| panic!("{}: {error}", path.display()));
        let again = Document::parse(document.to_cif_string()).unwrap();
        assert!(
            again == document,
            "{} changed after write then parse",
            path.display()
        );
    }
}

#[test]
#[ignore = "needs the downloaded corpus"]
fn every_binary_entry_decodes_and_round_trips() {
    for path in corpus_files("bcif") {
        let document = Document::open_bcif(&path)
            .unwrap_or_else(|error| panic!("{}: {error}", path.display()));
        let again = Document::parse(document.to_cif_string()).unwrap();
        assert!(
            again == document,
            "{} changed after write then parse",
            path.display()
        );
    }
}
