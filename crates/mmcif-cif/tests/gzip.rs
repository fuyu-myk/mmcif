use std::path::PathBuf;

use mmcif_cif::{Document, Error};

const TEXT: &str = "data_x\n_cell.length_a 10.0\n";

fn scratch(name: &str, bytes: &[u8]) -> PathBuf {
    let path = std::env::temp_dir().join(format!("mmcif-cif-gzip-{}-{name}", std::process::id()));
    std::fs::write(&path, bytes).unwrap();
    path
}

#[test]
fn uncompressed_files_are_read_as_they_are() {
    let path = scratch("plain.cif", TEXT.as_bytes());
    assert_eq!(mmcif_cif::read_file(&path).unwrap(), TEXT.as_bytes());
    assert_eq!(
        Document::open(&path).unwrap().to_cif_string(),
        Document::parse(TEXT).unwrap().to_cif_string()
    );
}

#[cfg(feature = "gzip")]
mod with_gzip {
    use std::io::Write;

    use flate2::Compression;
    use flate2::write::GzEncoder;

    use super::*;

    fn compressed(bytes: &[u8]) -> Vec<u8> {
        let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
        encoder.write_all(bytes).unwrap();
        encoder.finish().unwrap()
    }

    #[test]
    fn compressed_text_and_binary_files_are_decompressed() {
        let path = scratch("text.cif.gz", &compressed(TEXT.as_bytes()));
        assert_eq!(mmcif_cif::read_file(&path).unwrap(), TEXT.as_bytes());
        let document = Document::open(&path).unwrap();
        assert_eq!(
            document.to_cif_string(),
            Document::parse(TEXT).unwrap().to_cif_string()
        );

        let binary = std::fs::read(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../fuzz/seeds/bcif/1crn.bcif"
        ))
        .unwrap();
        let path = scratch("crambin.bcif.gz", &compressed(&binary));
        let unpacked = Document::open_bcif(&path).unwrap();
        assert_eq!(
            unpacked.to_cif_string(),
            Document::from_bcif(&binary).unwrap().to_cif_string()
        );
    }

    #[test]
    fn a_corrupt_compressed_file_is_an_io_error() {
        let mut bytes = compressed(TEXT.as_bytes());
        let middle = bytes.len() / 2;
        bytes.truncate(middle);
        let path = scratch("truncated.cif.gz", &bytes);
        assert!(matches!(Document::open(&path), Err(Error::Io(_))));
    }
}

#[cfg(not(feature = "gzip"))]
#[test]
fn compressed_files_need_the_gzip_feature() {
    let path = scratch("text.cif.gz", &[0x1f, 0x8b, 0x08, 0x00]);
    assert!(matches!(
        mmcif_cif::read_file(&path),
        Err(Error::Compressed)
    ));
    assert!(matches!(Document::open(&path), Err(Error::Compressed)));
    assert!(matches!(Document::open_bcif(&path), Err(Error::Compressed)));
}
