#![no_main]

use libfuzzer_sys::fuzz_target;
use mmcif_cif::Document;

fuzz_target!(|data: &[u8]| {
    if let Ok(document) = Document::from_bytes(data.to_vec()) {
        let written = document.to_cif_string();
        let reparsed = Document::parse(written).expect("written documents parse");
        assert_eq!(reparsed, document);
    }
});
