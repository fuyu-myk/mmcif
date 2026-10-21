#![no_main]

use libfuzzer_sys::fuzz_target;
use mmcif_cif::{BcifLimits, Document};

const LIMITS: BcifLimits = BcifLimits {
    max_elements: 1 << 20,
    max_text_bytes: 1 << 24,
};

fuzz_target!(|data: &[u8]| {
    if let Ok(document) = Document::from_bcif_with_limits(data, LIMITS) {
        let reparsed = Document::parse(document.to_cif_string()).expect("decoded documents write valid CIF");
        assert_eq!(reparsed, document);
        let typed = mmcif_model::Structure::from_document(&document);
        let via_text = mmcif_model::Structure::from_document(&reparsed);
        match (typed, via_text) {
            (Ok(typed), Ok(via_text)) => assert_eq!(format!("{typed:?}"), format!("{via_text:?}")),
            (Err(_), Err(_)) => {}
            (typed, via_text) => panic!("typed {:?} but text {:?}", typed.err(), via_text.err()),
        }
    }
});
