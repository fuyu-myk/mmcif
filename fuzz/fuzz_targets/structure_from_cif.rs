#![no_main]

use libfuzzer_sys::fuzz_target;
use mmcif_model::{AltlocPolicy, parse_mmcif};

fuzz_target!(|data: &[u8]| {
    let Ok(text) = std::str::from_utf8(data) else {
        return;
    };
    if let Ok(structure) = parse_mmcif(text) {
        if let Some(model) = structure.first_model() {
            model.single_conformer(AltlocPolicy::HighestOccupancy);
            for assembly in structure.assemblies() {
                let _ = structure.assemble(assembly.id(), model);
            }
        }
    }
});
