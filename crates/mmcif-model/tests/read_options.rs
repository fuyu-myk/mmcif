use mmcif_model::{
    AssemblyError, AssemblyId, AuthResidueKey, AuthSeqId, Format, InsCode, LabelSeqId, ReadError,
    ReadOptions, Structure, read_with,
};

const FIXTURE: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/tests/fixtures/pathologies.cif"
);
const LEGACY: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/fixtures/legacy.pdb");
const CRAMBIN_BCIF: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../fuzz/seeds/bcif/1crn.bcif"
);

const BROKEN_ENTITIES: (&str, &str) = ("1 2 SER y\n", "1 x SER y\n");
const BROKEN_SCHEME: (&str, &str) = (
    "AA 1 4 LYS 101 101 LYS H1 . n\n",
    "AA 1 4 LYS x 101 LYS H1 . n\n",
);
const BROKEN_OPERATORS: (&str, &str) = (
    "1 'identity operation' 1 0 0 0 0 1 0 0 0 0 1 0\n",
    "1 'identity operation' 1 0 0 0 0 1 0 0 0 0 1 ?\n",
);
const BROKEN_ASSEMBLIES: (&str, &str) = ("1 1 AA,B,C,D\n", "1 ? AA,B,C,D\n");

fn broken(breaks: &[(&str, &str)]) -> String {
    let mut text = std::fs::read_to_string(FIXTURE).unwrap();
    for (from, to) in breaks {
        assert!(text.contains(from), "the fixture lacks {from:?}");
        text = text.replace(from, to);
    }
    text
}

fn lenient(text: &str) -> mmcif_model::Loaded {
    let document = mmcif_cif::Document::parse(text).unwrap();
    Structure::from_block_with(
        &document.block(0).unwrap(),
        &ReadOptions::new().lenient(true),
    )
    .unwrap()
}

#[test]
fn reads_only_the_first_model_when_asked() {
    let options = ReadOptions::new().first_model_only(true);
    let loaded = read_with(FIXTURE, &options).unwrap();
    let numbers: Vec<i32> = loaded
        .structure
        .models()
        .iter()
        .map(|m| m.number().get())
        .collect();
    assert_eq!(numbers, [1]);
    assert_eq!(loaded.structure.first_model().unwrap().atom_count(), 16);

    let pdb = read_with(LEGACY, &options).unwrap();
    assert_eq!(pdb.structure.models().len(), 1);
    assert_eq!(mmcif_model::read(LEGACY).unwrap().models().len(), 2);
}

#[test]
fn a_strict_read_fails_on_any_malformed_optional_category() {
    for breaks in [
        BROKEN_ENTITIES,
        BROKEN_SCHEME,
        BROKEN_OPERATORS,
        BROKEN_ASSEMBLIES,
    ] {
        let text = broken(&[breaks]);
        assert!(mmcif_model::parse_mmcif(&text).is_err(), "{breaks:?}");
        let document = mmcif_cif::Document::parse(&text).unwrap();
        assert!(
            Structure::from_block_with(&document.block(0).unwrap(), &ReadOptions::new()).is_err()
        );
    }
}

#[test]
fn a_lenient_read_skips_malformed_entities_and_reports_them() {
    let loaded = lenient(&broken(&[BROKEN_ENTITIES]));
    assert!(loaded.structure.entities().is_empty());
    assert_eq!(loaded.structure.first_model().unwrap().atom_count(), 16);
    assert_eq!(loaded.structure.assemblies().len(), 3);
    assert!(matches!(
        loaded.skipped.as_slice(),
        [ReadError::InvalidValue {
            category: "entity_poly_seq",
            item: "num",
            ..
        }]
    ));
}

#[test]
fn a_lenient_read_numbers_residues_from_atoms_when_a_scheme_is_malformed() {
    let loaded = lenient(&broken(&[BROKEN_SCHEME]));
    let auth = AuthResidueKey {
        chain: "H1".into(),
        seq: AuthSeqId::new(101),
        ins_code: InsCode::NONE,
    };
    let label = loaded.structure.numbering().to_label(&auth).unwrap();
    assert_eq!(label.seq, LabelSeqId::new(4).unwrap());
    assert!(matches!(
        loaded.skipped.as_slice(),
        [ReadError::InvalidValue {
            category: "pdbx_poly_seq_scheme",
            ..
        }]
    ));
}

#[test]
fn a_lenient_read_skips_malformed_operators_and_assemblies() {
    let loaded = lenient(&broken(&[BROKEN_OPERATORS]));
    assert!(loaded.structure.operators().is_empty());
    assert_eq!(loaded.structure.assemblies().len(), 3);
    let model = loaded.structure.first_model().unwrap();
    assert!(matches!(
        loaded.structure.assemble(&AssemblyId::new("1"), model),
        Err(AssemblyError::UnknownOperator { .. })
    ));
    assert_eq!(loaded.skipped.len(), 1);

    let loaded = lenient(&broken(&[BROKEN_ASSEMBLIES]));
    assert!(loaded.structure.assemblies().is_empty());
    assert_eq!(loaded.structure.operators().len(), 3);
    assert_eq!(loaded.skipped.len(), 1);
}

#[test]
fn a_lenient_read_reports_every_skipped_group_in_order() {
    let loaded = lenient(&broken(&[
        BROKEN_ENTITIES,
        BROKEN_SCHEME,
        BROKEN_OPERATORS,
        BROKEN_ASSEMBLIES,
    ]));
    let categories: Vec<&str> = loaded
        .skipped
        .iter()
        .map(|error| match error {
            ReadError::InvalidValue { category, .. } => *category,
            other => panic!("unexpected {other}"),
        })
        .collect();
    assert_eq!(
        categories,
        [
            "entity_poly_seq",
            "pdbx_poly_seq_scheme",
            "pdbx_struct_oper_list",
            "pdbx_struct_assembly_gen"
        ]
    );
    assert!(lenient(&broken(&[])).skipped.is_empty());
}

#[test]
fn read_with_keeps_the_document_of_cif_and_binary_files() {
    let options = ReadOptions::new();
    let cif = read_with(FIXTURE, &options).unwrap();
    let document = cif.document.unwrap();
    assert_eq!(document.block(0).unwrap().name(), "TEST");
    assert!(cif.skipped.is_empty());

    let bcif = read_with(CRAMBIN_BCIF, &options).unwrap();
    assert!(
        bcif.document
            .unwrap()
            .block(0)
            .unwrap()
            .category("atom_site")
            .is_some()
    );
    assert_eq!(bcif.structure.first_model().unwrap().atom_count(), 327);

    assert!(read_with(LEGACY, &options).unwrap().document.is_none());
    assert!(matches!(
        read_with("structure.xyz", &options),
        Err(ReadError::UnknownFormat(_))
    ));
}

#[test]
fn formats_follow_the_extension() {
    assert_eq!(Format::from_path("1abc.cif"), Some(Format::Mmcif));
    assert_eq!(Format::from_path("1ABC.MMCIF"), Some(Format::Mmcif));
    assert_eq!(Format::from_path("1abc.bcif"), Some(Format::Bcif));
    assert_eq!(Format::from_path("pdb1abc.ent"), Some(Format::Pdb));
    assert_eq!(Format::from_path("1abc.pdb"), Some(Format::Pdb));
    assert_eq!(Format::from_path("1abc.xyz"), None);
    assert_eq!(Format::from_path("1abc"), None);
}
