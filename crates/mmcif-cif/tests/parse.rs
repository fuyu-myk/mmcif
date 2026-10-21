use mmcif_cif::{Document, Entry, ParseErrorKind, ValueKind, stream_blocks};

const SAMPLE: &str = "\
data_1ABC
#
_entry.id 1ABC
_struct.title 'A title with spaces'
_struct.pdbx_descriptor
;Multi-line
description
;
#
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.pdbx_PDB_ins_code
ATOM 1 N . ?
ATOM 2 \"O5'\" A ?
#
data_second
_cell.length_a 10.0
";

fn parse(text: &str) -> Document {
    Document::parse(text).unwrap()
}

fn error_kind(text: &str) -> ParseErrorKind {
    Document::parse(text).unwrap_err().kind
}

#[test]
fn reads_blocks_in_order() {
    let document = parse(SAMPLE);
    let names: Vec<_> = document.blocks().map(|b| b.name()).collect();
    assert_eq!(names, ["1ABC", "second"]);
    assert_eq!(document.block_named("1abc").unwrap().name(), "1ABC");
}

#[test]
fn single_row_categories_come_from_pairs() {
    let document = parse(SAMPLE);
    let block = document.block(0).unwrap();
    let structure = block.category("struct").unwrap();
    assert_eq!(structure.row_count(), 1);
    assert_eq!(structure.items(), ["title", "pdbx_descriptor"]);
    let title = structure.column("title").unwrap().get(0).unwrap();
    assert_eq!(title.as_str(), Some("A title with spaces"));
    assert_eq!(title.kind(), ValueKind::SingleQuoted);
    let description = structure.column("pdbx_descriptor").unwrap().get(0).unwrap();
    assert_eq!(description.as_str(), Some("Multi-line\ndescription"));
}

#[test]
fn loop_columns_are_case_insensitive_and_keep_nulls() {
    let document = parse(SAMPLE);
    let atoms = document.block(0).unwrap().category("_ATOM_SITE").unwrap();
    assert_eq!(atoms.name(), "atom_site");
    assert_eq!(atoms.row_count(), 2);
    let names: Vec<_> = atoms
        .column("LABEL_ATOM_ID")
        .unwrap()
        .iter()
        .map(|v| v.raw())
        .collect();
    assert_eq!(names, ["N", "O5'"]);
    let alt = atoms.column("label_alt_id").unwrap();
    assert!(alt.get(0).unwrap().is_inapplicable());
    assert_eq!(alt.get(1).unwrap().as_str(), Some("A"));
    assert!(
        atoms
            .column("pdbx_PDB_ins_code")
            .unwrap()
            .get(1)
            .unwrap()
            .is_unknown()
    );
    assert!(atoms.column("missing").is_none());
    assert!(atoms.column("id").unwrap().get(2).is_none());
}

#[test]
fn missing_categories_are_none() {
    let document = parse(SAMPLE);
    assert!(document.block(0).unwrap().category("cell").is_none());
    assert!(document.block(1).unwrap().category("cell").is_some());
}

#[test]
fn lists_category_names_in_order() {
    let document = parse(SAMPLE);
    assert_eq!(
        document.block(0).unwrap().category_names(),
        ["entry", "struct", "atom_site"]
    );
}

#[test]
fn keeps_save_frames() {
    let document = parse("data_dict\nsave_frame\n_item.name x\nsave_\n_top.level 1\n");
    let block = document.block(0).unwrap();
    let entries: Vec<_> = block.entries().collect();
    assert_eq!(entries.len(), 2);
    let Entry::Frame(frame) = entries[0] else {
        panic!("expected a save frame, got {:?}", entries[0]);
    };
    assert_eq!(frame.name(), "frame");
    assert_eq!(frame.entries().len(), 1);
}

#[test]
fn rejects_malformed_structure() {
    assert_eq!(error_kind("_a.b c"), ParseErrorKind::OutsideBlock);
    assert_eq!(error_kind("data_\n"), ParseErrorKind::EmptyBlockName);
    assert_eq!(
        error_kind("data_x\n_a.b\n_a.c d"),
        ParseErrorKind::MissingValue
    );
    assert_eq!(
        error_kind("data_x\n_a.b c d"),
        ParseErrorKind::UnexpectedValue
    );
    assert_eq!(error_kind("data_x\nloop_\n1 2"), ParseErrorKind::EmptyLoop);
    assert_eq!(
        error_kind("data_x\nloop_\n_a.x\n_a.y\n1 2 3"),
        ParseErrorKind::RaggedLoop { values: 3, tags: 2 }
    );
    assert_eq!(
        error_kind("data_x\nsave_a\nsave_b\n"),
        ParseErrorKind::NestedSaveFrame
    );
    assert_eq!(
        error_kind("data_x\nsave_\n"),
        ParseErrorKind::UnmatchedSaveEnd
    );
    assert_eq!(
        error_kind("data_x\nsave_a\n_a.b c\n"),
        ParseErrorKind::UnterminatedSaveFrame
    );
    assert_eq!(
        error_kind("data_x\nglobal_\n"),
        ParseErrorKind::Reserved("global_")
    );
}

#[test]
fn rejects_invalid_utf8_with_position() {
    let error = Document::from_bytes(b"data_x\n_a.b \xff\n".to_vec()).unwrap_err();
    assert_eq!(error.kind, ParseErrorKind::InvalidUtf8);
    assert_eq!((error.line, error.column), (2, 6));
}

#[test]
fn streams_blocks_one_at_a_time() {
    let names: Vec<String> = stream_blocks(SAMPLE)
        .map(|block| block.unwrap().block().name().to_owned())
        .collect();
    assert_eq!(names, ["1ABC", "second"]);
}

#[test]
fn stream_stops_after_an_error() {
    let mut stream = stream_blocks("data_a\n_x.y 1\ndata_b\n_x.y\ndata_c\n");
    assert!(stream.next().unwrap().is_ok());
    assert_eq!(
        stream.next().unwrap().err().unwrap().kind,
        ParseErrorKind::MissingValue
    );
    assert!(stream.next().is_none());
}

#[test]
#[allow(unsafe_code)]
fn opens_files_by_reading_or_mapping() {
    let path = std::env::temp_dir().join(format!("mmcif-cif-open-{}.cif", std::process::id()));
    std::fs::write(&path, SAMPLE).unwrap();
    let read = Document::open(&path).unwrap();
    // SAFETY: the test owns the file and does not modify it while the document is alive.
    let mapped = unsafe { Document::open_mapped(&path) }.unwrap();
    std::fs::remove_file(&path).unwrap();
    assert_eq!(read, parse(SAMPLE));
    assert_eq!(mapped, read);
}

#[test]
fn empty_input_has_no_blocks() {
    assert!(parse("").is_empty());
    assert!(parse("# only a comment\n").is_empty());
}

#[test]
fn rejects_duplicate_tags_case_insensitively() {
    assert_eq!(
        error_kind("data_x\n_a.b 1\n_A.B 2\n"),
        ParseErrorKind::DuplicateTag("_A.B".into())
    );
    assert_eq!(
        error_kind("data_x\nloop_\n_a.x\n1\nloop_\n_a.x\n2\n"),
        ParseErrorKind::DuplicateTag("_a.x".into())
    );
    assert!(Document::parse("data_x\n_a.b 1\nsave_f\n_a.b 2\nsave_\ndata_y\n_a.b 3\n").is_ok());
}

#[test]
fn loop_columns_belong_to_their_own_category() {
    let document = parse("data_x\nloop_\n_c.x\n_d.y\n1 2\n3 4\n");
    let block = document.block(0).unwrap();
    let c = block.category("c").unwrap();
    assert_eq!(c.items(), ["x"]);
    assert!(c.column("y").is_none());
    let d = block.category("d").unwrap();
    let ys: Vec<_> = d.column("y").unwrap().iter().map(|v| v.raw()).collect();
    assert_eq!(ys, ["2", "4"]);
    assert_eq!(block.category_names(), ["c", "d"]);
}
