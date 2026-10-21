use mmcif_cif::{Document, DocumentBuilder, Value, ValueKind};
use proptest::prelude::*;

fn round_trip(document: &Document) -> Document {
    let text = document.to_cif_string();
    Document::parse(text.clone()).unwrap_or_else(|error| panic!("{error}\n{text}"))
}

#[test]
fn preserves_quoting_and_layout() {
    let text = "data_x\n_a.plain v\n_a.single 'v'\n_a.double \"v\"\n_a.null .\n_a.unknown ?\n\
                _a.text\n;v\n;\nloop_\n_b.x\n_b.y\n1 ;mid\n'q r' \"s' t\"\n";
    let document = Document::parse(text).unwrap();
    let again = round_trip(&document);
    assert_eq!(again, document);
    let pairs = again.block(0).unwrap().category("a").unwrap();
    let kinds: Vec<_> = pairs
        .items()
        .into_iter()
        .map(|item| pairs.column(item).unwrap().get(0).unwrap().kind())
        .collect();
    assert_eq!(
        kinds,
        [
            ValueKind::Bare,
            ValueKind::SingleQuoted,
            ValueKind::DoubleQuoted,
            ValueKind::Inapplicable,
            ValueKind::Unknown,
            ValueKind::TextField,
        ]
    );
}

#[test]
fn preserves_carriage_returns_inside_text_fields() {
    let document = Document::parse("data_x\n_a.t\n;line\r\r\n;\n").unwrap();
    let value = document
        .block(0)
        .unwrap()
        .category("a")
        .unwrap()
        .column("t")
        .unwrap()
        .get(0)
        .unwrap();
    assert_eq!(value.raw(), "line\r");
    assert_eq!(round_trip(&document), document);
}

fn value_text() -> impl Strategy<Value = String> {
    prop_oneof![
        "[A-Za-z0-9.'\"_;#$\\[\\] -]{0,12}",
        "[A-Za-z0-9 ;'\"\r\n]{0,24}",
        Just(".".to_owned()),
        Just("?".to_owned()),
        Just("loop_".to_owned()),
        Just("data_x".to_owned()),
    ]
}

fn value_kind() -> impl Strategy<Value = Option<ValueKind>> {
    prop_oneof![
        Just(None),
        Just(Some(ValueKind::Bare)),
        Just(Some(ValueKind::SingleQuoted)),
        Just(Some(ValueKind::DoubleQuoted)),
        Just(Some(ValueKind::TextField)),
        Just(Some(ValueKind::Inapplicable)),
        Just(Some(ValueKind::Unknown)),
    ]
}

fn make_value(text: &str, kind: Option<ValueKind>) -> Value<'_> {
    match kind {
        Some(ValueKind::Inapplicable) => Value::INAPPLICABLE,
        Some(ValueKind::Unknown) => Value::UNKNOWN,
        Some(kind) => Value::new(kind, text).unwrap_or_else(|_| Value::text(text)),
        None => Value::text(text),
    }
}

#[derive(Debug, Clone)]
enum Item {
    Pair(String, (String, Option<ValueKind>)),
    Loop(Vec<String>, Vec<(String, Option<ValueKind>)>),
}

fn item() -> impl Strategy<Value = Item> {
    let tag = "_[a-z]{1,6}\\.[a-z_\\[\\]0-9]{1,8}";
    prop_oneof![
        (tag, (value_text(), value_kind())).prop_map(|(tag, cell)| Item::Pair(tag, cell)),
        (prop::collection::vec(tag, 1..4), 1usize..5).prop_flat_map(move |(tags, rows)| {
            let count = tags.len() * rows;
            prop::collection::vec((value_text(), value_kind()), count)
                .prop_map(move |cells| Item::Loop(tags.clone(), cells))
        }),
    ]
}

fn build(blocks: &[(String, Vec<Item>)]) -> Option<Document> {
    let mut builder = DocumentBuilder::new();
    for (name, items) in blocks {
        builder.start_block(name).ok()?;
        for item in items {
            match item {
                Item::Pair(tag, (text, kind)) => builder.pair(tag, make_value(text, *kind)).ok()?,
                Item::Loop(tags, cells) => {
                    builder.start_loop(tags.iter().map(String::as_str)).ok()?;
                    for (text, kind) in cells {
                        builder.push(make_value(text, *kind)).ok()?;
                    }
                }
            }
        }
    }
    builder.finish().ok()
}

proptest! {
    #[test]
    fn built_documents_survive_write_then_parse(
        blocks in prop::collection::vec(("[A-Za-z0-9_]{1,8}", prop::collection::vec(item(), 0..6)), 1..3)
    ) {
        if let Some(document) = build(&blocks) {
            prop_assert_eq!(round_trip(&document), document);
        }
    }

    #[test]
    fn parsed_text_survives_write_then_parse(text in "(data_[a-z]{1,3}\n|_[a-z]{1,2}\\.[a-z]{1,2} |loop_ |[a-z0-9.?]{1,3} |'[a-z ']{0,4}' |\n;[a-z\n]{0,6}\n;\n|# c\n){0,24}") {
        if let Ok(document) = Document::parse(text) {
            prop_assert_eq!(round_trip(&document), document);
        }
    }
}
