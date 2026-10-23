use std::collections::HashMap;
use std::collections::hash_map::Entry;

use super::array::{Array, FloatWidth};
use super::encoding::decode;
use super::msgpack::{self, MsgValue};
use super::record::Record;
use super::{BcifError, BcifLimits, Budget};
use crate::builder::{BuildError, DocumentBuilder};
use crate::columns::{Numbers, NumericColumn, RawColumn};
use crate::document::{Cell, Document};
use crate::value::{Value, ValueKind};

pub(crate) fn decode_document(bytes: &[u8], limits: BcifLimits) -> Result<Document, BcifError> {
    let mut budget = Budget::new(limits);
    let root = msgpack::parse(bytes, &mut budget)?;
    let root = Record::new(&root, "file")?;
    let mut builder = DocumentBuilder::new();
    for block in root.list("dataBlocks")? {
        let block = Record::new(block, "dataBlocks")?;
        builder
            .start_block(block.str("header")?)
            .map_err(document_error)?;
        for category in block.list("categories")? {
            add_category(
                &mut builder,
                Record::new(category, "categories")?,
                &mut budget,
                limits.max_text_bytes,
            )?;
        }
    }
    builder.finish().map_err(document_error)
}

fn document_error(error: BuildError) -> BcifError {
    match error {
        BuildError::InvalidBlockName(name) => BcifError::InvalidBlockName(name),
        BuildError::InvalidTag(tag) => BcifError::InvalidTag(tag),
        BuildError::DuplicateTag(tag) => BcifError::DuplicateTag(tag),
        BuildError::TooLarge => BcifError::TooLarge,
        error @ (BuildError::NoBlock
        | BuildError::NoLoop
        | BuildError::InvalidValue(_)
        | BuildError::EmptyLoop
        | BuildError::LoopWithoutRows
        | BuildError::RaggedLoop { .. }) => unreachable!(
            "BinaryCIF decoding starts a block first, never opens a loop, adds whole \
             columns of equal length and writes values whose kind holds their text: {error}"
        ),
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Presence {
    Present,
    Inapplicable,
    Unknown,
}

impl Presence {
    fn from_mask(value: i64) -> Result<Self, BcifError> {
        match value {
            0 => Ok(Presence::Present),
            1 => Ok(Presence::Inapplicable),
            2 => Ok(Presence::Unknown),
            _ => Err(BcifError::InvalidMask(value)),
        }
    }
}

struct DecodedColumn<'v> {
    tag: String,
    values: Array<'v>,
    mask: Option<Vec<Presence>>,
}

impl<'v> DecodedColumn<'v> {
    fn decode(
        category: &str,
        column: &'v MsgValue<'v>,
        rows: usize,
        budget: &mut Budget,
    ) -> Result<Self, BcifError> {
        let column = Record::new(column, "columns")?;
        let tag = format!("_{category}.{}", column.str("name")?);
        let values = decode_encoded(column.record("data")?, budget)?;
        check_rows("column data", values.len(), rows)?;
        let mask = column
            .optional("mask")
            .map(|mask| decode_mask(Record::new(mask, "mask")?, rows, budget))
            .transpose()?;
        Ok(DecodedColumn { tag, values, mask })
    }

    fn value<'s>(&'s self, row: usize, buffer: &'s mut String) -> Result<Value<'s>, BcifError> {
        match self.mask.as_ref().and_then(|mask| mask.get(row)) {
            Some(Presence::Inapplicable) => return Ok(Value::INAPPLICABLE),
            Some(Presence::Unknown) => return Ok(Value::UNKNOWN),
            Some(Presence::Present) | None => {}
        }
        let text = self
            .values
            .text(row, buffer)
            .ok_or(BcifError::LengthMismatch {
                what: "column data",
                expected: row + 1,
                actual: self.values.len(),
            })?;
        if text.contains("\n;") {
            return Err(BcifError::UnrepresentableValue {
                tag: self.tag.clone(),
            });
        }
        Ok(Value::text(text))
    }
}

fn check_rows(what: &'static str, actual: usize, expected: usize) -> Result<(), BcifError> {
    if actual == expected {
        Ok(())
    } else {
        Err(BcifError::LengthMismatch {
            what,
            expected,
            actual,
        })
    }
}

fn decode_encoded<'v>(encoded: Record<'v>, budget: &mut Budget) -> Result<Array<'v>, BcifError> {
    decode(encoded.bin("data")?, encoded.list("encoding")?, budget)
}

fn decode_mask(
    mask: Record<'_>,
    rows: usize,
    budget: &mut Budget,
) -> Result<Vec<Presence>, BcifError> {
    let values = decode_encoded(mask, budget)?.into_ints("mask")?;
    check_rows("mask", values.len(), rows)?;
    values.into_iter().map(Presence::from_mask).collect()
}

fn check_text(builder: &DocumentBuilder, limit: usize) -> Result<(), BcifError> {
    if builder.text_len() > limit {
        return Err(BcifError::TextLimitExceeded { limit });
    }
    Ok(())
}

fn add_category(
    builder: &mut DocumentBuilder,
    category: Record<'_>,
    budget: &mut Budget,
    max_text_bytes: usize,
) -> Result<(), BcifError> {
    let name = category.str("name")?;
    let name = name.strip_prefix('_').unwrap_or(name);
    let rows = category.count("rowCount")?;
    let columns = category
        .list("columns")?
        .iter()
        .map(|column| DecodedColumn::decode(name, column, rows, budget))
        .collect::<Result<Vec<_>, _>>()?;
    if rows == 0 || columns.is_empty() {
        return Ok(());
    }
    let mut buffer = String::new();
    if rows == 1 {
        for column in &columns {
            builder
                .pair(&column.tag, column.value(0, &mut buffer)?)
                .map_err(document_error)?;
            check_text(builder, max_text_bytes)?;
        }
        return Ok(());
    }
    let tags: Vec<String> = columns.iter().map(|column| column.tag.clone()).collect();
    let raw_columns = columns
        .into_iter()
        .map(|column| column.into_raw_column(builder))
        .collect::<Result<Vec<_>, _>>()?;
    builder
        .push_table(
            &tags.iter().map(String::as_str).collect::<Vec<_>>(),
            raw_columns,
        )
        .map_err(document_error)?;
    check_text(builder, max_text_bytes)
}

impl Presence {
    fn kind(self) -> ValueKind {
        match self {
            Presence::Present => ValueKind::Bare,
            Presence::Inapplicable => ValueKind::Inapplicable,
            Presence::Unknown => ValueKind::Unknown,
        }
    }
}

impl DecodedColumn<'_> {
    fn into_raw_column(self, builder: &mut DocumentBuilder) -> Result<RawColumn, BcifError> {
        let DecodedColumn { tag, values, mask } = self;
        let numbers = match values {
            Array::Ints(values) => Numbers::Int(values),
            Array::Bytes(bytes) => {
                Numbers::Int(bytes.iter().map(|&byte| i64::from(byte)).collect())
            }
            Array::Decimals(scaled, decimals) => Numbers::Decimal {
                scaled,
                decimals: u8::try_from(decimals)
                    .map_err(|_| BcifError::Overflow("fixed-point decimals"))?,
            },
            Array::Floats(values, FloatWidth::Single) => {
                Numbers::Single(values.iter().map(|&value| value as f32).collect())
            }
            Array::Floats(values, FloatWidth::Double) => Numbers::Double(values),
            Array::Strings(strings) => {
                return text_column(builder, &tag, &strings, mask.as_deref());
            }
        };
        let mask = mask.map(|mask| mask.into_iter().map(Presence::kind).collect());
        Ok(RawColumn::Numeric(NumericColumn::new(numbers, mask)))
    }
}

fn text_column(
    builder: &mut DocumentBuilder,
    tag: &str,
    strings: &[&str],
    mask: Option<&[Presence]>,
) -> Result<RawColumn, BcifError> {
    let mut stored: HashMap<(usize, usize), Cell> = HashMap::new();
    let mut nulls: [Option<Cell>; 2] = [None, None];
    let mut cells = Vec::with_capacity(strings.len());
    for (row, text) in strings.iter().enumerate() {
        let presence = mask
            .and_then(|mask| mask.get(row))
            .copied()
            .unwrap_or(Presence::Present);
        let cell = match presence {
            Presence::Inapplicable | Presence::Unknown => {
                let (slot, value) = if presence == Presence::Inapplicable {
                    (0, Value::INAPPLICABLE)
                } else {
                    (1, Value::UNKNOWN)
                };
                match nulls[slot] {
                    Some(cell) => cell,
                    None => *nulls[slot].insert(builder.text_cell(value).map_err(document_error)?),
                }
            }
            Presence::Present => match stored.entry((text.as_ptr() as usize, text.len())) {
                Entry::Occupied(entry) => *entry.get(),
                Entry::Vacant(entry) => {
                    if text.contains("\n;") {
                        return Err(BcifError::UnrepresentableValue {
                            tag: tag.to_owned(),
                        });
                    }
                    *entry.insert(
                        builder
                            .text_cell(Value::text(text))
                            .map_err(document_error)?,
                    )
                }
            },
        };
        cells.push(cell);
    }
    Ok(RawColumn::Text(cells))
}

#[cfg(test)]
mod tests {
    use crate::bcif::msgpack::MAX_DEPTH;
    use crate::bcif::test_files::{category, column, encoded, file, sample, strings};
    use crate::bcif::test_pack::{INT8, Pack, byte_array, i8_bytes, obj, text};
    use crate::{BcifError, BcifLimits, BuildError, Document, Error, ValueKind};

    fn texts(document: &Document, category: &str, item: &str) -> Vec<(String, ValueKind)> {
        let block = document.block(0).unwrap();
        let column = block.category(category).unwrap().column(item).unwrap();
        column
            .iter()
            .map(|value| (value.raw().to_owned(), value.kind()))
            .collect()
    }

    fn raws(document: &Document, category: &str, item: &str) -> Vec<String> {
        texts(document, category, item)
            .into_iter()
            .map(|(raw, _)| raw)
            .collect()
    }

    #[test]
    fn builds_documents_with_pairs_loops_and_masks() {
        let document = Document::from_bcif(&sample()).unwrap();
        assert_eq!(document.len(), 1);
        let block = document.block(0).unwrap();
        assert_eq!(block.name(), "TEST");
        assert_eq!(block.category_names(), ["entry", "atom_site"]);
        assert_eq!(
            texts(&document, "entry", "id"),
            [("1ABC".into(), ValueKind::Bare)]
        );
        assert_eq!(block.category("atom_site").unwrap().row_count(), 4);
        assert_eq!(raws(&document, "atom_site", "id"), ["1", "2", "3", "4"]);
        assert_eq!(
            texts(&document, "atom_site", "label_atom_id")[2..],
            [
                ("O5'".into(), ValueKind::Bare),
                ("two words".into(), ValueKind::SingleQuoted)
            ]
        );
        let kinds: Vec<_> = texts(&document, "atom_site", "label_alt_id")
            .into_iter()
            .map(|(_, kind)| kind)
            .collect();
        assert_eq!(
            kinds,
            [
                ValueKind::Inapplicable,
                ValueKind::Inapplicable,
                ValueKind::Bare,
                ValueKind::Unknown
            ]
        );
        assert_eq!(
            raws(&document, "atom_site", "Cartn_x"),
            ["12.345", "-0.500", "0.000", "0.007"]
        );
        assert_eq!(
            raws(&document, "atom_site", "occupancy"),
            ["1", "0.5", "0.25", "0.1"]
        );
        assert_eq!(
            raws(&document, "atom_site", "B_iso_or_equiv"),
            ["10.5", "0.1", "3", "0.0000001"]
        );
        assert_eq!(
            raws(&document, "atom_site", "pdbx_formal_charge"),
            ["0", "?", "-1", "1"]
        );
        let again = Document::parse(document.to_cif_string()).unwrap();
        assert_eq!(again, document);
    }

    #[test]
    fn every_truncation_and_byte_mutation_fails_cleanly() {
        let bytes = sample();
        for len in 0..bytes.len() {
            assert!(Document::from_bcif(&bytes[..len]).is_err(), "prefix {len}");
        }
        let mut decoded = 0;
        for position in 0..bytes.len() {
            for replacement in [0x00, 0x7f, 0xc1, 0xcf, 0xdd, 0xff] {
                let mut mutated = bytes.clone();
                mutated[position] = replacement;
                decoded += usize::from(Document::from_bcif(&mutated).is_ok());
            }
        }
        assert!(decoded < bytes.len() * 6);
    }

    #[test]
    fn rejects_invalid_masks_and_lengths() {
        let data = || encoded(i8_bytes(&[1, 2]), vec![byte_array(INT8)]);
        let bad_mask = file(
            "x",
            vec![category("_a", 2, vec![column("b", data(), Some(&[0, 3]))])],
        );
        assert_eq!(
            Document::from_bcif(&bad_mask).unwrap_err(),
            BcifError::InvalidMask(3)
        );
        let short_mask = file(
            "x",
            vec![category("_a", 2, vec![column("b", data(), Some(&[0]))])],
        );
        assert_eq!(
            Document::from_bcif(&short_mask).unwrap_err(),
            BcifError::LengthMismatch {
                what: "mask",
                expected: 2,
                actual: 1
            }
        );
        let ragged = file(
            "x",
            vec![category("_a", 3, vec![column("b", data(), None)])],
        );
        assert_eq!(
            Document::from_bcif(&ragged).unwrap_err(),
            BcifError::LengthMismatch {
                what: "column data",
                expected: 3,
                actual: 2
            }
        );
    }

    #[test]
    fn rejects_text_that_cif_cannot_hold() {
        let bytes = file(
            "x",
            vec![category(
                "_a",
                1,
                vec![column("b", strings(&["x\n;y"]), None)],
            )],
        );
        assert_eq!(
            Document::from_bcif(&bytes).unwrap_err(),
            BcifError::UnrepresentableValue { tag: "_a.b".into() }
        );
    }

    #[test]
    fn reports_an_oversized_document_as_too_large() {
        assert_eq!(
            super::document_error(BuildError::TooLarge),
            BcifError::TooLarge
        );
    }

    #[test]
    fn rejects_names_and_tags_that_do_not_form_a_document() {
        assert_eq!(
            Document::from_bcif(&file("two words", Vec::new())).unwrap_err(),
            BcifError::InvalidBlockName("two words".into())
        );
        assert_eq!(
            Document::from_bcif(&file("", Vec::new())).unwrap_err(),
            BcifError::InvalidBlockName(String::new())
        );
        let spaced = file(
            "x",
            vec![category(
                "_a",
                1,
                vec![column("b c", strings(&["1"]), None)],
            )],
        );
        assert_eq!(
            Document::from_bcif(&spaced).unwrap_err(),
            BcifError::InvalidTag("_a.b c".into())
        );
        let repeated = file(
            "x",
            vec![category(
                "_a",
                2,
                vec![
                    column("b", strings(&["1", "2"]), None),
                    column("B", strings(&["3", "4"]), None),
                ],
            )],
        );
        assert_eq!(
            Document::from_bcif(&repeated).unwrap_err(),
            BcifError::DuplicateTag("_a.B".into())
        );
        let twice = file(
            "x",
            vec![
                category("_a", 1, vec![column("b", strings(&["1"]), None)]),
                category("_a", 1, vec![column("b", strings(&["2"]), None)]),
            ],
        );
        assert_eq!(
            Document::from_bcif(&twice).unwrap_err(),
            BcifError::DuplicateTag("_a.b".into())
        );
    }

    #[test]
    fn rejects_wrong_structure() {
        assert_eq!(
            Document::from_bcif(&obj([("version", text("1"))]).encode()).unwrap_err(),
            BcifError::MissingField("dataBlocks")
        );
        assert!(matches!(
            Document::from_bcif(&Pack::Array(Vec::new()).encode()).unwrap_err(),
            BcifError::UnexpectedType { field: "file", .. }
        ));
        let header = obj([(
            "dataBlocks",
            Pack::Array(vec![obj([("header", Pack::Int(1))])]),
        )]);
        assert!(matches!(
            Document::from_bcif(&header.encode()).unwrap_err(),
            BcifError::UnexpectedType {
                field: "header",
                ..
            }
        ));
        let mut nested = Pack::Nil;
        for _ in 0..=MAX_DEPTH {
            nested = Pack::Array(vec![nested]);
        }
        let deep = obj([("dataBlocks", nested)]).encode();
        assert_eq!(
            Document::from_bcif(&deep).unwrap_err(),
            BcifError::TooDeep { limit: MAX_DEPTH }
        );
        assert!(matches!(
            Document::from_bcif(&[0xdf, 0xff, 0xff, 0xff, 0xff]).unwrap_err(),
            BcifError::Truncated { .. }
        ));
    }

    #[test]
    fn honours_custom_limits_and_reports_io_errors() {
        let limits = BcifLimits {
            max_elements: 10,
            ..BcifLimits::default()
        };
        assert_eq!(
            Document::from_bcif_with_limits(&sample(), limits).unwrap_err(),
            BcifError::LimitExceeded { limit: 10 }
        );
        assert!(matches!(
            Document::open_bcif("/nonexistent/file.bcif"),
            Err(Error::Io(_))
        ));
    }

    #[test]
    fn bounds_decoded_text() {
        let limits = BcifLimits {
            max_text_bytes: 8,
            ..BcifLimits::default()
        };
        assert_eq!(
            Document::from_bcif_with_limits(&sample(), limits).unwrap_err(),
            BcifError::TextLimitExceeded { limit: 8 }
        );
    }
}
