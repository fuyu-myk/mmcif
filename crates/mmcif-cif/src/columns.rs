use std::fmt::Write;
use std::sync::OnceLock;

use crate::document::{Cell, Span};
use crate::value::{Value, ValueKind};

/// A value of a numeric BinaryCIF column, as returned by
/// [`Column::typed`](crate::Column::typed).
///
/// Integer columns, and fixed-point columns with a factor of 1, give
/// [`Int`](TypedValue::Int). Other numeric columns give
/// [`Float`](TypedValue::Float).
///
/// # Examples
///
/// ```
/// use mmcif_cif::{Document, TypedValue};
///
/// # let bytes = b"\x81\xaadataBlocks\x91\x82\xa6header\xa41ABC\xaacategories\x91\x83\
/// #     \xa4name\xa5_atom\xa8rowCount\x02\xa7columns\x91\x82\xa4name\xa2id\xa4data\
/// #     \x82\xa4data\xc4\x02\x01\x02\xa8encoding\x91\x82\xa4kind\xa9ByteArray\xa4type\x01";
/// let document = Document::from_bcif(bytes)?;
/// let atoms = document.block(0).unwrap().category("atom").unwrap();
/// let ids = atoms.column("id").unwrap();
///
/// let total: i64 = (0..ids.len())
///     .filter_map(|row| match ids.typed(row) {
///         Some(TypedValue::Int(id)) => Some(id),
///         _ => None,
///     })
///     .sum();
/// assert_eq!(total, 3);
/// # Ok::<(), mmcif_cif::BcifError>(())
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum TypedValue {
    /// An integer.
    Int(i64),
    /// A floating-point number.
    Float(f64),
    /// The inapplicable null, written `.`.
    Inapplicable,
    /// The unknown null, written `?`.
    Unknown,
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) enum Numbers {
    Int(Vec<i64>),
    Decimal { scaled: Vec<i64>, decimals: u8 },
    Double(Vec<f64>),
    Single(Vec<f32>),
}

const EXACT_POWERS_OF_TEN: [f64; 23] = [
    1e0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9, 1e10, 1e11, 1e12, 1e13, 1e14, 1e15, 1e16,
    1e17, 1e18, 1e19, 1e20, 1e21, 1e22,
];
const EXACT_INTEGER_LIMIT: u64 = 1 << 53;

impl Numbers {
    pub fn len(&self) -> usize {
        match self {
            Numbers::Int(values) => values.len(),
            Numbers::Decimal { scaled, .. } => scaled.len(),
            Numbers::Double(values) => values.len(),
            Numbers::Single(values) => values.len(),
        }
    }

    fn write(&self, row: usize, out: &mut String) {
        match self {
            Numbers::Int(values) => write_display(out, &values[row]),
            Numbers::Decimal { scaled, decimals } => {
                write_decimal(out, scaled[row], *decimals as usize)
            }
            Numbers::Double(values) => write_display(out, &values[row]),
            Numbers::Single(values) => write_display(out, &values[row]),
        }
    }

    fn typed(&self, row: usize) -> TypedValue {
        match self {
            Numbers::Int(values) => TypedValue::Int(values[row]),
            Numbers::Decimal {
                scaled,
                decimals: 0,
            } => TypedValue::Int(scaled[row]),
            Numbers::Decimal { scaled, decimals } => {
                TypedValue::Float(decimal_value(scaled[row], *decimals))
            }
            Numbers::Double(values) => TypedValue::Float(values[row]),
            Numbers::Single(values) => TypedValue::Float(parse_written(&values[row])),
        }
    }
}

fn write_display(out: &mut String, value: &impl std::fmt::Display) {
    write!(out, "{value}").expect("writing to a String cannot fail");
}

pub(crate) fn write_decimal(out: &mut String, scaled: i64, decimals: usize) {
    if scaled < 0 {
        out.push('-');
    }
    write!(
        out,
        "{:0>width$}",
        scaled.unsigned_abs(),
        width = decimals + 1
    )
    .expect("writing to a String cannot fail");
    if decimals > 0 {
        out.insert(out.len() - decimals, '.');
    }
}

fn parse_written(value: &impl std::fmt::Display) -> f64 {
    let mut text = String::new();
    write_display(&mut text, value);
    text.parse().expect("formatted numbers parse")
}

fn decimal_value(scaled: i64, decimals: u8) -> f64 {
    match EXACT_POWERS_OF_TEN.get(decimals as usize) {
        Some(power) if scaled.unsigned_abs() <= EXACT_INTEGER_LIMIT => scaled as f64 / power,
        _ => {
            let mut text = String::new();
            write_decimal(&mut text, scaled, decimals as usize);
            text.parse().expect("formatted decimals parse")
        }
    }
}

#[derive(Debug)]
struct TextCache {
    arena: String,
    ends: Vec<usize>,
}

#[derive(Debug)]
pub(crate) struct NumericColumn {
    numbers: Numbers,
    mask: Option<Vec<ValueKind>>,
    text: OnceLock<TextCache>,
}

impl Clone for NumericColumn {
    fn clone(&self) -> Self {
        NumericColumn {
            numbers: self.numbers.clone(),
            mask: self.mask.clone(),
            text: OnceLock::new(),
        }
    }
}

impl NumericColumn {
    pub fn new(numbers: Numbers, mask: Option<Vec<ValueKind>>) -> Self {
        NumericColumn {
            numbers,
            mask,
            text: OnceLock::new(),
        }
    }

    pub fn len(&self) -> usize {
        self.numbers.len()
    }

    fn null_kind(&self, row: usize) -> Option<ValueKind> {
        self.mask
            .as_ref()
            .and_then(|mask| mask.get(row).copied())
            .filter(|kind| kind.is_null())
    }

    pub fn typed(&self, row: usize) -> Option<TypedValue> {
        if row >= self.len() {
            return None;
        }
        Some(match self.null_kind(row) {
            Some(ValueKind::Inapplicable) => TypedValue::Inapplicable,
            Some(_) => TypedValue::Unknown,
            None => self.numbers.typed(row),
        })
    }

    fn text_cache(&self) -> &TextCache {
        self.text.get_or_init(|| {
            let mut arena = String::new();
            let mut ends = Vec::with_capacity(self.len());
            for row in 0..self.len() {
                if self.null_kind(row).is_none() {
                    self.numbers.write(row, &mut arena);
                }
                ends.push(arena.len());
            }
            TextCache { arena, ends }
        })
    }

    pub fn value(&self, row: usize) -> Option<Value<'_>> {
        if row >= self.len() {
            return None;
        }
        Some(match self.null_kind(row) {
            Some(ValueKind::Inapplicable) => Value::INAPPLICABLE,
            Some(_) => Value::UNKNOWN,
            None => {
                let cache = self.text_cache();
                let start = row
                    .checked_sub(1)
                    .map_or(0, |previous| cache.ends[previous]);
                Value::from_parts(ValueKind::Bare, &cache.arena[start..cache.ends[row]])
            }
        })
    }
}

#[derive(Debug, Clone)]
pub(crate) enum RawColumn {
    Text(Vec<Cell>),
    Numeric(NumericColumn),
}

impl RawColumn {
    pub fn len(&self) -> usize {
        match self {
            RawColumn::Text(cells) => cells.len(),
            RawColumn::Numeric(column) => column.len(),
        }
    }
}

#[derive(Debug, Clone)]
pub(crate) struct RawTable {
    pub tags: Vec<Span>,
    pub columns: Vec<RawColumn>,
    pub rows: usize,
}

#[cfg(test)]
mod tests {
    use super::*;

    fn texts(column: &NumericColumn) -> Vec<String> {
        (0..column.len())
            .map(|row| column.value(row).unwrap().raw().to_owned())
            .collect()
    }

    #[test]
    fn formats_like_the_text_path() {
        let decimals = NumericColumn::new(
            Numbers::Decimal {
                scaled: vec![12345, -500, 0, 7],
                decimals: 3,
            },
            None,
        );
        assert_eq!(texts(&decimals), ["12.345", "-0.500", "0.000", "0.007"]);
        let single = NumericColumn::new(Numbers::Single(vec![0.1, 1.0, -2.5]), None);
        assert_eq!(texts(&single), ["0.1", "1", "-2.5"]);
        let double = NumericColumn::new(Numbers::Double(vec![0.1, f64::from(0.1f32)]), None);
        assert_eq!(texts(&double), ["0.1", "0.10000000149011612"]);
    }

    #[test]
    fn typed_values_equal_parsing_the_formatted_text() {
        let cases = [
            Numbers::Decimal {
                scaled: vec![12345, -500, 1, i64::MAX, i64::MIN + 1],
                decimals: 3,
            },
            Numbers::Decimal {
                scaled: vec![5, 9_007_199_254_740_993],
                decimals: 25,
            },
            Numbers::Single(vec![0.1, 1.0e-7, 3.4e38, -0.0]),
            Numbers::Double(vec![0.1, 1e300, -0.0]),
        ];
        for numbers in cases {
            let column = NumericColumn::new(numbers, None);
            for row in 0..column.len() {
                let parsed: f64 = column.value(row).unwrap().raw().parse().unwrap();
                match column.typed(row).unwrap() {
                    TypedValue::Float(value) => assert_eq!(value.to_bits(), parsed.to_bits()),
                    other => panic!("expected a float, got {other:?}"),
                }
            }
        }
    }

    #[test]
    fn integral_decimals_are_integers() {
        let column = NumericColumn::new(
            Numbers::Decimal {
                scaled: vec![-42],
                decimals: 0,
            },
            None,
        );
        assert_eq!(column.typed(0), Some(TypedValue::Int(-42)));
        assert_eq!(column.value(0).unwrap().raw(), "-42");
    }

    #[test]
    fn masked_rows_are_null() {
        let mask = vec![ValueKind::Bare, ValueKind::Inapplicable, ValueKind::Unknown];
        let column = NumericColumn::new(Numbers::Int(vec![1, 2, 3]), Some(mask));
        assert_eq!(column.typed(0), Some(TypedValue::Int(1)));
        assert_eq!(column.typed(1), Some(TypedValue::Inapplicable));
        assert_eq!(column.typed(2), Some(TypedValue::Unknown));
        assert!(column.value(1).unwrap().is_inapplicable());
        assert!(column.value(2).unwrap().is_unknown());
        assert_eq!(column.typed(3), None);
    }
}
