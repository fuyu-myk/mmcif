use std::fmt::Write;

use super::BcifError;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum FloatWidth {
    Single,
    Double,
}

impl FloatWidth {
    pub fn round(self, value: f64) -> f64 {
        match self {
            FloatWidth::Single => f64::from(value as f32),
            FloatWidth::Double => value,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) enum Array<'v> {
    Bytes(&'v [u8]),
    Ints(Vec<i64>),
    Floats(Vec<f64>, FloatWidth),
    Decimals(Vec<i64>, usize),
    Strings(Vec<&'v str>),
}

impl<'v> Array<'v> {
    pub fn len(&self) -> usize {
        match self {
            Array::Bytes(bytes) => bytes.len(),
            Array::Ints(values) | Array::Decimals(values, _) => values.len(),
            Array::Floats(values, _) => values.len(),
            Array::Strings(values) => values.len(),
        }
    }

    pub fn kind_name(&self) -> &'static str {
        match self {
            Array::Bytes(_) => "raw byte",
            Array::Ints(_) => "integer",
            Array::Floats(..) => "floating-point",
            Array::Decimals(..) => "fixed-point",
            Array::Strings(_) => "string",
        }
    }

    pub fn into_ints(self, encoding: &'static str) -> Result<Vec<i64>, BcifError> {
        match self {
            Array::Ints(values) => Ok(values),
            Array::Bytes(bytes) => Ok(bytes.iter().map(|&byte| i64::from(byte)).collect()),
            other => Err(BcifError::UnexpectedArray {
                encoding,
                found: other.kind_name(),
            }),
        }
    }

    pub fn text<'s>(&'s self, row: usize, buffer: &'s mut String) -> Option<&'s str> {
        buffer.clear();
        match self {
            Array::Strings(values) => return values.get(row).copied(),
            Array::Bytes(bytes) => write_number(buffer, bytes.get(row)?),
            Array::Ints(values) => write_number(buffer, values.get(row)?),
            Array::Floats(values, FloatWidth::Single) => {
                write_number(buffer, &(*values.get(row)? as f32));
            }
            Array::Floats(values, FloatWidth::Double) => write_number(buffer, values.get(row)?),
            Array::Decimals(values, decimals) => {
                write_decimal(buffer, *values.get(row)?, *decimals)
            }
        }
        Some(buffer)
    }
}

fn write_number(buffer: &mut String, value: &impl std::fmt::Display) {
    write!(buffer, "{value}").expect("writing to a String cannot fail");
}

fn write_decimal(buffer: &mut String, scaled: i64, decimals: usize) {
    if scaled < 0 {
        buffer.push('-');
    }
    write!(
        buffer,
        "{:0>width$}",
        scaled.unsigned_abs(),
        width = decimals + 1
    )
    .expect("writing to a String cannot fail");
    if decimals > 0 {
        buffer.insert(buffer.len() - decimals, '.');
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn texts(array: &Array<'_>) -> Vec<String> {
        let mut buffer = String::new();
        (0..array.len())
            .map(|row| array.text(row, &mut buffer).unwrap().to_owned())
            .collect()
    }

    #[test]
    fn writes_fixed_decimals_from_scaled_integers() {
        let array = Array::Decimals(vec![12345, -500, 0, 7, -7, 1_000_000], 3);
        assert_eq!(
            texts(&array),
            ["12.345", "-0.500", "0.000", "0.007", "-0.007", "1000.000"]
        );
        assert_eq!(texts(&Array::Decimals(vec![-42, 5], 0)), ["-42", "5"]);
        assert_eq!(
            texts(&Array::Decimals(vec![i64::MIN], 2)),
            ["-92233720368547758.08"]
        );
    }

    #[test]
    fn writes_shortest_round_trip_floats_at_source_width() {
        let single = Array::Floats(vec![f64::from(0.1f32), 1.0, -2.5], FloatWidth::Single);
        assert_eq!(texts(&single), ["0.1", "1", "-2.5"]);
        let double = Array::Floats(vec![0.1, f64::from(0.1f32)], FloatWidth::Double);
        assert_eq!(texts(&double), ["0.1", "0.10000000149011612"]);
        assert_eq!(FloatWidth::Single.round(0.1), f64::from(0.1f32));
    }

    #[test]
    fn writes_integers_bytes_and_strings() {
        assert_eq!(texts(&Array::Ints(vec![-3, 40])), ["-3", "40"]);
        assert_eq!(texts(&Array::Bytes(&[0, 255])), ["0", "255"]);
        assert_eq!(texts(&Array::Strings(vec!["a", ""])), ["a", ""]);
        assert_eq!(Array::Ints(vec![1]).text(1, &mut String::new()), None);
    }

    #[test]
    fn converts_only_integer_like_arrays() {
        assert_eq!(Array::Bytes(&[7]).into_ints("Delta").unwrap(), [7]);
        assert_eq!(
            Array::Strings(vec!["x"]).into_ints("Delta"),
            Err(BcifError::UnexpectedArray {
                encoding: "Delta",
                found: "string"
            })
        );
    }
}
