use super::array::{Array, FloatWidth};
use super::msgpack::MsgValue;
use super::numeric::{
    Sentinels, delta, fixed_point, integer_packing, interval_quantization, run_length,
};
use super::record::Record;
use super::strings::StringArray;
use super::typed::{DataType, byte_array};
use super::{BcifError, Budget};

fn data_type(record: &Record<'_>, key: &'static str) -> Result<DataType, BcifError> {
    DataType::from_code(record.int(key)?)
}

fn float_width(record: &Record<'_>) -> Result<FloatWidth, BcifError> {
    let code = record.int("srcType")?;
    match DataType::from_code(code)? {
        DataType::Float32 => Ok(FloatWidth::Single),
        DataType::Float64 => Ok(FloatWidth::Double),
        _ => Err(BcifError::UnsupportedDataType(code)),
    }
}

enum Encoding<'v> {
    ByteArray(DataType),
    FixedPoint {
        factor: f64,
        width: FloatWidth,
    },
    IntervalQuantization {
        min: f64,
        max: f64,
        steps: i64,
        width: FloatWidth,
    },
    RunLength {
        size: usize,
    },
    Delta {
        origin: i64,
    },
    IntegerPacking {
        sentinels: Sentinels,
        size: usize,
    },
    StringArray(StringArray<'v>),
}

impl<'v> Encoding<'v> {
    fn parse(value: &'v MsgValue<'v>) -> Result<Self, BcifError> {
        let record = Record::new(value, "encoding")?;
        let encoding = match record.str("kind")? {
            "ByteArray" => Encoding::ByteArray(data_type(&record, "type")?),
            "FixedPoint" => Encoding::FixedPoint {
                factor: record.number("factor")?,
                width: float_width(&record)?,
            },
            "IntervalQuantization" => Encoding::IntervalQuantization {
                min: record.number("min")?,
                max: record.number("max")?,
                steps: record.int("numSteps")?,
                width: float_width(&record)?,
            },
            "RunLength" => {
                data_type(&record, "srcType")?;
                Encoding::RunLength {
                    size: record.count("srcSize")?,
                }
            }
            "Delta" => {
                data_type(&record, "srcType")?;
                Encoding::Delta {
                    origin: record.int("origin")?,
                }
            }
            "IntegerPacking" => Encoding::IntegerPacking {
                sentinels: Sentinels::new(record.int("byteCount")?, record.bool("isUnsigned")?)?,
                size: record.count("srcSize")?,
            },
            "StringArray" => Encoding::StringArray(StringArray::parse(&record)?),
            other => return Err(BcifError::UnknownEncoding(other.to_owned())),
        };
        Ok(encoding)
    }

    fn apply(self, input: Array<'v>, budget: &mut Budget) -> Result<Array<'v>, BcifError> {
        match self {
            Encoding::ByteArray(data_type) => byte_array(input, data_type, budget),
            Encoding::FixedPoint { factor, width } => {
                fixed_point(input.into_ints("FixedPoint")?, factor, width)
            }
            Encoding::IntervalQuantization {
                min,
                max,
                steps,
                width,
            } => interval_quantization(
                input.into_ints("IntervalQuantization")?,
                min,
                max,
                steps,
                width,
            ),
            Encoding::RunLength { size } => run_length(input.into_ints("RunLength")?, size, budget),
            Encoding::Delta { origin } => delta(input.into_ints("Delta")?, origin),
            Encoding::IntegerPacking { sentinels, size } => {
                integer_packing(input.into_ints("IntegerPacking")?, sentinels, size, budget)
            }
            Encoding::StringArray(strings) => strings.apply(input, budget),
        }
    }
}

pub(crate) fn decode<'v>(
    data: &'v [u8],
    encodings: &'v [MsgValue<'v>],
    budget: &mut Budget,
) -> Result<Array<'v>, BcifError> {
    encodings
        .iter()
        .rev()
        .try_fold(Array::Bytes(data), |current, encoding| {
            Encoding::parse(encoding)?.apply(current, budget)
        })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bcif::BcifLimits;
    use crate::bcif::test_pack::{
        FLOAT32, INT8, INT16, INT32, Pack, byte_array, decode_texts, decode_texts_with, i8_bytes,
        i16_bytes, i32_bytes, obj, text,
    };
    use crate::bcif::test_pack::{
        delta as delta_step, fixed_point as fixed_point_step, integer_packing as packing_step,
        run_length as run_length_step,
    };

    #[test]
    fn decodes_encodings_in_reverse_order() {
        let encodings = vec![
            delta_step(1),
            run_length_step(5),
            packing_step(1, false, 4),
            byte_array(INT8),
        ];
        assert_eq!(
            decode_texts(&i8_bytes(&[0, 1, 1, 4]), encodings).unwrap(),
            ["1", "2", "3", "4", "5"]
        );
    }

    #[test]
    fn composes_delta_run_length_and_packing() {
        let encodings = vec![
            delta_step(10),
            run_length_step(301),
            packing_step(1, false, 4),
            byte_array(INT8),
        ];
        let values = decode_texts(&i8_bytes(&[0, 1, 1, 127, 127, 46]), encodings).unwrap();
        assert_eq!(values.len(), 301);
        assert_eq!(values[..3], ["10", "11", "12"]);
        assert_eq!(values[300], "310");
    }

    #[test]
    fn decodes_coordinates_like_molstar() {
        let encodings = vec![
            fixed_point_step(Pack::Int(1000), FLOAT32),
            delta_step(17047),
            packing_step(2, false, 3),
            byte_array(INT16),
        ];
        let packed = i16_bytes(&[0, 32767, 32767, 2, -1000]);
        assert_eq!(
            decode_texts(&packed, encodings).unwrap(),
            ["17.047", "82.583", "81.583"]
        );
    }

    #[test]
    fn accepts_float_factors_and_raw_bytes() {
        let encodings = vec![
            fixed_point_step(Pack::Float(100.0), FLOAT32),
            byte_array(INT32),
        ];
        assert_eq!(
            decode_texts(&i32_bytes(&[100, -5]), encodings).unwrap(),
            ["1.00", "-0.05"]
        );
        assert_eq!(decode_texts(&[3, 4], Vec::new()).unwrap(), ["3", "4"]);
    }

    #[test]
    fn rejects_unknown_kinds_and_types() {
        let unknown = obj([("kind", text("Zstd"))]);
        assert_eq!(
            decode_texts(&[], vec![unknown]),
            Err(BcifError::UnknownEncoding("Zstd".into()))
        );
        assert_eq!(
            decode_texts(&[], vec![byte_array(9)]),
            Err(BcifError::UnsupportedDataType(9))
        );
        assert_eq!(
            decode_texts(&[], vec![fixed_point_step(Pack::Int(10), INT32)]),
            Err(BcifError::UnsupportedDataType(INT32))
        );
        let missing = obj([
            ("kind", text("FixedPoint")),
            ("srcType", Pack::Int(FLOAT32)),
        ]);
        assert_eq!(
            decode_texts(&[], vec![missing]),
            Err(BcifError::MissingField("factor"))
        );
        assert!(matches!(
            decode_texts(&[], vec![Pack::Int(1)]),
            Err(BcifError::UnexpectedType {
                field: "encoding",
                ..
            })
        ));
        let mistyped = obj([
            ("kind", text("Delta")),
            ("srcType", Pack::Int(INT32)),
            ("origin", text("0")),
        ]);
        assert!(matches!(
            decode_texts(&[], vec![mistyped]),
            Err(BcifError::UnexpectedType {
                field: "origin",
                ..
            })
        ));
    }

    #[test]
    fn rejects_steps_applied_to_the_wrong_array() {
        let encodings = vec![delta_step(0), byte_array(FLOAT32)];
        assert_eq!(
            decode_texts(&1.0f32.to_le_bytes(), encodings),
            Err(BcifError::UnexpectedArray {
                encoding: "Delta",
                found: "floating-point"
            })
        );
    }

    #[test]
    fn enforces_the_element_limit() {
        let limits = BcifLimits {
            max_elements: 5,
            ..BcifLimits::default()
        };
        assert_eq!(
            decode_texts_with(&i32_bytes(&[0; 6]), vec![byte_array(INT32)], limits),
            Err(BcifError::LimitExceeded { limit: 5 })
        );
        let encodings = vec![run_length_step(2 * i64::from(i32::MAX)), byte_array(INT32)];
        assert_eq!(
            decode_texts(&i32_bytes(&[1, i32::MAX, 2, i32::MAX]), encodings),
            Err(BcifError::LimitExceeded {
                limit: BcifLimits::default().max_elements
            })
        );
    }
}
