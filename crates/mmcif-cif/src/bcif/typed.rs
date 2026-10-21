use super::array::{Array, FloatWidth};
use super::{BcifError, Budget};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum DataType {
    Int8,
    Int16,
    Int32,
    Uint8,
    Uint16,
    Uint32,
    Float32,
    Float64,
}

impl DataType {
    pub fn from_code(code: i64) -> Result<Self, BcifError> {
        Ok(match code {
            1 => DataType::Int8,
            2 => DataType::Int16,
            3 => DataType::Int32,
            4 => DataType::Uint8,
            5 => DataType::Uint16,
            6 => DataType::Uint32,
            32 => DataType::Float32,
            33 => DataType::Float64,
            _ => return Err(BcifError::UnsupportedDataType(code)),
        })
    }

    fn size(self) -> usize {
        match self {
            DataType::Int8 | DataType::Uint8 => 1,
            DataType::Int16 | DataType::Uint16 => 2,
            DataType::Int32 | DataType::Uint32 | DataType::Float32 => 4,
            DataType::Float64 => 8,
        }
    }
}

fn read_le<const N: usize, T>(bytes: &[u8], convert: impl Fn([u8; N]) -> T) -> Vec<T> {
    let (chunks, _) = bytes.as_chunks::<N>();
    chunks.iter().map(|&chunk| convert(chunk)).collect()
}

pub(crate) fn byte_array<'v>(
    input: Array<'v>,
    data_type: DataType,
    budget: &mut Budget,
) -> Result<Array<'v>, BcifError> {
    let Array::Bytes(bytes) = input else {
        return Err(BcifError::UnexpectedArray {
            encoding: "ByteArray",
            found: input.kind_name(),
        });
    };
    let element_size = data_type.size();
    if bytes.len() % element_size != 0 {
        return Err(BcifError::MisalignedBytes {
            len: bytes.len(),
            element_size,
        });
    }
    budget.take(bytes.len() / element_size)?;
    Ok(match data_type {
        DataType::Int8 => Array::Ints(read_le(bytes, |b| i64::from(i8::from_le_bytes(b)))),
        DataType::Int16 => Array::Ints(read_le(bytes, |b| i64::from(i16::from_le_bytes(b)))),
        DataType::Int32 => Array::Ints(read_le(bytes, |b| i64::from(i32::from_le_bytes(b)))),
        DataType::Uint8 => Array::Ints(read_le(bytes, |b| i64::from(u8::from_le_bytes(b)))),
        DataType::Uint16 => Array::Ints(read_le(bytes, |b| i64::from(u16::from_le_bytes(b)))),
        DataType::Uint32 => Array::Ints(read_le(bytes, |b| i64::from(u32::from_le_bytes(b)))),
        DataType::Float32 => Array::Floats(
            read_le(bytes, |b| f64::from(f32::from_le_bytes(b))),
            FloatWidth::Single,
        ),
        DataType::Float64 => Array::Floats(read_le(bytes, f64::from_le_bytes), FloatWidth::Double),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bcif::BcifLimits;

    fn ints(array: Result<Array<'_>, BcifError>) -> Vec<i64> {
        match array.unwrap() {
            Array::Ints(values) => values,
            other => panic!("expected integers, got {other:?}"),
        }
    }

    fn typed(bytes: &[u8], code: i64) -> Result<Array<'_>, BcifError> {
        byte_array(
            Array::Bytes(bytes),
            DataType::from_code(code)?,
            &mut Budget::new(BcifLimits::default()),
        )
    }

    #[test]
    fn reads_little_endian_typed_arrays() {
        assert_eq!(ints(typed(&[0xff, 0x02], 1)), [-1, 2]);
        assert_eq!(ints(typed(&[0x00, 0x80], 2)), [-32768]);
        assert_eq!(ints(typed(&[0xfe, 0xff, 0xff, 0xff], 3)), [-2]);
        assert_eq!(ints(typed(&[0xff], 4)), [255]);
        assert_eq!(ints(typed(&[0xff, 0xff], 5)), [65535]);
        assert_eq!(ints(typed(&[0xff; 4], 6)), [4_294_967_295]);
        assert_eq!(
            typed(&1.5f32.to_le_bytes(), 32),
            Ok(Array::Floats(vec![1.5], FloatWidth::Single))
        );
        assert_eq!(
            typed(&0.1f64.to_le_bytes(), 33),
            Ok(Array::Floats(vec![0.1], FloatWidth::Double))
        );
    }

    #[test]
    fn rejects_misaligned_and_non_byte_input() {
        assert_eq!(
            typed(&[1, 2, 3], 2),
            Err(BcifError::MisalignedBytes {
                len: 3,
                element_size: 2
            })
        );
        assert!(matches!(
            byte_array(
                Array::Ints(vec![1]),
                DataType::Int8,
                &mut Budget::new(BcifLimits::default())
            ),
            Err(BcifError::UnexpectedArray {
                encoding: "ByteArray",
                ..
            })
        ));
        assert_eq!(typed(&[], 7), Err(BcifError::UnsupportedDataType(7)));
    }
}
