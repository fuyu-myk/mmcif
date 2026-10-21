use super::array::Array;
use super::encoding::decode;
use super::msgpack::MsgValue;
use super::record::Record;
use super::{BcifError, Budget};

pub(crate) struct StringArray<'v> {
    data_encoding: &'v [MsgValue<'v>],
    string_data: &'v str,
    offset_encoding: &'v [MsgValue<'v>],
    offsets: &'v [u8],
}

impl<'v> StringArray<'v> {
    pub fn parse(record: &Record<'v>) -> Result<Self, BcifError> {
        Ok(StringArray {
            data_encoding: record.list("dataEncoding")?,
            string_data: record.str("stringData")?,
            offset_encoding: record.list("offsetEncoding")?,
            offsets: record.bin("offsets")?,
        })
    }

    pub fn apply(self, input: Array<'v>, budget: &mut Budget) -> Result<Array<'v>, BcifError> {
        let Array::Bytes(data) = input else {
            return Err(BcifError::UnexpectedArray {
                encoding: "StringArray",
                found: input.kind_name(),
            });
        };
        let offsets =
            decode(self.offsets, self.offset_encoding, budget)?.into_ints("StringArray offsets")?;
        let indices = decode(data, self.data_encoding, budget)?.into_ints("StringArray indices")?;
        let strings = split(self.string_data, &offsets, budget)?;
        budget.take(indices.len())?;
        indices
            .into_iter()
            .map(|index| lookup(&strings, index))
            .collect::<Result<_, _>>()
            .map(Array::Strings)
    }
}

enum Positions {
    Bytes,
    Table(Vec<Option<usize>>),
}

impl Positions {
    fn new(text: &str, last_offset: Option<i64>) -> Self {
        if text.is_ascii() {
            return Positions::Bytes;
        }
        let chars = text.chars().count();
        let counts_code_points =
            chars != text.encode_utf16().count() && last_offset == i64::try_from(chars).ok();
        let mut table = Vec::with_capacity(text.len() + 1);
        for (byte, ch) in text.char_indices() {
            table.push(Some(byte));
            if !counts_code_points && ch.len_utf16() == 2 {
                table.push(None);
            }
        }
        table.push(Some(text.len()));
        Positions::Table(table)
    }

    fn byte(&self, offset: i64) -> Option<usize> {
        let offset = usize::try_from(offset).ok()?;
        match self {
            Positions::Bytes => Some(offset),
            Positions::Table(table) => table.get(offset).copied().flatten(),
        }
    }
}

fn split<'v>(
    text: &'v str,
    offsets: &[i64],
    budget: &mut Budget,
) -> Result<Vec<&'v str>, BcifError> {
    budget.take(offsets.len().saturating_sub(1))?;
    let positions = Positions::new(text, offsets.last().copied());
    offsets
        .windows(2)
        .map(|pair| {
            let (start, end) = (pair[0], pair[1]);
            positions
                .byte(start)
                .zip(positions.byte(end))
                .and_then(|(from, to)| text.get(from..to))
                .ok_or(BcifError::InvalidStringOffsets { start, end })
        })
        .collect()
}

fn lookup<'v>(strings: &[&'v str], index: i64) -> Result<&'v str, BcifError> {
    if index < 0 {
        return Ok("");
    }
    usize::try_from(index)
        .ok()
        .and_then(|index| strings.get(index).copied())
        .ok_or(BcifError::InvalidStringIndex(index))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bcif::test_pack::{INT32, Pack, byte_array, decode_texts, i32_bytes, string_array};

    fn strings(text: &str, offsets: &[i32], indices: &[i32]) -> Result<Vec<String>, BcifError> {
        decode_texts(
            &i32_bytes(indices),
            vec![string_array(text, offsets, vec![byte_array(INT32)])],
        )
    }

    #[test]
    fn looks_up_strings_by_index() {
        assert_eq!(
            strings("ALAGLYHOH", &[0, 3, 6, 9], &[0, 1, 0, 2, -1]).unwrap(),
            ["ALA", "GLY", "ALA", "HOH", ""]
        );
        assert_eq!(strings("", &[0], &[-1, -1]).unwrap(), ["", ""]);
    }

    #[test]
    fn reads_utf16_offsets_like_molstar() {
        assert_eq!(
            strings("é𝔸x", &[0, 1, 3, 4], &[2, 1, 0]).unwrap(),
            ["x", "𝔸", "é"]
        );
    }

    #[test]
    fn reads_code_point_offsets_like_python_encoders() {
        assert_eq!(
            strings("é𝔸x", &[0, 1, 2, 3], &[2, 1, 0]).unwrap(),
            ["x", "𝔸", "é"]
        );
    }

    #[test]
    fn rejects_offsets_that_do_not_fit() {
        let cases: &[(&str, &[i32], (i64, i64))] = &[
            ("abc", &[0, 5], (0, 5)),
            ("abc", &[0, 2, 1], (2, 1)),
            ("abc", &[-1, 2], (-1, 2)),
            ("𝔸x", &[0, 1, 3], (0, 1)),
            ("éa", &[0, 2, 3], (2, 3)),
        ];
        for &(text, offsets, (start, end)) in cases {
            assert_eq!(
                strings(text, offsets, &[0]),
                Err(BcifError::InvalidStringOffsets { start, end }),
                "{text:?} {offsets:?}"
            );
        }
    }

    #[test]
    fn rejects_indices_past_the_last_string() {
        assert_eq!(
            strings("ab", &[0, 1, 2], &[0, 2]),
            Err(BcifError::InvalidStringIndex(2))
        );
    }

    #[test]
    fn requires_raw_bytes_as_input() {
        let encodings = vec![
            string_array("a", &[0, 1], vec![byte_array(INT32)]),
            byte_array(INT32),
        ];
        assert_eq!(
            decode_texts(&i32_bytes(&[0]), encodings),
            Err(BcifError::UnexpectedArray {
                encoding: "StringArray",
                found: "integer"
            })
        );
        let missing = vec![crate::bcif::test_pack::obj([(
            "kind",
            Pack::Str("StringArray".into()),
        )])];
        assert_eq!(
            decode_texts(&[], missing),
            Err(BcifError::MissingField("dataEncoding"))
        );
    }
}
