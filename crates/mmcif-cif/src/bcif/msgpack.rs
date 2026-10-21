use super::{BcifError, Budget};

pub(crate) const MAX_DEPTH: usize = 64;

#[derive(Debug, Clone, PartialEq)]
pub(crate) enum MsgValue<'a> {
    Nil,
    Bool(bool),
    Integer(i128),
    Float(f64),
    Str(&'a str),
    Bin(&'a [u8]),
    Array(Vec<MsgValue<'a>>),
    Map(Vec<(MsgValue<'a>, MsgValue<'a>)>),
}

pub(crate) fn parse<'a>(bytes: &'a [u8], budget: &mut Budget) -> Result<MsgValue<'a>, BcifError> {
    let mut reader = Reader {
        bytes,
        offset: 0,
        budget,
    };
    let value = reader.value(0)?;
    match reader.remaining() {
        0 => Ok(value),
        count => Err(BcifError::TrailingBytes { count }),
    }
}

struct Reader<'a, 'b> {
    bytes: &'a [u8],
    offset: usize,
    budget: &'b mut Budget,
}

impl<'a> Reader<'a, '_> {
    fn remaining(&self) -> usize {
        self.bytes.len() - self.offset
    }

    fn truncated(&self) -> BcifError {
        BcifError::Truncated {
            offset: self.offset,
        }
    }

    fn take(&mut self, len: usize) -> Result<&'a [u8], BcifError> {
        let end = self
            .offset
            .checked_add(len)
            .ok_or_else(|| self.truncated())?;
        let slice = self
            .bytes
            .get(self.offset..end)
            .ok_or_else(|| self.truncated())?;
        self.offset = end;
        Ok(slice)
    }

    fn fixed<const N: usize>(&mut self) -> Result<[u8; N], BcifError> {
        let mut array = [0; N];
        array.copy_from_slice(self.take(N)?);
        Ok(array)
    }

    fn unsigned(&mut self, width: usize) -> Result<u64, BcifError> {
        let bytes = self.take(width)?;
        Ok(bytes
            .iter()
            .fold(0, |acc, &byte| (acc << 8) | u64::from(byte)))
    }

    fn length(&mut self, width: usize) -> Result<usize, BcifError> {
        let len = self.unsigned(width)?;
        usize::try_from(len).map_err(|_| self.truncated())
    }

    fn value(&mut self, depth: usize) -> Result<MsgValue<'a>, BcifError> {
        let offset = self.offset;
        let [marker] = self.fixed::<1>()?;
        let value = match marker {
            0x00..=0x7f => MsgValue::Integer(i128::from(marker)),
            0x80..=0x8f => self.map(usize::from(marker & 0x0f), depth)?,
            0x90..=0x9f => self.array(usize::from(marker & 0x0f), depth)?,
            0xa0..=0xbf => self.str(usize::from(marker & 0x1f))?,
            0xc0 => MsgValue::Nil,
            0xc2 => MsgValue::Bool(false),
            0xc3 => MsgValue::Bool(true),
            0xc4..=0xc6 => {
                let len = self.length(1 << (marker - 0xc4))?;
                MsgValue::Bin(self.take(len)?)
            }
            0xca => MsgValue::Float(f64::from(f32::from_be_bytes(self.fixed()?))),
            0xcb => MsgValue::Float(f64::from_be_bytes(self.fixed()?)),
            0xcc..=0xcf => MsgValue::Integer(i128::from(self.unsigned(1 << (marker - 0xcc))?)),
            0xd0 => MsgValue::Integer(i128::from(i8::from_be_bytes(self.fixed()?))),
            0xd1 => MsgValue::Integer(i128::from(i16::from_be_bytes(self.fixed()?))),
            0xd2 => MsgValue::Integer(i128::from(i32::from_be_bytes(self.fixed()?))),
            0xd3 => MsgValue::Integer(i128::from(i64::from_be_bytes(self.fixed()?))),
            0xd9..=0xdb => {
                let len = self.length(1 << (marker - 0xd9))?;
                self.str(len)?
            }
            0xdc | 0xdd => {
                let len = self.length(2 << (marker - 0xdc))?;
                self.array(len, depth)?
            }
            0xde | 0xdf => {
                let len = self.length(2 << (marker - 0xde))?;
                self.map(len, depth)?
            }
            0xe0..=0xff => MsgValue::Integer(i128::from(i8::from_be_bytes([marker]))),
            _ => return Err(BcifError::UnsupportedMarker { marker, offset }),
        };
        Ok(value)
    }

    fn str(&mut self, len: usize) -> Result<MsgValue<'a>, BcifError> {
        let bytes = self.take(len)?;
        std::str::from_utf8(bytes)
            .map(MsgValue::Str)
            .map_err(|_| BcifError::InvalidUtf8)
    }

    fn enter(&mut self, depth: usize, count: usize, bytes_each: usize) -> Result<usize, BcifError> {
        if depth >= MAX_DEPTH {
            return Err(BcifError::TooDeep { limit: MAX_DEPTH });
        }
        let needed = count.checked_mul(bytes_each);
        if needed.is_none_or(|needed| needed > self.remaining()) {
            return Err(self.truncated());
        }
        self.budget.take(count)?;
        Ok(depth + 1)
    }

    fn array(&mut self, len: usize, depth: usize) -> Result<MsgValue<'a>, BcifError> {
        let child = self.enter(depth, len, 1)?;
        let mut items = Vec::with_capacity(len);
        for _ in 0..len {
            items.push(self.value(child)?);
        }
        Ok(MsgValue::Array(items))
    }

    fn map(&mut self, len: usize, depth: usize) -> Result<MsgValue<'a>, BcifError> {
        let child = self.enter(depth, len, 2)?;
        let mut entries = Vec::with_capacity(len);
        for _ in 0..len {
            let key = self.value(child)?;
            entries.push((key, self.value(child)?));
        }
        Ok(MsgValue::Map(entries))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bcif::BcifLimits;
    use crate::bcif::test_pack::Pack;

    fn decode(bytes: &[u8]) -> Result<MsgValue<'_>, BcifError> {
        parse(bytes, &mut Budget::new(BcifLimits::default()))
    }

    #[test]
    fn decodes_scalars_of_every_width() {
        let cases: &[(&[u8], MsgValue<'_>)] = &[
            (&[0xc0], MsgValue::Nil),
            (&[0xc3], MsgValue::Bool(true)),
            (&[0x7f], MsgValue::Integer(127)),
            (&[0xe0], MsgValue::Integer(-32)),
            (&[0xcc, 0xff], MsgValue::Integer(255)),
            (&[0xcd, 0x01, 0x00], MsgValue::Integer(256)),
            (
                &[0xce, 0xff, 0xff, 0xff, 0xff],
                MsgValue::Integer(4_294_967_295),
            ),
            (
                &[0xcf, 0xff, 0, 0, 0, 0, 0, 0, 0],
                MsgValue::Integer(0xff << 56),
            ),
            (&[0xd0, 0x80], MsgValue::Integer(-128)),
            (&[0xd1, 0xff, 0x00], MsgValue::Integer(-256)),
            (
                &[0xd2, 0x80, 0, 0, 0],
                MsgValue::Integer(i128::from(i32::MIN)),
            ),
            (
                &[0xd3, 0x80, 0, 0, 0, 0, 0, 0, 0],
                MsgValue::Integer(i128::from(i64::MIN)),
            ),
            (&[0xca, 0x3f, 0xc0, 0, 0], MsgValue::Float(1.5)),
            (&[0xcb, 0x3f, 0xf8, 0, 0, 0, 0, 0, 0], MsgValue::Float(1.5)),
            (&[0xa2, b'h', b'i'], MsgValue::Str("hi")),
            (&[0xd9, 0x01, b'x'], MsgValue::Str("x")),
            (&[0xda, 0x00, 0x01, b'x'], MsgValue::Str("x")),
            (&[0xdb, 0, 0, 0, 0x01, b'x'], MsgValue::Str("x")),
            (&[0xc4, 0x01, 0x09], MsgValue::Bin(&[9])),
            (&[0xc5, 0x00, 0x01, 0x09], MsgValue::Bin(&[9])),
            (&[0xc6, 0, 0, 0, 0x01, 0x09], MsgValue::Bin(&[9])),
        ];
        for (bytes, expected) in cases {
            assert_eq!(&decode(bytes).unwrap(), expected, "{bytes:02x?}");
        }
    }

    #[test]
    fn decodes_containers() {
        let bytes = [0xdc, 0x00, 0x02, 0x01, 0xde, 0x00, 0x01, 0xa1, b'k', 0xc2];
        let expected = MsgValue::Array(vec![
            MsgValue::Integer(1),
            MsgValue::Map(vec![(MsgValue::Str("k"), MsgValue::Bool(false))]),
        ]);
        assert_eq!(decode(&bytes).unwrap(), expected);
        let bytes = [0xdd, 0, 0, 0, 0, 0xdf, 0, 0, 0, 0, 0x90, 0x80];
        assert_eq!(decode(&bytes[..5]).unwrap(), MsgValue::Array(Vec::new()));
        assert_eq!(decode(&bytes[5..10]).unwrap(), MsgValue::Map(Vec::new()));
    }

    #[test]
    fn rejects_truncated_and_oversized_lengths() {
        for bytes in [
            &[][..],
            &[0xcd, 0x01],
            &[0xa3, b'a'],
            &[0x92, 0x01],
            &[0xdd, 0xff, 0xff, 0xff, 0xff],
            &[0xdf, 0xff, 0xff, 0xff, 0xff, 0x01],
            &[0xc6, 0xff, 0xff, 0xff, 0xff],
            &[0xdb, 0xff, 0xff, 0xff, 0xff],
        ] {
            assert!(
                matches!(decode(bytes), Err(BcifError::Truncated { .. })),
                "{bytes:02x?}"
            );
        }
    }

    #[test]
    fn rejects_ext_types_and_trailing_bytes() {
        for bytes in [&[0xc1][..], &[0xd4, 0x01, 0x02], &[0xc7, 0x01, 0x01, 0x00]] {
            assert!(matches!(
                decode(bytes),
                Err(BcifError::UnsupportedMarker { offset: 0, .. })
            ));
        }
        assert_eq!(
            decode(&[0xc0, 0xc0]),
            Err(BcifError::TrailingBytes { count: 1 })
        );
        assert_eq!(decode(&[0xa1, 0xff]), Err(BcifError::InvalidUtf8));
    }

    #[test]
    fn caps_nesting_depth() {
        let mut bytes = vec![0x91; MAX_DEPTH];
        bytes.push(0xc0);
        assert!(decode(&bytes).is_ok());
        let mut bytes = vec![0x91; 10_000];
        bytes.push(0xc0);
        assert_eq!(decode(&bytes), Err(BcifError::TooDeep { limit: MAX_DEPTH }));
    }

    #[test]
    fn counts_container_elements_against_the_limit() {
        let bytes = Pack::Array(vec![Pack::Nil; 10]).encode();
        let mut budget = Budget::new(BcifLimits {
            max_elements: 9,
            ..BcifLimits::default()
        });
        assert_eq!(
            parse(&bytes, &mut budget),
            Err(BcifError::LimitExceeded { limit: 9 })
        );
    }
}
