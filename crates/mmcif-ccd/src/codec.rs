use crate::error::{CcdError, corrupt};

const TRUNCATED: &str = "data ends early";

pub(crate) fn put_varint(out: &mut Vec<u8>, mut value: u64) {
    while value >= 0x80 {
        out.push((value as u8) | 0x80);
        value >>= 7;
    }
    out.push(value as u8);
}

pub(crate) fn varint_len(value: u64) -> usize {
    let bits = (u64::BITS - value.leading_zeros()).max(1) as usize;
    bits.div_ceil(7)
}

pub(crate) fn put_str(out: &mut Vec<u8>, text: &str) {
    put_varint(out, text.len() as u64);
    out.extend_from_slice(text.as_bytes());
}

pub(crate) fn put_opt_str(out: &mut Vec<u8>, text: Option<&str>) {
    out.push(u8::from(text.is_some()));
    if let Some(text) = text {
        put_str(out, text);
    }
}

pub(crate) struct Reader<'a> {
    bytes: &'a [u8],
}

impl<'a> Reader<'a> {
    pub fn new(bytes: &'a [u8]) -> Self {
        Reader { bytes }
    }

    pub fn remaining(&self) -> usize {
        self.bytes.len()
    }

    pub fn take(&mut self, len: usize) -> Result<&'a [u8], CcdError> {
        let (taken, rest) = self.bytes.split_at_checked(len).ok_or(corrupt(TRUNCATED))?;
        self.bytes = rest;
        Ok(taken)
    }

    pub fn array<const N: usize>(&mut self) -> Result<[u8; N], CcdError> {
        let mut array = [0; N];
        array.copy_from_slice(self.take(N)?);
        Ok(array)
    }

    pub fn u8(&mut self) -> Result<u8, CcdError> {
        Ok(self.array::<1>()?[0])
    }

    pub fn i8(&mut self) -> Result<i8, CcdError> {
        Ok(i8::from_le_bytes(self.array()?))
    }

    pub fn u16(&mut self) -> Result<u16, CcdError> {
        Ok(u16::from_le_bytes(self.array()?))
    }

    pub fn u32(&mut self) -> Result<u32, CcdError> {
        Ok(u32::from_le_bytes(self.array()?))
    }

    pub fn i32(&mut self) -> Result<i32, CcdError> {
        Ok(i32::from_le_bytes(self.array()?))
    }

    pub fn u64(&mut self) -> Result<u64, CcdError> {
        Ok(u64::from_le_bytes(self.array()?))
    }

    pub fn bool(&mut self) -> Result<bool, CcdError> {
        match self.u8()? {
            0 => Ok(false),
            1 => Ok(true),
            _ => Err(corrupt("invalid presence byte")),
        }
    }

    pub fn varint(&mut self) -> Result<u64, CcdError> {
        let mut value = 0u64;
        for shift in (0..u64::BITS).step_by(7) {
            let byte = self.u8()?;
            let part = u64::from(byte & 0x7F);
            if (part << shift) >> shift != part {
                return Err(corrupt("varint overflows 64 bits"));
            }
            value |= part << shift;
            if byte & 0x80 == 0 {
                return Ok(value);
            }
        }
        Err(corrupt("varint overflows 64 bits"))
    }

    pub fn length(&mut self) -> Result<usize, CcdError> {
        usize::try_from(self.varint()?).map_err(|_| corrupt(TRUNCATED))
    }

    pub fn count(&mut self, min_item_len: usize) -> Result<usize, CcdError> {
        let count = self.length()?;
        if count > self.remaining() / min_item_len {
            return Err(corrupt("count exceeds remaining data"));
        }
        Ok(count)
    }

    pub fn str(&mut self) -> Result<&'a str, CcdError> {
        let len = self.length()?;
        std::str::from_utf8(self.take(len)?).map_err(|_| corrupt("text is not UTF-8"))
    }

    pub fn opt_str(&mut self) -> Result<Option<&'a str>, CcdError> {
        if self.bool()? {
            Ok(Some(self.str()?))
        } else {
            Ok(None)
        }
    }

    pub fn finish(self) -> Result<(), CcdError> {
        if self.bytes.is_empty() {
            Ok(())
        } else {
            Err(corrupt("unexpected trailing bytes"))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn encoded(value: u64) -> Vec<u8> {
        let mut out = Vec::new();
        put_varint(&mut out, value);
        out
    }

    #[test]
    fn varints_round_trip() {
        let values = [
            0,
            1,
            127,
            128,
            300,
            16_383,
            16_384,
            u32::MAX.into(),
            u64::MAX,
        ];
        for value in values {
            let bytes = encoded(value);
            assert_eq!(bytes.len(), varint_len(value), "{value}");
            let mut reader = Reader::new(&bytes);
            assert_eq!(reader.varint().unwrap(), value);
            reader.finish().unwrap();
        }
    }

    #[test]
    fn rejects_overlong_varints() {
        let mut too_wide = vec![0xFF; 9];
        too_wide.push(0x02);
        assert!(Reader::new(&too_wide).varint().is_err());
        assert!(Reader::new(&[0x80; 11]).varint().is_err());
        assert!(Reader::new(&[0x80, 0x80]).varint().is_err());
    }

    #[test]
    fn rejects_reads_past_the_end() {
        let mut reader = Reader::new(&[1, 2, 3]);
        assert!(reader.u32().is_err());
        assert_eq!(reader.u16().unwrap(), 0x0201);
        assert!(reader.take(2).is_err());
        assert!(matches!(
            reader.finish(),
            Err(CcdError::CorruptIndex { .. })
        ));
    }

    #[test]
    fn rejects_counts_larger_than_the_data() {
        let mut bytes = encoded(u64::MAX);
        bytes.extend_from_slice(&[0; 4]);
        assert!(Reader::new(&bytes).count(1).is_err());
        let mut bytes = encoded(3);
        bytes.extend_from_slice(&[0; 5]);
        assert!(Reader::new(&bytes).count(2).is_err());
        assert_eq!(Reader::new(&bytes).count(1).unwrap(), 3);
    }

    #[test]
    fn strings_round_trip() {
        let mut out = Vec::new();
        put_str(&mut out, "O5'");
        put_opt_str(&mut out, None);
        put_opt_str(&mut out, Some(""));
        let mut reader = Reader::new(&out);
        assert_eq!(reader.str().unwrap(), "O5'");
        assert_eq!(reader.opt_str().unwrap(), None);
        assert_eq!(reader.opt_str().unwrap(), Some(""));
        reader.finish().unwrap();
    }

    #[test]
    fn rejects_invalid_text() {
        assert!(Reader::new(&[2, 0xC3, 0x28]).str().is_err());
        assert!(Reader::new(&[2]).opt_str().is_err());
    }
}
