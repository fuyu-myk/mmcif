use super::BcifError;
use super::msgpack::MsgValue;

impl<'a> MsgValue<'a> {
    fn as_text(&self) -> Option<Result<&'a str, BcifError>> {
        match *self {
            MsgValue::Str(text) => Some(Ok(text)),
            MsgValue::Bin(bytes) => {
                Some(std::str::from_utf8(bytes).map_err(|_| BcifError::InvalidUtf8))
            }
            _ => None,
        }
    }

    fn get(&self, key: &str) -> Option<&MsgValue<'a>> {
        let MsgValue::Map(entries) = self else {
            return None;
        };
        entries
            .iter()
            .find(|(name, _)| matches!(name.as_text(), Some(Ok(name)) if name == key))
            .map(|(_, value)| value)
    }
}

fn unexpected(field: &'static str, expected: &'static str) -> BcifError {
    BcifError::UnexpectedType { field, expected }
}

#[derive(Debug, Clone, Copy)]
pub(crate) struct Record<'v> {
    value: &'v MsgValue<'v>,
}

impl<'v> Record<'v> {
    pub fn new(value: &'v MsgValue<'v>, what: &'static str) -> Result<Self, BcifError> {
        match value {
            MsgValue::Map(_) => Ok(Record { value }),
            _ => Err(unexpected(what, "a map")),
        }
    }

    pub fn field(&self, key: &'static str) -> Result<&'v MsgValue<'v>, BcifError> {
        self.value.get(key).ok_or(BcifError::MissingField(key))
    }

    pub fn optional(&self, key: &'static str) -> Option<&'v MsgValue<'v>> {
        self.value
            .get(key)
            .filter(|value| !matches!(value, MsgValue::Nil))
    }

    pub fn record(&self, key: &'static str) -> Result<Record<'v>, BcifError> {
        Record::new(self.field(key)?, key)
    }

    pub fn str(&self, key: &'static str) -> Result<&'v str, BcifError> {
        self.field(key)?
            .as_text()
            .unwrap_or(Err(unexpected(key, "a string")))
    }

    pub fn bin(&self, key: &'static str) -> Result<&'v [u8], BcifError> {
        match *self.field(key)? {
            MsgValue::Bin(bytes) => Ok(bytes),
            _ => Err(unexpected(key, "binary data")),
        }
    }

    pub fn int(&self, key: &'static str) -> Result<i64, BcifError> {
        match *self.field(key)? {
            MsgValue::Integer(value) => {
                i64::try_from(value).map_err(|_| unexpected(key, "a 64-bit integer"))
            }
            _ => Err(unexpected(key, "an integer")),
        }
    }

    pub fn count(&self, key: &'static str) -> Result<usize, BcifError> {
        usize::try_from(self.int(key)?).map_err(|_| unexpected(key, "a non-negative integer"))
    }

    pub fn number(&self, key: &'static str) -> Result<f64, BcifError> {
        match *self.field(key)? {
            MsgValue::Integer(value) => Ok(value as f64),
            MsgValue::Float(value) => Ok(value),
            _ => Err(unexpected(key, "a number")),
        }
    }

    pub fn bool(&self, key: &'static str) -> Result<bool, BcifError> {
        match *self.field(key)? {
            MsgValue::Bool(value) => Ok(value),
            _ => Err(unexpected(key, "a boolean")),
        }
    }

    pub fn list(&self, key: &'static str) -> Result<&'v [MsgValue<'v>], BcifError> {
        match self.field(key)? {
            MsgValue::Array(items) => Ok(items),
            _ => Err(unexpected(key, "an array")),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bcif::msgpack::parse;
    use crate::bcif::test_pack::{Pack, obj};
    use crate::bcif::{BcifLimits, Budget};

    fn decode(bytes: &[u8]) -> Result<MsgValue<'_>, BcifError> {
        parse(bytes, &mut Budget::new(BcifLimits::default()))
    }

    #[test]
    fn round_trips_the_test_encoder() {
        let pack = obj([
            ("n", Pack::Int(-70_000)),
            ("f", Pack::Float(0.25)),
            ("s", Pack::Str("x".repeat(300))),
            ("b", Pack::Bin(vec![1; 70_000])),
            ("a", Pack::Array(vec![Pack::Nil; 20])),
        ]);
        let bytes = pack.encode();
        let value = decode(&bytes).unwrap();
        let record = Record::new(&value, "root").unwrap();
        assert_eq!(record.int("n").unwrap(), -70_000);
        assert_eq!(record.number("f").unwrap(), 0.25);
        assert_eq!(record.str("s").unwrap().len(), 300);
        assert_eq!(record.bin("b").unwrap().len(), 70_000);
        assert_eq!(record.list("a").unwrap().len(), 20);
    }

    #[test]
    fn reports_missing_and_mistyped_fields() {
        let bytes = obj([("k", Pack::Int(1)), ("big", Pack::UInt(u64::MAX))]).encode();
        let value = decode(&bytes).unwrap();
        let record = Record::new(&value, "root").unwrap();
        assert_eq!(record.str("x"), Err(BcifError::MissingField("x")));
        assert!(matches!(
            record.str("k"),
            Err(BcifError::UnexpectedType { field: "k", .. })
        ));
        assert!(matches!(
            record.int("big"),
            Err(BcifError::UnexpectedType { field: "big", .. })
        ));
        assert!(record.optional("x").is_none());
        assert!(Record::new(&MsgValue::Nil, "root").is_err());
    }

    #[test]
    fn accepts_binary_keys_and_strings() {
        let bytes = [
            0x82, 0xc4, 0x01, b'k', 0xc4, 0x01, b'v', 0xa1, b'x', 0xc4, 0x01, 0xff,
        ];
        let value = decode(&bytes).unwrap();
        let record = Record::new(&value, "root").unwrap();
        assert_eq!(record.str("k"), Ok("v"));
        assert_eq!(record.str("x"), Err(BcifError::InvalidUtf8));
    }
}
