use std::borrow::Borrow;
use std::fmt;

macro_rules! string_id {
    ($name:ident) => {
        #[derive(Debug, Clone, PartialEq, Eq, Hash, PartialOrd, Ord)]
        pub struct $name(Box<str>);

        impl $name {
            pub fn new(value: impl Into<Box<str>>) -> Self {
                $name(value.into())
            }

            pub fn as_str(&self) -> &str {
                &self.0
            }
        }

        impl fmt::Display for $name {
            fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                f.write_str(&self.0)
            }
        }

        impl From<&str> for $name {
            fn from(value: &str) -> Self {
                $name::new(value)
            }
        }

        impl Borrow<str> for $name {
            fn borrow(&self) -> &str {
                &self.0
            }
        }
    };
}

string_id!(LabelAsymId);
string_id!(AuthAsymId);
string_id!(EntityId);
string_id!(CompId);
string_id!(AtomName);
string_id!(OperatorId);
string_id!(AssemblyId);

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct LabelSeqId(u32);

impl LabelSeqId {
    pub fn new(value: u32) -> Option<Self> {
        (value > 0).then_some(LabelSeqId(value))
    }

    pub fn get(self) -> u32 {
        self.0
    }
}

impl fmt::Display for LabelSeqId {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.0.fmt(f)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct AuthSeqId(i32);

impl AuthSeqId {
    pub fn new(value: i32) -> Self {
        AuthSeqId(value)
    }

    pub fn get(self) -> i32 {
        self.0
    }
}

impl fmt::Display for AuthSeqId {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.0.fmt(f)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Default)]
pub struct InsCode(Option<char>);

impl InsCode {
    pub const NONE: InsCode = InsCode(None);

    pub fn new(code: char) -> Self {
        InsCode(Some(code))
    }

    pub fn get(self) -> Option<char> {
        self.0
    }
}

impl fmt::Display for InsCode {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self.0 {
            Some(code) => write!(f, "{code}"),
            None => Ok(()),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Default)]
pub struct AltLoc(Option<char>);

impl AltLoc {
    pub const NONE: AltLoc = AltLoc(None);

    pub fn new(code: char) -> Self {
        AltLoc(Some(code))
    }

    pub fn get(self) -> Option<char> {
        self.0
    }

    pub fn is_none(self) -> bool {
        self.0.is_none()
    }
}

impl fmt::Display for AltLoc {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self.0 {
            Some(code) => write!(f, "{code}"),
            None => Ok(()),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct ModelNumber(i32);

impl ModelNumber {
    pub fn new(value: i32) -> Self {
        ModelNumber(value)
    }

    pub fn get(self) -> i32 {
        self.0
    }
}

impl fmt::Display for ModelNumber {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.0.fmt(f)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct Element([u8; 2]);

impl Element {
    pub fn parse(symbol: &str) -> Option<Element> {
        match symbol.as_bytes() {
            [first] if first.is_ascii_alphabetic() => {
                Some(Element([first.to_ascii_uppercase(), 0]))
            }
            [first, second] if first.is_ascii_alphabetic() && second.is_ascii_alphabetic() => {
                Some(Element([
                    first.to_ascii_uppercase(),
                    second.to_ascii_lowercase(),
                ]))
            }
            _ => None,
        }
    }

    pub fn symbol(&self) -> &str {
        let len = if self.0[1] == 0 { 1 } else { 2 };
        std::str::from_utf8(&self.0[..len]).expect("element symbols are ASCII letters")
    }

    pub fn is_hydrogen(&self) -> bool {
        matches!(self.symbol(), "H" | "D" | "T")
    }
}

impl fmt::Display for Element {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.symbol())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn label_seq_ids_start_at_one() {
        assert!(LabelSeqId::new(0).is_none());
        assert_eq!(LabelSeqId::new(7).unwrap().get(), 7);
    }

    #[test]
    fn elements_normalise_case() {
        assert_eq!(Element::parse("FE").unwrap().symbol(), "Fe");
        assert_eq!(Element::parse("c").unwrap().symbol(), "C");
        assert!(Element::parse("D").unwrap().is_hydrogen());
        assert!(Element::parse("").is_none());
        assert!(Element::parse("C1").is_none());
        assert!(Element::parse("ABC").is_none());
    }

    #[test]
    fn missing_codes_display_empty() {
        assert_eq!(InsCode::NONE.to_string(), "");
        assert_eq!(InsCode::new('A').to_string(), "A");
        assert_eq!(AltLoc::new('B').get(), Some('B'));
    }
}
