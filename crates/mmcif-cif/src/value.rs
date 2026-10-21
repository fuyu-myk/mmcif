#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ValueKind {
    Bare,
    SingleQuoted,
    DoubleQuoted,
    TextField,
    Inapplicable,
    Unknown,
}

impl ValueKind {
    pub fn for_text(text: &str) -> ValueKind {
        if text.contains('\n') {
            ValueKind::TextField
        } else if is_bare_safe(text) {
            ValueKind::Bare
        } else if is_quotable(text, '\'') {
            ValueKind::SingleQuoted
        } else if is_quotable(text, '"') {
            ValueKind::DoubleQuoted
        } else {
            ValueKind::TextField
        }
    }

    pub fn can_hold(self, text: &str) -> bool {
        match self {
            ValueKind::Bare => is_bare_writable(text),
            ValueKind::SingleQuoted => is_quotable(text, '\''),
            ValueKind::DoubleQuoted => is_quotable(text, '"'),
            ValueKind::TextField => !text.contains("\n;"),
            ValueKind::Inapplicable => text == ".",
            ValueKind::Unknown => text == "?",
        }
    }

    pub fn is_null(self) -> bool {
        matches!(self, ValueKind::Inapplicable | ValueKind::Unknown)
    }
}

fn is_cif_whitespace(byte: u8) -> bool {
    matches!(byte, b' ' | b'\t' | b'\n' | b'\r')
}

pub(crate) fn is_reserved_word(text: &str) -> bool {
    let bytes = text.as_bytes();
    let starts_with = |prefix: &[u8]| {
        bytes.len() >= prefix.len() && bytes[..prefix.len()].eq_ignore_ascii_case(prefix)
    };
    starts_with(b"data_")
        || starts_with(b"save_")
        || text.eq_ignore_ascii_case("loop_")
        || text.eq_ignore_ascii_case("global_")
        || text.eq_ignore_ascii_case("stop_")
}

fn is_bare_body(text: &str) -> bool {
    !text.bytes().any(is_cif_whitespace) && !is_reserved_word(text)
}

fn is_bare_writable(text: &str) -> bool {
    match text.as_bytes().first() {
        None | Some(b'_' | b'#' | b'\'' | b'"') => false,
        Some(_) => text != "." && text != "?" && is_bare_body(text),
    }
}

fn is_bare_safe(text: &str) -> bool {
    !text.starts_with(['$', '[', ']', ';']) && is_bare_writable(text)
}

fn is_quotable(text: &str, quote: char) -> bool {
    if text.contains('\n') {
        return false;
    }
    let bytes = text.as_bytes();
    !bytes
        .windows(2)
        .any(|pair| pair[0] == quote as u8 && is_cif_whitespace(pair[1]))
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Value<'a> {
    kind: ValueKind,
    raw: &'a str,
}

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
#[error("{text:?} cannot be written as {kind:?}")]
pub struct InvalidValue {
    pub kind: ValueKind,
    pub text: String,
}

impl<'a> Value<'a> {
    pub const INAPPLICABLE: Value<'static> = Value {
        kind: ValueKind::Inapplicable,
        raw: ".",
    };
    pub const UNKNOWN: Value<'static> = Value {
        kind: ValueKind::Unknown,
        raw: "?",
    };

    pub fn new(kind: ValueKind, raw: &'a str) -> Result<Self, InvalidValue> {
        if kind.can_hold(raw) {
            Ok(Value { kind, raw })
        } else {
            Err(InvalidValue {
                kind,
                text: raw.to_owned(),
            })
        }
    }

    pub fn text(text: &'a str) -> Self {
        Value {
            kind: ValueKind::for_text(text),
            raw: text,
        }
    }

    pub(crate) fn from_parts(kind: ValueKind, raw: &'a str) -> Self {
        Value { kind, raw }
    }

    pub fn kind(&self) -> ValueKind {
        self.kind
    }

    pub fn raw(&self) -> &'a str {
        self.raw
    }

    pub fn as_str(&self) -> Option<&'a str> {
        (!self.kind.is_null()).then_some(self.raw)
    }

    pub fn is_null(&self) -> bool {
        self.kind.is_null()
    }

    pub fn is_unknown(&self) -> bool {
        self.kind == ValueKind::Unknown
    }

    pub fn is_inapplicable(&self) -> bool {
        self.kind == ValueKind::Inapplicable
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn chooses_bare_for_plain_tokens() {
        assert_eq!(ValueKind::for_text("ALA"), ValueKind::Bare);
        assert_eq!(ValueKind::for_text("1.234"), ValueKind::Bare);
        assert_eq!(ValueKind::for_text("O5'"), ValueKind::Bare);
    }

    #[test]
    fn quotes_values_that_would_be_misread() {
        assert_eq!(ValueKind::for_text(""), ValueKind::SingleQuoted);
        assert_eq!(ValueKind::for_text("."), ValueKind::SingleQuoted);
        assert_eq!(ValueKind::for_text("?"), ValueKind::SingleQuoted);
        assert_eq!(ValueKind::for_text("_x"), ValueKind::SingleQuoted);
        assert_eq!(ValueKind::for_text("loop_"), ValueKind::SingleQuoted);
        assert_eq!(ValueKind::for_text("DATA_x"), ValueKind::SingleQuoted);
        assert_eq!(ValueKind::for_text("two words"), ValueKind::SingleQuoted);
        assert_eq!(ValueKind::for_text("it' s"), ValueKind::DoubleQuoted);
        assert_eq!(ValueKind::for_text("a' b\" c"), ValueKind::TextField);
        assert_eq!(ValueKind::for_text("line\nbreak"), ValueKind::TextField);
    }

    #[test]
    fn rejects_kinds_that_cannot_round_trip() {
        assert!(Value::new(ValueKind::Bare, "two words").is_err());
        assert!(Value::new(ValueKind::SingleQuoted, "a' b").is_err());
        assert!(Value::new(ValueKind::TextField, "a\n;b").is_err());
        assert!(Value::new(ValueKind::Unknown, ".").is_err());
        assert!(Value::new(ValueKind::Bare, ";mid").is_ok());
        assert_eq!(ValueKind::for_text(";mid"), ValueKind::SingleQuoted);
    }

    #[test]
    fn null_values_have_no_text() {
        assert_eq!(Value::UNKNOWN.as_str(), None);
        assert_eq!(Value::INAPPLICABLE.as_str(), None);
        assert_eq!(Value::text("x").as_str(), Some("x"));
    }
}
