/// The form a value takes in CIF text.
///
/// A parsed value has the kind it was read as, and is written back in the
/// same form.
///
/// # Examples
///
/// ```
/// use mmcif_cif::{Document, ValueKind};
///
/// let document = Document::parse("data_x\nloop_\n_a.v\nALA 'two words' \"O5'\" . ?\n")?;
/// let category = document.block(0).unwrap().category("a").unwrap();
/// let kinds: Vec<_> = category.column("v").unwrap().iter().map(|v| v.kind()).collect();
///
/// assert_eq!(
///     kinds,
///     [
///         ValueKind::Bare,
///         ValueKind::SingleQuoted,
///         ValueKind::DoubleQuoted,
///         ValueKind::Inapplicable,
///         ValueKind::Unknown,
///     ]
/// );
/// # Ok::<(), mmcif_cif::ParseError>(())
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ValueKind {
    /// An unquoted value, such as `ALA`.
    Bare,
    /// A value in single quotes, such as `'two words'`.
    SingleQuoted,
    /// A value in double quotes, such as `"O5'"`.
    DoubleQuoted,
    /// A text field, from a `;` at the start of a line to the next line starting with `;`.
    TextField,
    /// The inapplicable null, written `.`.
    Inapplicable,
    /// The unknown null, written `?`.
    Unknown,
}

impl ValueKind {
    /// Returns the kind [`Value::text`] gives to `text`.
    ///
    /// Text with a line break is a [`TextField`](ValueKind::TextField).
    /// Other text is [`Bare`](ValueKind::Bare) unless it is empty, is `.` or
    /// `?`, contains whitespace, starts with one of `_ # ' " $ [ ] ;`, or is a
    /// reserved word (`loop_`, `global_`, `stop_`, or a word starting with
    /// `data_` or `save_`, in any ASCII case). Text that is not bare is
    /// [`SingleQuoted`](ValueKind::SingleQuoted) if no `'` in it is followed
    /// by whitespace, otherwise [`DoubleQuoted`](ValueKind::DoubleQuoted) if
    /// no `"` in it is followed by whitespace, and otherwise a
    /// [`TextField`](ValueKind::TextField). The result is never a null kind.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_cif::ValueKind;
    ///
    /// assert_eq!(ValueKind::for_text("ALA"), ValueKind::Bare);
    /// assert_eq!(ValueKind::for_text("two words"), ValueKind::SingleQuoted);
    /// assert_eq!(ValueKind::for_text("?"), ValueKind::SingleQuoted);
    /// assert_eq!(ValueKind::for_text("it' s"), ValueKind::DoubleQuoted);
    /// assert_eq!(ValueKind::for_text("line\nbreak"), ValueKind::TextField);
    /// ```
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

    /// Returns `true` if `text` can be written as a value of this kind.
    ///
    /// [`Bare`](ValueKind::Bare) holds the text [`ValueKind::for_text`]
    /// writes bare, and also text starting with `$`, `[`, `]` or `;`. The
    /// quoted kinds hold text with no line break in which their quote is never
    /// followed by whitespace. [`TextField`](ValueKind::TextField) holds text
    /// with no line that starts with `;`. [`Inapplicable`](ValueKind::Inapplicable)
    /// holds only `.` and [`Unknown`](ValueKind::Unknown) only `?`.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_cif::ValueKind;
    ///
    /// assert!(ValueKind::Bare.can_hold("ALA"));
    /// assert!(!ValueKind::Bare.can_hold("two words"));
    /// assert!(!ValueKind::SingleQuoted.can_hold("it' s"));
    /// assert!(!ValueKind::TextField.can_hold("a\n;b"));
    /// assert!(ValueKind::Unknown.can_hold("?"));
    /// ```
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

    /// Returns `true` for [`Inapplicable`](ValueKind::Inapplicable) and
    /// [`Unknown`](ValueKind::Unknown).
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

/// A CIF value: its text and the [`ValueKind`] it is written as.
///
/// The text excludes quotes and text field delimiters. Two values are equal
/// when both their text and their kind are equal.
///
/// # Examples
///
/// ```
/// use mmcif_cif::{Document, ValueKind};
///
/// let document = Document::parse("data_x\n_struct.title 'two words'\n_struct.type ?\n")?;
/// let category = document.block(0).unwrap().category("struct").unwrap();
///
/// let title = category.column("title").unwrap().get(0).unwrap();
/// assert_eq!(title.as_str(), Some("two words"));
/// assert_eq!(title.kind(), ValueKind::SingleQuoted);
///
/// let kind = category.column("type").unwrap().get(0).unwrap();
/// assert!(kind.is_unknown());
/// assert_eq!(kind.as_str(), None);
/// # Ok::<(), mmcif_cif::ParseError>(())
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Value<'a> {
    kind: ValueKind,
    raw: &'a str,
}

/// An error from creating a [`Value`] whose kind cannot hold its text.
///
/// # Examples
///
/// ```
/// use mmcif_cif::{Value, ValueKind};
///
/// let error = Value::new(ValueKind::Bare, "two words").unwrap_err();
///
/// assert_eq!(error.kind, ValueKind::Bare);
/// assert_eq!(error.text, "two words");
/// ```
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
#[error("{text:?} cannot be written as {kind:?}")]
pub struct InvalidValue {
    /// The kind the value was given.
    pub kind: ValueKind,
    /// The text that `kind` cannot hold.
    pub text: String,
}

impl<'a> Value<'a> {
    /// The inapplicable null, `.`.
    pub const INAPPLICABLE: Value<'static> = Value {
        kind: ValueKind::Inapplicable,
        raw: ".",
    };
    /// The unknown null, `?`.
    pub const UNKNOWN: Value<'static> = Value {
        kind: ValueKind::Unknown,
        raw: "?",
    };

    /// Creates a value of the given kind with text `raw`.
    ///
    /// # Errors
    ///
    /// Returns an [`InvalidValue`] if `kind` cannot hold `raw`, as decided by
    /// [`ValueKind::can_hold`].
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_cif::{Value, ValueKind};
    ///
    /// let atom = Value::new(ValueKind::DoubleQuoted, "O5'")?;
    /// assert_eq!(atom.as_str(), Some("O5'"));
    ///
    /// assert!(Value::new(ValueKind::Unknown, ".").is_err());
    /// # Ok::<(), mmcif_cif::InvalidValue>(())
    /// ```
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

    /// Creates a value holding `text`, with the kind chosen by
    /// [`ValueKind::for_text`].
    ///
    /// The value is never null: `Value::text("?")` is the quoted text `'?'`.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_cif::{Value, ValueKind};
    ///
    /// assert_eq!(Value::text("ALA").kind(), ValueKind::Bare);
    /// assert_eq!(Value::text("two words").kind(), ValueKind::SingleQuoted);
    /// assert!(!Value::text("?").is_null());
    /// ```
    pub fn text(text: &'a str) -> Self {
        Value {
            kind: ValueKind::for_text(text),
            raw: text,
        }
    }

    pub(crate) fn from_parts(kind: ValueKind, raw: &'a str) -> Self {
        Value { kind, raw }
    }

    /// Returns the kind of the value.
    pub fn kind(&self) -> ValueKind {
        self.kind
    }

    /// Returns the text of the value, which is `.` or `?` for a null.
    ///
    /// Quotes are not included. For a text field, the text runs from after the
    /// opening `;` to the line break before the closing `;`, without that line
    /// break.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_cif::Document;
    ///
    /// let document = Document::parse("data_x\n_a.note\n;first\nsecond\n;\n_a.flag .\n")?;
    /// let category = document.block(0).unwrap().category("a").unwrap();
    ///
    /// assert_eq!(category.column("note").unwrap().get(0).unwrap().raw(), "first\nsecond");
    /// assert_eq!(category.column("flag").unwrap().get(0).unwrap().raw(), ".");
    /// # Ok::<(), mmcif_cif::ParseError>(())
    /// ```
    pub fn raw(&self) -> &'a str {
        self.raw
    }

    /// Returns the text of the value, or `None` if the value is null.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_cif::Value;
    ///
    /// assert_eq!(Value::text("ALA").as_str(), Some("ALA"));
    /// assert_eq!(Value::UNKNOWN.as_str(), None);
    /// ```
    pub fn as_str(&self) -> Option<&'a str> {
        (!self.kind.is_null()).then_some(self.raw)
    }

    /// Returns `true` if the value is the inapplicable null `.` or the unknown
    /// null `?`.
    ///
    /// The quoted values `'.'` and `'?'` are not null.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_cif::{Value, ValueKind};
    ///
    /// assert!(Value::INAPPLICABLE.is_null());
    /// assert!(!Value::new(ValueKind::SingleQuoted, "?")?.is_null());
    /// # Ok::<(), mmcif_cif::InvalidValue>(())
    /// ```
    pub fn is_null(&self) -> bool {
        self.kind.is_null()
    }

    /// Returns `true` if the value is the unknown null `?`.
    pub fn is_unknown(&self) -> bool {
        self.kind == ValueKind::Unknown
    }

    /// Returns `true` if the value is the inapplicable null `.`.
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
