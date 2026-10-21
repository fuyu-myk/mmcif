use crate::error::{ParseError, ParseErrorKind};
use crate::value::ValueKind;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum TokenKind {
    Data,
    Save,
    Global,
    Stop,
    Loop,
    Tag,
    Value(ValueKind),
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct Token {
    pub kind: TokenKind,
    pub start: usize,
    pub end: usize,
}

pub(crate) struct Tokenizer<'a> {
    src: &'a [u8],
    pos: usize,
}

const UTF8_BOM: &[u8] = b"\xEF\xBB\xBF";

fn is_whitespace(byte: u8) -> bool {
    matches!(byte, b' ' | b'\t' | b'\n' | b'\r')
}

impl<'a> Tokenizer<'a> {
    pub fn new(src: &'a [u8]) -> Self {
        let pos = if src.starts_with(UTF8_BOM) {
            UTF8_BOM.len()
        } else {
            0
        };
        Tokenizer { src, pos }
    }

    fn error(&self, offset: usize, kind: ParseErrorKind) -> ParseError {
        ParseError::at(self.src, offset, kind)
    }

    fn skip_whitespace_and_comments(&mut self) {
        while let Some(&byte) = self.src.get(self.pos) {
            if is_whitespace(byte) {
                self.pos += 1;
            } else if byte == b'#' {
                self.pos = memchr::memchr(b'\n', &self.src[self.pos..])
                    .map_or(self.src.len(), |i| self.pos + i + 1);
            } else {
                break;
            }
        }
    }

    fn at_line_start(&self) -> bool {
        self.pos == 0 || self.src[self.pos - 1] == b'\n'
    }

    pub fn next_token(&mut self) -> Result<Option<Token>, ParseError> {
        self.skip_whitespace_and_comments();
        let Some(&first) = self.src.get(self.pos) else {
            return Ok(None);
        };
        let token = match first {
            b'\'' | b'"' => self.quoted(first)?,
            b';' if self.at_line_start() => self.text_field()?,
            _ => self.bare(),
        };
        Ok(Some(token))
    }

    fn quoted(&mut self, quote: u8) -> Result<Token, ParseError> {
        let open = self.pos;
        let mut search = open + 1;
        loop {
            let found = memchr::memchr2(quote, b'\n', &self.src[search..]).map(|i| search + i);
            match found {
                Some(i) if self.src[i] == quote => {
                    let closes = self.src.get(i + 1).is_none_or(|&b| is_whitespace(b));
                    if closes {
                        self.pos = i + 1;
                        let kind = if quote == b'\'' {
                            ValueKind::SingleQuoted
                        } else {
                            ValueKind::DoubleQuoted
                        };
                        return Ok(Token {
                            kind: TokenKind::Value(kind),
                            start: open + 1,
                            end: i,
                        });
                    }
                    search = i + 1;
                }
                _ => return Err(self.error(open, ParseErrorKind::UnterminatedQuote)),
            }
        }
    }

    fn text_field(&mut self) -> Result<Token, ParseError> {
        let open = self.pos;
        let start = open + 1;
        let newline = memchr::memmem::find(&self.src[start..], b"\n;")
            .map(|i| start + i)
            .ok_or_else(|| self.error(open, ParseErrorKind::UnterminatedTextField))?;
        let end = if newline > start && self.src[newline - 1] == b'\r' {
            newline - 1
        } else {
            newline
        };
        self.pos = newline + 2;
        Ok(Token {
            kind: TokenKind::Value(ValueKind::TextField),
            start,
            end,
        })
    }

    fn bare(&mut self) -> Token {
        let start = self.pos;
        let end = self.src[start..]
            .iter()
            .position(|&b| is_whitespace(b))
            .map_or(self.src.len(), |i| start + i);
        self.pos = end;
        let word = &self.src[start..end];
        let has_prefix = |prefix: &[u8]| {
            word.len() >= prefix.len() && word[..prefix.len()].eq_ignore_ascii_case(prefix)
        };
        let (kind, start) = if word[0] == b'_' {
            (TokenKind::Tag, start)
        } else if has_prefix(b"data_") {
            (TokenKind::Data, start + 5)
        } else if has_prefix(b"save_") {
            (TokenKind::Save, start + 5)
        } else if word.eq_ignore_ascii_case(b"loop_") {
            (TokenKind::Loop, start)
        } else if word.eq_ignore_ascii_case(b"global_") {
            (TokenKind::Global, start)
        } else if word.eq_ignore_ascii_case(b"stop_") {
            (TokenKind::Stop, start)
        } else if word == b"." {
            (TokenKind::Value(ValueKind::Inapplicable), start)
        } else if word == b"?" {
            (TokenKind::Value(ValueKind::Unknown), start)
        } else {
            (TokenKind::Value(ValueKind::Bare), start)
        };
        Token { kind, start, end }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tokens(src: &str) -> Vec<(TokenKind, &str)> {
        let mut tokenizer = Tokenizer::new(src.as_bytes());
        let mut out = Vec::new();
        while let Some(token) = tokenizer.next_token().unwrap() {
            out.push((token.kind, &src[token.start..token.end]));
        }
        out
    }

    fn error_kind(src: &str) -> ParseErrorKind {
        let mut tokenizer = Tokenizer::new(src.as_bytes());
        loop {
            match tokenizer.next_token() {
                Ok(Some(_)) => continue,
                Ok(None) => panic!("expected an error for {src:?}"),
                Err(error) => return error.kind,
            }
        }
    }

    use TokenKind::*;
    use ValueKind::*;

    #[test]
    fn reads_headers_tags_and_values() {
        assert_eq!(
            tokens("data_1ABC\n_entry.id 1ABC # comment\nloop_ _a.x . ?"),
            vec![
                (Data, "1ABC"),
                (Tag, "_entry.id"),
                (Value(Bare), "1ABC"),
                (Loop, "loop_"),
                (Tag, "_a.x"),
                (Value(Inapplicable), "."),
                (Value(Unknown), "?"),
            ]
        );
    }

    #[test]
    fn quote_closes_only_before_whitespace() {
        assert_eq!(
            tokens("'a'b' \"O5'\" 'end'"),
            vec![
                (Value(SingleQuoted), "a'b"),
                (Value(DoubleQuoted), "O5'"),
                (Value(SingleQuoted), "end"),
            ]
        );
    }

    #[test]
    fn text_field_spans_lines_and_drops_final_newline() {
        assert_eq!(
            tokens(";first\nsecond\n;\nnext"),
            vec![(Value(TextField), "first\nsecond"), (Value(Bare), "next")]
        );
        assert_eq!(tokens(";\n;"), vec![(Value(TextField), "")]);
        assert_eq!(tokens(";a\r\nb\r\n;"), vec![(Value(TextField), "a\r\nb")]);
    }

    #[test]
    fn semicolon_mid_line_is_bare() {
        assert_eq!(
            tokens("x ;y"),
            vec![(Value(Bare), "x"), (Value(Bare), ";y")]
        );
    }

    #[test]
    fn skips_a_leading_byte_order_mark() {
        assert_eq!(tokens("\u{FEFF}data_x"), vec![(Data, "x")]);
    }

    #[test]
    fn hash_inside_token_is_not_a_comment() {
        assert_eq!(tokens("a#b"), vec![(Value(Bare), "a#b")]);
    }

    #[test]
    fn reserved_words_are_case_insensitive() {
        assert_eq!(
            tokens("DATA_x LOOP_ Save_f save_ global_ STOP_"),
            vec![
                (Data, "x"),
                (Loop, "LOOP_"),
                (Save, "f"),
                (Save, ""),
                (Global, "global_"),
                (Stop, "STOP_"),
            ]
        );
    }

    #[test]
    fn reports_unterminated_values() {
        assert_eq!(error_kind("'abc"), ParseErrorKind::UnterminatedQuote);
        assert_eq!(error_kind("'abc\n'"), ParseErrorKind::UnterminatedQuote);
        assert_eq!(error_kind(";abc\n"), ParseErrorKind::UnterminatedTextField);
    }

    #[test]
    fn error_positions_are_one_based() {
        let mut tokenizer = Tokenizer::new(b"data_x\n  'oops");
        tokenizer.next_token().unwrap();
        let error = tokenizer.next_token().unwrap_err();
        assert_eq!((error.line, error.column, error.offset), (2, 3, 9));
    }
}
