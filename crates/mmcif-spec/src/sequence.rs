use mmcif_model::CompId;

use crate::diagnostic::Diagnostic;
use crate::span::{SourceId, Span};

const ONE_LETTER: [(char, &str); 21] = [
    ('A', "ALA"),
    ('C', "CYS"),
    ('D', "ASP"),
    ('E', "GLU"),
    ('F', "PHE"),
    ('G', "GLY"),
    ('H', "HIS"),
    ('I', "ILE"),
    ('K', "LYS"),
    ('L', "LEU"),
    ('M', "MET"),
    ('N', "ASN"),
    ('P', "PRO"),
    ('Q', "GLN"),
    ('R', "ARG"),
    ('S', "SER"),
    ('T', "THR"),
    ('V', "VAL"),
    ('W', "TRP"),
    ('Y', "TYR"),
    ('X', "UNK"),
];

/// Returns the one-letter code of a standard amino acid component, or `None`
/// for any other component.
///
/// `UNK` maps to `X`.
///
/// # Examples
///
/// ```
/// use mmcif_model::CompId;
/// use mmcif_spec::one_letter;
///
/// assert_eq!(one_letter(&CompId::new("TRP")), Some('W'));
/// assert_eq!(one_letter(&CompId::new("UNK")), Some('X'));
/// assert_eq!(one_letter(&CompId::new("HOH")), None);
/// ```
pub fn one_letter(component: &CompId) -> Option<char> {
    ONE_LETTER
        .iter()
        .find(|(_, code)| *code == component.as_str())
        .map(|(letter, _)| *letter)
}

/// Returns the standard amino acid component for an uppercase one-letter
/// code, or `None` if `letter` is not one.
///
/// `X` maps to `UNK`.
///
/// # Examples
///
/// ```
/// use mmcif_model::CompId;
/// use mmcif_spec::component;
///
/// assert_eq!(component('C'), Some(CompId::new("CYS")));
/// assert_eq!(component('X'), Some(CompId::new("UNK")));
/// assert_eq!(component('B'), None);
/// assert_eq!(component('c'), None);
/// ```
pub fn component(letter: char) -> Option<CompId> {
    ONE_LETTER
        .iter()
        .find(|(known, _)| *known == letter)
        .map(|(_, code)| CompId::new(*code))
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Token {
    Designed { min: u32, max: u32 },
    Residue(CompId),
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Piece {
    pub span: Span,
    pub token: Token,
}

pub fn parse(source: SourceId, scalar: Span, text: &str) -> Result<Vec<Piece>, Diagnostic> {
    let mut pieces = Vec::new();
    let mut rest = text.char_indices().peekable();
    while let Some((start, letter)) = rest.next() {
        if letter.is_ascii_digit() {
            let mut end = start + letter.len_utf8();
            while let Some((offset, digit)) = rest.peek() {
                if !digit.is_ascii_digit() && *digit != '.' {
                    break;
                }
                end = offset + digit.len_utf8();
                rest.next();
            }
            let span = inner(scalar, text, start, end);
            pieces.push(Piece {
                span,
                token: run(source, span, &text[start..end])?,
            });
        } else if letter.is_ascii_uppercase() {
            let span = inner(scalar, text, start, start + letter.len_utf8());
            match component(letter) {
                Some(comp) => pieces.push(Piece {
                    span,
                    token: Token::Residue(comp),
                }),
                None => {
                    return Err(Diagnostic::error(
                        "unknown-residue-letter",
                        format!("`{letter}` is not a residue code"),
                    )
                    .at(source, span, "not one of the 20 standard residues")
                    .help("use a one-letter code, or a run length such as 12"));
                }
            }
        } else {
            return Err(Diagnostic::error(
                "bad-sequence",
                format!("`{letter}` cannot appear in a sequence"),
            )
            .at(
                source,
                inner(scalar, text, start, start + letter.len_utf8()),
                "expected a residue code or a run length",
            ));
        }
    }
    if pieces.is_empty() {
        return Err(
            Diagnostic::error("bad-sequence", "the sequence is empty").at(
                source,
                scalar,
                "a chain needs at least one residue",
            ),
        );
    }
    Ok(pieces)
}

fn run(source: SourceId, span: Span, text: &str) -> Result<Token, Diagnostic> {
    let bad = |message: &str, help: &str| {
        Diagnostic::error("bad-sequence", format!("`{text}` is not a run length"))
            .at(source, span, message.to_string())
            .help(help.to_string())
    };
    let (min, max) = match text.split_once("..") {
        None => (text, text),
        Some((min, max)) => (min, max),
    };
    if min.is_empty() || max.is_empty() {
        return Err(bad(
            "a run needs both ends",
            "write a count such as 12, or a range such as 80..140",
        ));
    }
    let number = |part: &str| {
        part.parse::<u32>()
            .map_err(|_| bad("not a number", "write a count such as 12"))
    };
    let (min, max) = (number(min)?, number(max)?);
    if min > max {
        return Err(bad(
            "the range counts down",
            "write the shorter length first",
        ));
    }
    if max == 0 {
        return Err(bad("a run of zero residues", "use a count of 1 or more"));
    }
    Ok(Token::Designed { min, max })
}

fn inner(scalar: Span, text: &str, start: usize, end: usize) -> Span {
    if scalar.end() - scalar.start() == text.len() {
        Span::new(scalar.start() + start, scalar.start() + end)
    } else {
        scalar
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::span::Sources;

    fn source() -> SourceId {
        Sources::new().add("design.yaml", "")
    }

    fn pieces(text: &str) -> Vec<Piece> {
        parse(source(), Span::new(10, 10 + text.len()), text).unwrap()
    }

    fn designed(min: u32, max: u32) -> Token {
        Token::Designed { min, max }
    }

    fn residue(code: &str) -> Token {
        Token::Residue(CompId::new(code))
    }

    #[test]
    fn a_range_becomes_one_designed_run() {
        let parsed = pieces("80..140");
        assert_eq!(parsed.len(), 1);
        assert_eq!(parsed[0].token, designed(80, 140));
        assert_eq!(parsed[0].span, Span::new(10, 17));
    }

    #[test]
    fn a_bare_count_is_an_exact_run() {
        assert_eq!(pieces("17")[0].token, designed(17, 17));
    }

    #[test]
    fn literal_residues_become_one_piece_each() {
        let parsed = pieces("AAVT");
        assert_eq!(parsed.len(), 4);
        assert_eq!(parsed[0].token, residue("ALA"));
        assert_eq!(parsed[2].token, residue("VAL"));
        assert_eq!(parsed[2].span, Span::new(12, 13));
    }

    #[test]
    fn runs_and_residues_interleave_in_order() {
        let parsed = pieces("3..5C6C3");
        let tokens: Vec<&Token> = parsed.iter().map(|piece| &piece.token).collect();
        assert_eq!(
            tokens,
            [
                &designed(3, 5),
                &residue("CYS"),
                &designed(6, 6),
                &residue("CYS"),
                &designed(3, 3)
            ]
        );
        assert_eq!(parsed[1].span, Span::new(14, 15));
        assert_eq!(parsed[4].span, Span::new(17, 18));
    }

    #[test]
    fn a_quoted_scalar_falls_back_to_the_whole_span() {
        let text = "3..5C";
        let parsed = parse(source(), Span::new(10, 17), text).unwrap();
        assert!(parsed.iter().all(|piece| piece.span == Span::new(10, 17)));
    }

    #[test]
    fn malformed_runs_are_refused_with_their_reason() {
        let fail = |text: &str| {
            parse(source(), Span::new(0, text.len()), text)
                .unwrap_err()
                .code
        };
        assert_eq!(fail("20..15"), "bad-sequence");
        assert_eq!(fail("0"), "bad-sequence");
        assert_eq!(fail("55.."), "bad-sequence");
        assert_eq!(fail("1..2..3"), "bad-sequence");
        assert_eq!(fail(""), "bad-sequence");
        assert_eq!(fail("ac"), "bad-sequence");
        assert_eq!(fail("B"), "unknown-residue-letter");
    }

    #[test]
    fn a_descending_range_says_which_way_round_to_write_it() {
        let error = parse(source(), Span::new(0, 6), "20..15").unwrap_err();
        assert_eq!(error.labels[0].message, "the range counts down");
        assert_eq!(
            error.help.as_deref(),
            Some("write the shorter length first")
        );
    }

    #[test]
    fn every_standard_letter_maps_to_a_component() {
        assert_eq!(component('W'), Some(CompId::new("TRP")));
        assert_eq!(component('X'), Some(CompId::new("UNK")));
        assert_eq!(component('B'), None);
        assert_eq!(ONE_LETTER.len(), 21);
    }
}
