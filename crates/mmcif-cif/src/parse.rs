use std::collections::HashSet;

use crate::document::{Cell, RawBlock, RawEntry, RawFrame, RawLoop, Span, StreamedBlock};
use crate::error::{ParseError, ParseErrorKind};
use crate::token::{Token, TokenKind, Tokenizer};

pub(crate) struct Parser<'a> {
    src: &'a str,
    tokens: Tokenizer<'a>,
    peeked: Option<Token>,
    failed: bool,
}

impl<'a> Parser<'a> {
    pub fn new(src: &'a str) -> Result<Self, ParseError> {
        if u32::try_from(src.len()).is_err() {
            return Err(ParseError::at(src.as_bytes(), 0, ParseErrorKind::TooLarge));
        }
        Ok(Parser {
            src,
            tokens: Tokenizer::new(src.as_bytes()),
            peeked: None,
            failed: false,
        })
    }

    fn error(&self, offset: usize, kind: ParseErrorKind) -> ParseError {
        ParseError::at(self.src.as_bytes(), offset, kind)
    }

    fn peek(&mut self) -> Result<Option<Token>, ParseError> {
        if self.peeked.is_none() {
            self.peeked = self.tokens.next_token()?;
        }
        Ok(self.peeked)
    }

    fn advance(&mut self) -> Result<Option<Token>, ParseError> {
        match self.peeked.take() {
            Some(token) => Ok(Some(token)),
            None => self.tokens.next_token(),
        }
    }

    fn span(token: Token) -> Span {
        Span {
            start: token.start as u32,
            len: (token.end - token.start) as u32,
        }
    }

    fn cell(&self, token: Token) -> Result<Cell, ParseError> {
        let TokenKind::Value(kind) = token.kind else {
            return Err(self.error(token.start, ParseErrorKind::MissingValue));
        };
        Cell::new(token.start, token.end, kind)
            .ok_or_else(|| self.error(token.start, ParseErrorKind::TooLarge))
    }

    fn next_block(&mut self) -> Result<Option<RawBlock>, ParseError> {
        let Some(header) = self.advance()? else {
            return Ok(None);
        };
        if header.kind != TokenKind::Data {
            return Err(self.error(header.start, ParseErrorKind::OutsideBlock));
        }
        if header.start == header.end {
            return Err(self.error(header.start, ParseErrorKind::EmptyBlockName));
        }
        let mut entries = Vec::new();
        let mut frame: Option<(Token, Vec<RawEntry>)> = None;
        let mut block_tags = HashSet::new();
        let mut frame_tags = HashSet::new();
        while let Some(token) = self.peek()? {
            if token.kind == TokenKind::Data {
                break;
            }
            self.advance()?;
            match token.kind {
                TokenKind::Tag => {
                    let value = match self.advance()? {
                        Some(value) if matches!(value.kind, TokenKind::Value(_)) => {
                            self.cell(value)?
                        }
                        _ => return Err(self.error(token.start, ParseErrorKind::MissingValue)),
                    };
                    let tags = if frame.is_some() {
                        &mut frame_tags
                    } else {
                        &mut block_tags
                    };
                    self.claim(tags, Self::span(token))?;
                    let target = frame
                        .as_mut()
                        .map_or(&mut entries, |(_, frame_entries)| frame_entries);
                    target.push(RawEntry::Pair {
                        tag: Self::span(token),
                        value,
                    });
                }
                TokenKind::Loop => {
                    let raw_loop = self.parse_loop(token)?;
                    let tags = if frame.is_some() {
                        &mut frame_tags
                    } else {
                        &mut block_tags
                    };
                    for &tag in &raw_loop.tags {
                        self.claim(tags, tag)?;
                    }
                    let target = frame
                        .as_mut()
                        .map_or(&mut entries, |(_, frame_entries)| frame_entries);
                    target.push(RawEntry::Loop(raw_loop));
                }
                TokenKind::Save if token.start == token.end => match frame.take() {
                    Some((open, frame_entries)) => entries.push(RawEntry::Frame(RawFrame {
                        name: Self::span(open),
                        entries: frame_entries,
                    })),
                    None => return Err(self.error(token.start, ParseErrorKind::UnmatchedSaveEnd)),
                },
                TokenKind::Save => {
                    if frame.is_some() {
                        return Err(self.error(token.start, ParseErrorKind::NestedSaveFrame));
                    }
                    frame = Some((token, Vec::new()));
                    frame_tags.clear();
                }
                TokenKind::Value(_) => {
                    return Err(self.error(token.start, ParseErrorKind::UnexpectedValue));
                }
                TokenKind::Global => {
                    return Err(self.error(token.start, ParseErrorKind::Reserved("global_")));
                }
                TokenKind::Stop => {
                    return Err(self.error(token.start, ParseErrorKind::Reserved("stop_")));
                }
                TokenKind::Data => unreachable!("data headers end the block before advancing"),
            }
        }
        if let Some((open, _)) = frame {
            return Err(self.error(open.start, ParseErrorKind::UnterminatedSaveFrame));
        }
        Ok(Some(RawBlock {
            name: Self::span(header),
            entries,
        }))
    }

    fn claim(&self, seen: &mut HashSet<String>, tag: Span) -> Result<(), ParseError> {
        let text = tag.resolve(self.src);
        if seen.insert(text.to_ascii_lowercase()) {
            Ok(())
        } else {
            Err(self.error(
                tag.start as usize,
                ParseErrorKind::DuplicateTag(text.to_owned()),
            ))
        }
    }

    fn parse_loop(&mut self, keyword: Token) -> Result<RawLoop, ParseError> {
        let mut tags = Vec::new();
        while let Some(token) = self.peek()?.filter(|t| t.kind == TokenKind::Tag) {
            self.advance()?;
            tags.push(Self::span(token));
        }
        if tags.is_empty() {
            return Err(self.error(keyword.start, ParseErrorKind::EmptyLoop));
        }
        let mut cells = Vec::new();
        while let Some(token) = self
            .peek()?
            .filter(|t| matches!(t.kind, TokenKind::Value(_)))
        {
            self.advance()?;
            cells.push(self.cell(token)?);
        }
        if cells.len() % tags.len() != 0 {
            return Err(self.error(
                keyword.start,
                ParseErrorKind::RaggedLoop {
                    values: cells.len(),
                    tags: tags.len(),
                },
            ));
        }
        Ok(RawLoop { tags, cells })
    }
}

impl Iterator for Parser<'_> {
    type Item = Result<RawBlock, ParseError>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.failed {
            return None;
        }
        let result = self.next_block().transpose();
        self.failed = matches!(result, Some(Err(_)));
        result
    }
}

pub struct BlockStream<'a> {
    parser: Option<Parser<'a>>,
    error: Option<ParseError>,
}

impl<'a> BlockStream<'a> {
    pub fn new(src: &'a str) -> Self {
        match Parser::new(src) {
            Ok(parser) => BlockStream {
                parser: Some(parser),
                error: None,
            },
            Err(error) => BlockStream {
                parser: None,
                error: Some(error),
            },
        }
    }
}

impl<'a> Iterator for BlockStream<'a> {
    type Item = Result<StreamedBlock<'a>, ParseError>;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(error) = self.error.take() {
            return Some(Err(error));
        }
        let parser = self.parser.as_mut()?;
        let src = parser.src;
        Some(parser.next()?.map(|raw| StreamedBlock::new(src, raw)))
    }
}
