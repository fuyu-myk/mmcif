use saphyr_parser::{Event, Parser, Span as YamlSpan};

use crate::diagnostic::Diagnostic;
use crate::span::{SourceId, Span};

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Value {
    Scalar(String),
    Sequence(Vec<Node>),
    Mapping(Mapping),
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Node {
    pub span: Span,
    pub value: Value,
}

impl Node {
    pub fn as_str(&self) -> Option<&str> {
        match &self.value {
            Value::Scalar(text) => Some(text),
            _ => None,
        }
    }

    pub fn as_sequence(&self) -> Option<&[Node]> {
        match &self.value {
            Value::Sequence(items) => Some(items),
            _ => None,
        }
    }

    pub fn as_mapping(&self) -> Option<&Mapping> {
        match &self.value {
            Value::Mapping(mapping) => Some(mapping),
            _ => None,
        }
    }

    pub fn kind(&self) -> &'static str {
        match &self.value {
            Value::Scalar(_) => "a value",
            Value::Sequence(_) => "a list",
            Value::Mapping(_) => "a mapping",
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct Mapping {
    entries: Vec<(Node, Node)>,
}

impl Mapping {
    pub fn entries(&self) -> &[(Node, Node)] {
        &self.entries
    }

    pub fn entry(&self, key: &str) -> Option<(&Node, &Node)> {
        self.entries
            .iter()
            .find(|(known, _)| known.as_str() == Some(key))
            .map(|(known, value)| (known, value))
    }

    pub fn get(&self, key: &str) -> Option<&Node> {
        self.entry(key).map(|(_, value)| value)
    }

    pub fn single(&self) -> Option<(&Node, &Node)> {
        match self.entries.as_slice() {
            [(key, value)] => Some((key, value)),
            _ => None,
        }
    }
}

enum Frame {
    Sequence { start: usize, items: Vec<Node> },
    Mapping(MappingFrame),
}

struct MappingFrame {
    start: usize,
    entries: Vec<(Node, Node)>,
    pending: Option<Node>,
}

pub fn parse(source: SourceId, text: &str) -> Result<Node, Diagnostic> {
    let mut parser = Parser::new_from_str(text);
    let mut stack: Vec<Frame> = Vec::new();
    let mut document: Option<Node> = None;

    while let Some(event) = parser.next_event() {
        let (event, span) = event.map_err(|error| {
            let at = error.marker().index().min(text.len());
            Diagnostic::error("yaml-syntax", error.info().to_string()).at(
                source,
                Span::new(at, (at + 1).min(text.len())),
                "could not be parsed as YAML",
            )
        })?;
        match event {
            Event::StreamStart | Event::StreamEnd | Event::Nothing | Event::DocumentEnd => {}
            Event::DocumentStart(_) if document.is_some() => {
                return Err(Diagnostic::error(
                    "yaml-multiple-documents",
                    "a specification holds one YAML document",
                )
                .at(source, at(span), "a second document starts here"));
            }
            Event::DocumentStart(_) => {}
            Event::Alias(_) => {
                return Err(Diagnostic::error(
                    "yaml-alias",
                    "YAML anchors and aliases are not supported",
                )
                .at(source, at(span), "an alias is used here")
                .help("write the value out instead"));
            }
            Event::Scalar(text, _, anchor, _) => {
                if anchor != 0 {
                    return Err(anchor_unsupported(source, span));
                }
                let node = Node {
                    span: at(span),
                    value: Value::Scalar(text.into_owned()),
                };
                place(&mut stack, &mut document, node, source)?;
            }
            Event::SequenceStart(anchor, _) => {
                if anchor != 0 {
                    return Err(anchor_unsupported(source, span));
                }
                stack.push(Frame::Sequence {
                    start: span.start.index(),
                    items: Vec::new(),
                });
            }
            Event::MappingStart(anchor, _) => {
                if anchor != 0 {
                    return Err(anchor_unsupported(source, span));
                }
                stack.push(Frame::Mapping(MappingFrame {
                    start: span.start.index(),
                    entries: Vec::new(),
                    pending: None,
                }));
            }
            Event::SequenceEnd | Event::MappingEnd => {
                let node = close(stack.pop(), span);
                place(&mut stack, &mut document, node, source)?;
            }
        }
    }

    document.ok_or_else(|| {
        Diagnostic::error("empty-spec", "the specification is empty").at(
            source,
            Span::new(0, text.len().min(1)),
            "no entities are declared",
        )
    })
}

fn at(span: YamlSpan) -> Span {
    Span::new(span.start.index(), span.end.index())
}

fn anchor_unsupported(source: SourceId, span: YamlSpan) -> Diagnostic {
    Diagnostic::error("yaml-alias", "YAML anchors and aliases are not supported")
        .at(source, at(span), "an anchor is declared here")
        .help("write the value out instead")
}

fn close(frame: Option<Frame>, end: YamlSpan) -> Node {
    match frame {
        None => Node {
            span: at(end),
            value: Value::Mapping(Mapping::default()),
        },
        Some(Frame::Sequence { start, items }) => Node {
            span: cover(start, items.last().map(|item| item.span), end),
            value: Value::Sequence(items),
        },
        Some(Frame::Mapping(frame)) => Node {
            span: cover(
                frame.start,
                frame.entries.last().map(|(_, value)| value.span),
                end,
            ),
            value: Value::Mapping(Mapping {
                entries: frame.entries,
            }),
        },
    }
}

fn cover(start: usize, last: Option<Span>, end: YamlSpan) -> Span {
    let stop = last.map(Span::end).unwrap_or_else(|| end.end.index());
    Span::new(start, stop)
}

fn place(
    stack: &mut [Frame],
    document: &mut Option<Node>,
    node: Node,
    source: SourceId,
) -> Result<(), Diagnostic> {
    match stack.last_mut() {
        None => {
            *document = Some(node);
            Ok(())
        }
        Some(Frame::Sequence { items, .. }) => {
            items.push(node);
            Ok(())
        }
        Some(Frame::Mapping(frame)) => {
            match frame.pending.take() {
                None => {
                    if node.as_str().is_none() {
                        return Err(Diagnostic::error("yaml-complex-key", "keys must be names")
                            .at(source, node.span, format!("this key is {}", node.kind())));
                    }
                    if let Some(first) = duplicate(&frame.entries, &node) {
                        let name = node.as_str().unwrap_or_default().to_string();
                        return Err(Diagnostic::error(
                            "duplicate-key",
                            format!("`{name}` is set twice"),
                        )
                        .at(source, first, "first set here")
                        .at(source, node.span, "set again here"));
                    }
                    frame.pending = Some(node);
                    Ok(())
                }
                Some(key) => {
                    frame.entries.push((key, node));
                    Ok(())
                }
            }
        }
    }
}

fn duplicate(entries: &[(Node, Node)], key: &Node) -> Option<Span> {
    entries
        .iter()
        .find(|(known, _)| known.as_str() == key.as_str())
        .map(|(known, _)| known.span)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::span::Sources;

    fn load(text: &str) -> Result<Node, Diagnostic> {
        let mut sources = Sources::new();
        let id = sources.add("design.yaml", text);
        parse(id, text)
    }

    #[test]
    fn a_nested_document_keeps_every_span() {
        let text = "entities:\n  - protein:\n      id: B\n      sequence: 80..140\n";
        let document = load(text).unwrap();
        let entities = document.as_mapping().unwrap().get("entities").unwrap();
        let entity = &entities.as_sequence().unwrap()[0];
        let (kind, body) = entity.as_mapping().unwrap().single().unwrap();
        assert_eq!(kind.as_str(), Some("protein"));
        let sequence = body.as_mapping().unwrap().get("sequence").unwrap();
        assert_eq!(sequence.as_str(), Some("80..140"));
        assert_eq!(&text[sequence.span.start()..sequence.span.end()], "80..140");
        let id = body.as_mapping().unwrap().entry("id").unwrap().0.span;
        assert_eq!(&text[id.start()..id.end()], "id");
    }

    #[test]
    fn a_mapping_span_covers_its_entries() {
        let text = "a:\n  b: 1\n  c: 2\n";
        let document = load(text).unwrap();
        let inner = document.as_mapping().unwrap().get("a").unwrap();
        assert_eq!(&text[inner.span.start()..inner.span.end()], "b: 1\n  c: 2");
    }

    #[test]
    fn flow_collections_parse_like_block_ones() {
        let document = load("constraints:\n  - bond:\n      atom1: [R, 4, SG]\n").unwrap();
        let bond = &document
            .as_mapping()
            .unwrap()
            .get("constraints")
            .unwrap()
            .as_sequence()
            .unwrap()[0];
        let atom = bond
            .as_mapping()
            .unwrap()
            .single()
            .unwrap()
            .1
            .as_mapping()
            .unwrap()
            .get("atom1")
            .unwrap();
        let parts: Vec<&str> = atom
            .as_sequence()
            .unwrap()
            .iter()
            .map(|item| item.as_str().unwrap())
            .collect();
        assert_eq!(parts, ["R", "4", "SG"]);
    }

    #[test]
    fn a_repeated_key_reports_both_places() {
        let error = load("entities: []\nentities: []\n").unwrap_err();
        assert_eq!(error.code, "duplicate-key");
        assert_eq!(error.labels.len(), 2);
        assert_eq!(error.labels[0].span, Span::new(0, 8));
        assert_eq!(error.labels[1].span, Span::new(13, 21));
    }

    #[test]
    fn broken_yaml_reports_where_it_broke() {
        let error = load("entities:\n  - protein:\n   id: B\n").unwrap_err();
        assert_eq!(error.code, "yaml-syntax");
        assert_eq!(error.labels.len(), 1);
    }

    #[test]
    fn anchors_aliases_and_extra_documents_are_refused() {
        assert_eq!(load("a: &x 1\nb: *x\n").unwrap_err().code, "yaml-alias");
        assert_eq!(
            load("---\na: 1\n---\nb: 2\n").unwrap_err().code,
            "yaml-multiple-documents"
        );
        assert_eq!(
            load("? [a, b]\n: 1\n").unwrap_err().code,
            "yaml-complex-key"
        );
    }

    #[test]
    fn an_empty_file_is_not_a_specification() {
        let error = load("").unwrap_err();
        assert_eq!(error.code, "empty-spec");
        assert_eq!(error.labels[0].span, Span::new(0, 0));
    }

    #[test]
    fn scalars_keep_their_text_without_interpretation() {
        let document = load("a: 17\nb: \"80..140\"\nc: all\n").unwrap();
        let mapping = document.as_mapping().unwrap();
        assert_eq!(mapping.get("a").unwrap().as_str(), Some("17"));
        assert_eq!(mapping.get("b").unwrap().as_str(), Some("80..140"));
        assert_eq!(mapping.get("c").unwrap().as_str(), Some("all"));
        let keys: Vec<&str> = mapping
            .entries()
            .iter()
            .map(|(key, _)| key.as_str().unwrap())
            .collect();
        assert_eq!(keys, ["a", "b", "c"]);
        assert_eq!(mapping.get("a").unwrap().kind(), "a value");
    }
}
