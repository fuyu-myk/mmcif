use crate::diagnostic::{Diagnostic, Report};
use crate::document::Spanned;
use crate::span::{SourceId, Span};
use crate::yaml::{Mapping, Node};

pub struct Reader<'a> {
    pub source: SourceId,
    pub report: &'a mut Report,
}

impl Reader<'_> {
    pub fn mismatch(&mut self, node: &Node, wanted: &str) {
        self.report.push(
            Diagnostic::error("wrong-shape", format!("expected {wanted}")).at(
                self.source,
                node.span,
                format!("this is {}", node.kind()),
            ),
        );
    }

    pub fn mapping<'t>(&mut self, node: &'t Node) -> Option<&'t Mapping> {
        match node.as_mapping() {
            Some(mapping) => Some(mapping),
            None => {
                self.mismatch(node, "a mapping");
                None
            }
        }
    }

    pub fn sequence<'t>(&mut self, node: &'t Node) -> Option<&'t [Node]> {
        match node.as_sequence() {
            Some(items) => Some(items),
            None => {
                self.mismatch(node, "a list");
                None
            }
        }
    }

    pub fn scalar(&mut self, node: &Node) -> Option<Spanned<String>> {
        match node.as_str() {
            Some(text) => Some(Spanned::new(node.span, text.to_string())),
            None => {
                self.mismatch(node, "a value");
                None
            }
        }
    }

    pub fn required<'t>(
        &mut self,
        mapping: &'t Mapping,
        key: &str,
        span: Span,
    ) -> Option<&'t Node> {
        match mapping.get(key) {
            Some(node) => Some(node),
            None => {
                self.report.push(
                    Diagnostic::error("missing-key", format!("`{key}` is missing")).at(
                        self.source,
                        span,
                        format!("`{key}` is required here"),
                    ),
                );
                None
            }
        }
    }

    pub fn keys(&mut self, mapping: &Mapping, known: &[&str], unsupported: &[&str]) {
        for (key, _) in mapping.entries() {
            let Some(name) = key.as_str() else { continue };
            if known.contains(&name) {
                continue;
            }
            if unsupported.contains(&name) {
                self.report.push(
                    Diagnostic::warning("unsupported-key", format!("`{name}` is not checked"))
                        .at(self.source, key.span, "parsed but not represented")
                        .help("the checker ignores this section"),
                );
                continue;
            }
            let mut diagnostic = Diagnostic::error(
                "unknown-key",
                format!("`{name}` is not a known key"),
            )
            .at(self.source, key.span, "not part of the schema");
            if let Some(closest) = closest(name, known) {
                diagnostic = diagnostic.help(format!("did you mean `{closest}`?"));
            }
            self.report.push(diagnostic);
        }
    }

    pub fn single<'t>(
        &mut self,
        node: &'t Node,
        wanted: &str,
    ) -> Option<(Spanned<String>, &'t Node)> {
        let mapping = self.mapping(node)?;
        match mapping.single() {
            Some((key, value)) => Some((self.scalar(key)?, value)),
            None => {
                self.report.push(
                    Diagnostic::error("wrong-shape", format!("expected one {wanted}")).at(
                        self.source,
                        node.span,
                        format!("{} keys are set here", mapping.entries().len()),
                    ),
                );
                None
            }
        }
    }
}

pub fn closest<'k>(name: &str, known: &[&'k str]) -> Option<&'k str> {
    known
        .iter()
        .map(|candidate| (distance(name, candidate), *candidate))
        .filter(|(score, _)| *score <= 2)
        .min_by_key(|(score, _)| *score)
        .map(|(_, candidate)| candidate)
}

fn distance(left: &str, right: &str) -> usize {
    let right: Vec<char> = right.chars().collect();
    let mut row: Vec<usize> = (0..=right.len()).collect();
    for (i, a) in left.chars().enumerate() {
        let mut previous = row[0];
        row[0] = i + 1;
        for (j, b) in right.iter().enumerate() {
            let cost = usize::from(a != *b);
            let next = (row[j] + 1).min(row[j + 1] + 1).min(previous + cost);
            previous = row[j + 1];
            row[j + 1] = next;
        }
    }
    row[right.len()]
}
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn edit_distance_only_suggests_near_misses() {
        assert_eq!(closest("sequance", &["sequence", "id"]), Some("sequence"));
        assert_eq!(closest("totally-different", &["sequence", "id"]), None);
        assert_eq!(distance("", "abc"), 3);
        assert_eq!(distance("abc", "abc"), 0);
    }
}
