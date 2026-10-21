use mmcif_ir::SecondaryStructure;

use crate::diagnostic::{Diagnostic, Report};
use crate::document::{
    AtomRef, Bond, Document, Entity, EntityKind, FileEntity, Ligand, Protein, Range, Ranges,
    Selection, Spanned, StructureSelection,
};
use crate::sequence;
use crate::span::{SourceId, Span};
use crate::yaml::{Mapping, Node};

const ROOT_KEYS: [&str; 2] = ["entities", "constraints"];
const FILE_KEYS: [&str; 4] = ["path", "include", "design", "secondary_structure"];
const FILE_UNSUPPORTED: [&str; 2] = ["binding_types", "structure_groups"];
const CHAIN_KEYS: [&str; 2] = ["id", "res_index"];
const CHAIN_UNSUPPORTED: [&str; 1] = ["symmetric_group"];
const ENTITY_KINDS: [&str; 3] = ["protein", "ligand", "file"];
const ENTITY_UNSUPPORTED: [&str; 2] = ["dna", "rna"];

pub fn document(source: SourceId, root: &Node, report: &mut Report) -> Document {
    let mut reader = Reader { source, report };
    let mut document = Document {
        entities: Vec::new(),
        bonds: Vec::new(),
    };
    let Some(mapping) = reader.mapping(root) else {
        return document;
    };
    reader.keys(mapping, &ROOT_KEYS, &[]);
    match mapping.entry("entities") {
        None => reader.report.push(
            Diagnostic::error("missing-entities", "`entities` is missing").at(
                source,
                root.span,
                "a specification declares at least one entity",
            ),
        ),
        Some((_, node)) => {
            for item in reader.sequence(node).unwrap_or_default() {
                if let Some(entity) = reader.entity(item) {
                    document.entities.push(entity);
                }
            }
        }
    }
    if let Some(node) = mapping.get("constraints") {
        for item in reader.sequence(node).unwrap_or_default() {
            if let Some(bond) = reader.constraint(item) {
                document.bonds.push(bond);
            }
        }
    }
    document
}

struct Reader<'a> {
    source: SourceId,
    report: &'a mut Report,
}

impl Reader<'_> {
    fn mismatch<'t>(&mut self, node: &'t Node, wanted: &str) -> Option<&'t Node> {
        self.report.push(
            Diagnostic::error("wrong-shape", format!("expected {wanted}")).at(
                self.source,
                node.span,
                format!("this is {}", node.kind()),
            ),
        );
        None
    }

    fn mapping<'t>(&mut self, node: &'t Node) -> Option<&'t Mapping> {
        match node.as_mapping() {
            Some(mapping) => Some(mapping),
            None => {
                self.mismatch(node, "a mapping")?;
                None
            }
        }
    }

    fn sequence<'t>(&mut self, node: &'t Node) -> Option<&'t [Node]> {
        match node.as_sequence() {
            Some(items) => Some(items),
            None => {
                self.mismatch(node, "a list")?;
                None
            }
        }
    }

    fn scalar(&mut self, node: &Node) -> Option<Spanned<String>> {
        match node.as_str() {
            Some(text) => Some(Spanned::new(node.span, text.to_string())),
            None => {
                self.mismatch(node, "a value")?;
                None
            }
        }
    }

    fn required<'t>(&mut self, mapping: &'t Mapping, key: &str, span: Span) -> Option<&'t Node> {
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

    fn keys(&mut self, mapping: &Mapping, known: &[&str], unsupported: &[&str]) {
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

    fn single<'t>(&mut self, node: &'t Node, wanted: &str) -> Option<(Spanned<String>, &'t Node)> {
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

    fn entity(&mut self, node: &Node) -> Option<Entity> {
        let (kind, body) = self.single(node, "entity kind")?;
        let name = kind.value.as_str();
        if ENTITY_UNSUPPORTED.contains(&name) {
            self.report.push(
                Diagnostic::warning(
                    "unsupported-entity",
                    format!("`{name}` entities are skipped"),
                )
                .at(self.source, kind.span, "not represented yet")
                .help("protein, ligand and file entities are checked"),
            );
            return None;
        }
        if !ENTITY_KINDS.contains(&name) {
            let mut diagnostic = Diagnostic::error(
                "unknown-entity",
                format!("`{name}` is not an entity kind"),
            )
            .at(self.source, kind.span, "not part of the schema");
            if let Some(closest) = closest(name, &ENTITY_KINDS) {
                diagnostic = diagnostic.help(format!("did you mean `{closest}`?"));
            }
            self.report.push(diagnostic);
            return None;
        }
        let mapping = self.mapping(body)?;
        let kind = match name {
            "protein" => {
                self.keys(mapping, &["id", "sequence"], &["msa"]);
                let id = self.required(mapping, "id", body.span)?;
                let id = self.scalar(id)?;
                let node = self.required(mapping, "sequence", body.span)?;
                let text = self.scalar(node)?;
                let pieces = match sequence::parse(self.source, node.span, &text.value) {
                    Ok(pieces) => pieces,
                    Err(diagnostic) => {
                        self.report.push(diagnostic);
                        return None;
                    }
                };
                EntityKind::Protein(Protein {
                    id,
                    sequence: Spanned::new(node.span, pieces),
                })
            }
            "ligand" => {
                self.keys(mapping, &["id", "ccd"], &["smiles"]);
                let id = self.required(mapping, "id", body.span)?;
                let id = self.scalar(id)?;
                let ccd = self.required(mapping, "ccd", body.span)?;
                let ccd = self.scalar(ccd)?;
                EntityKind::Ligand(Ligand { id, ccd })
            }
            _ => EntityKind::File(self.file(body, mapping)?),
        };
        Some(Entity {
            span: node.span,
            kind,
        })
    }

    fn file(&mut self, body: &Node, mapping: &Mapping) -> Option<FileEntity> {
        self.keys(mapping, &FILE_KEYS, &FILE_UNSUPPORTED);
        let path = self.required(mapping, "path", body.span)?;
        let path = self.scalar(path)?;
        Some(FileEntity {
            path,
            include: self.selections(mapping.get("include")),
            design: self.selections(mapping.get("design")),
            structure: self.structure(mapping.get("secondary_structure")),
        })
    }

    fn selections(&mut self, node: Option<&Node>) -> Vec<Selection> {
        let Some(node) = node else {
            return Vec::new();
        };
        let mut selections = Vec::new();
        for item in self.sequence(node).unwrap_or_default() {
            let Some((kind, body)) = self.single(item, "selector") else {
                continue;
            };
            if kind.value != "chain" {
                self.report.push(
                    Diagnostic::error("unknown-key", format!("`{}` is not a selector", kind.value))
                        .at(self.source, kind.span, "expected `chain`"),
                );
                continue;
            }
            let Some(mapping) = self.mapping(body) else {
                continue;
            };
            self.keys(mapping, &CHAIN_KEYS, &CHAIN_UNSUPPORTED);
            let Some(id) = self
                .required(mapping, "id", body.span)
                .and_then(|node| self.scalar(node))
            else {
                continue;
            };
            let residues = mapping.get("res_index").and_then(|node| self.ranges(node));
            selections.push(Selection {
                span: item.span,
                chain: id,
                residues,
            });
        }
        selections
    }

    fn structure(&mut self, node: Option<&Node>) -> Vec<StructureSelection> {
        let Some(node) = node else {
            return Vec::new();
        };
        let kinds = [
            ("loop", SecondaryStructure::Loop),
            ("helix", SecondaryStructure::Helix),
            ("sheet", SecondaryStructure::Strand),
        ];
        let known: Vec<&str> = std::iter::once("id")
            .chain(kinds.iter().map(|(name, _)| *name))
            .collect();
        let mut selections = Vec::new();
        for item in self.sequence(node).unwrap_or_default() {
            let Some((selector, body)) = self.single(item, "selector") else {
                continue;
            };
            if selector.value != "chain" {
                self.report.push(
                    Diagnostic::error(
                        "unknown-key",
                        format!("`{}` is not a selector", selector.value),
                    )
                    .at(self.source, selector.span, "expected `chain`"),
                );
                continue;
            }
            let Some(mapping) = self.mapping(body) else {
                continue;
            };
            self.keys(mapping, &known, &[]);
            let Some(chain) = self
                .required(mapping, "id", body.span)
                .and_then(|node| self.scalar(node))
            else {
                continue;
            };
            let mut regions = Vec::new();
            for (name, kind) in kinds {
                if let Some(node) = mapping.get(name)
                    && let Some(ranges) = self.ranges(node)
                {
                    regions.push((kind, ranges));
                }
            }
            selections.push(StructureSelection { chain, regions });
        }
        selections
    }

    fn ranges(&mut self, node: &Node) -> Option<Spanned<Ranges>> {
        let text = self.scalar(node)?;
        if text.value == "all" {
            return Some(Spanned::new(node.span, Ranges::All));
        }
        let mut ranges = Vec::new();
        for part in text.value.split(',') {
            let part = part.trim();
            let bad = |message: &str| {
                Diagnostic::error("bad-range", format!("`{part}` is not a residue range"))
                    .at(self.source, node.span, message.to_string())
                    .help("write indices such as 14, 2..50 or 55..")
            };
            let range = match part.split_once("..") {
                None => part.parse::<u32>().ok().map(|index| Range {
                    min: index,
                    max: Some(index),
                }),
                Some((min, "")) => min.parse::<u32>().ok().map(|min| Range { min, max: None }),
                Some((min, max)) => match (min.parse::<u32>(), max.parse::<u32>()) {
                    (Ok(min), Ok(max)) if min <= max => Some(Range {
                        min,
                        max: Some(max),
                    }),
                    _ => None,
                },
            };
            match range {
                Some(range) if range.min >= 1 => ranges.push(range),
                Some(_) => {
                    self.report.push(bad("residue indices start at 1"));
                    return None;
                }
                None => {
                    self.report.push(bad("could not be read as a range"));
                    return None;
                }
            }
        }
        Some(Spanned::new(node.span, Ranges::List(ranges)))
    }

    fn constraint(&mut self, node: &Node) -> Option<Bond> {
        let (kind, body) = self.single(node, "constraint")?;
        if kind.value != "bond" {
            self.report.push(
                Diagnostic::warning(
                    "unsupported-constraint",
                    format!("`{}` constraints are skipped", kind.value),
                )
                .at(self.source, kind.span, "not represented yet")
                .help("bond constraints are checked"),
            );
            return None;
        }
        let mapping = self.mapping(body)?;
        self.keys(mapping, &["atom1", "atom2"], &[]);
        let first = self.required(mapping, "atom1", body.span)?;
        let first = self.atom(first)?;
        let second = self.required(mapping, "atom2", body.span)?;
        let second = self.atom(second)?;
        Some(Bond {
            span: node.span,
            first,
            second,
        })
    }

    fn atom(&mut self, node: &Node) -> Option<AtomRef> {
        let items = self.sequence(node)?;
        let [chain, index, atom] = items else {
            self.report.push(
                Diagnostic::error("wrong-shape", "an atom is [chain, residue, atom]").at(
                    self.source,
                    node.span,
                    format!("{} values are given", items.len()),
                ),
            );
            return None;
        };
        let chain = self.scalar(chain)?;
        let index = self.scalar(index)?;
        let atom = self.scalar(atom)?;
        let Ok(number) = index.value.parse::<u32>() else {
            self.report.push(
                Diagnostic::error(
                    "bad-residue-index",
                    format!("`{}` is not a residue index", index.value),
                )
                .at(self.source, index.span, "expected a number")
                .help("indices start at 1 and count the shortest sampled length"),
            );
            return None;
        };
        if number == 0 {
            self.report.push(
                Diagnostic::error("bad-residue-index", "residue indices start at 1").at(
                    self.source,
                    index.span,
                    "index 0 does not exist",
                ),
            );
            return None;
        }
        Some(AtomRef {
            span: node.span,
            chain,
            index: Spanned::new(index.span, number),
            atom,
        })
    }
}

fn closest<'k>(name: &str, known: &[&'k str]) -> Option<&'k str> {
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
    use crate::document::Ranges;
    use crate::span::Sources;
    use crate::yaml;

    fn read(text: &str) -> (Document, Report) {
        let mut sources = Sources::new();
        let id = sources.add("design.yaml", text);
        let root = yaml::parse(id, text).expect("the fixture parses");
        let mut report = Report::new();
        let document = document(id, &root, &mut report);
        (document, report)
    }

    #[test]
    fn a_protein_entity_carries_its_sequence() {
        let (document, report) =
            read("entities:\n  - protein:\n      id: B\n      sequence: 80..140\n");
        assert!(report.is_empty(), "{:?}", report.codes());
        let EntityKind::Protein(protein) = &document.entities[0].kind else {
            panic!("expected a protein");
        };
        assert_eq!(protein.id.value, "B");
        assert_eq!(protein.sequence.value.len(), 1);
    }

    #[test]
    fn a_file_entity_collects_its_selections() {
        let (document, report) = read(concat!(
            "entities:\n",
            "  - file:\n",
            "      path: 8r3a.cif\n",
            "      include:\n",
            "        - chain:\n",
            "            id: A\n",
            "            res_index: 2..50,55..\n",
            "        - chain:\n",
            "            id: B\n",
            "      design:\n",
            "        - chain:\n",
            "            id: A\n",
            "            res_index: 14..19\n",
            "      secondary_structure:\n",
            "        - chain:\n",
            "            id: A\n",
            "            helix: 15..17\n",
            "            sheet: 19\n",
        ));
        assert!(report.is_empty(), "{:?}", report.codes());
        let EntityKind::File(file) = &document.entities[0].kind else {
            panic!("expected a file entity");
        };
        assert_eq!(file.path.value, "8r3a.cif");
        assert_eq!(file.include.len(), 2);
        assert!(file.include[1].residues.is_none());
        assert_eq!(
            file.include[0].residues.as_ref().unwrap().value,
            Ranges::List(vec![
                Range {
                    min: 2,
                    max: Some(50)
                },
                Range { min: 55, max: None }
            ])
        );
        assert_eq!(file.design[0].chain.value, "A");
        assert_eq!(file.structure[0].regions.len(), 2);
        assert_eq!(file.structure[0].regions[0].0, SecondaryStructure::Helix);
    }

    #[test]
    fn bonds_read_their_three_part_atom_references() {
        let (document, report) = read(concat!(
            "entities:\n  - protein:\n      id: S\n      sequence: 10..14C6C3\n",
            "constraints:\n  - bond:\n      atom1: [S, 11, SG]\n      atom2: [S, 18, SG]\n",
        ));
        assert!(report.is_empty(), "{:?}", report.codes());
        let bond = &document.bonds[0];
        assert_eq!(bond.first.chain.value, "S");
        assert_eq!(bond.first.index.value, 11);
        assert_eq!(bond.second.atom.value, "SG");
    }

    #[test]
    fn an_unknown_key_suggests_the_closest_one() {
        let (_, report) = read("entities:\n  - protein:\n      id: B\n      sequance: 17\n");
        assert_eq!(report.codes(), ["unknown-key", "missing-key"]);
        assert_eq!(
            report.diagnostics()[0].help.as_deref(),
            Some("did you mean `sequence`?")
        );
    }

    #[test]
    fn recognised_but_unchecked_sections_warn_instead_of_failing() {
        let (document, report) = read(concat!(
            "entities:\n",
            "  - file:\n",
            "      path: 8r3a.cif\n",
            "      binding_types:\n",
            "        - chain:\n",
            "            id: A\n",
            "            binding: 5..7\n",
            "  - dna:\n",
            "      id: D\n",
            "      sequence: ATGC\n",
        ));
        assert_eq!(report.codes(), ["unsupported-key", "unsupported-entity"]);
        assert!(!report.has_errors());
        assert_eq!(document.entities.len(), 1);
    }

    #[test]
    fn a_bad_entity_does_not_hide_the_next_one() {
        let (document, report) = read(concat!(
            "entities:\n",
            "  - protein:\n      id: B\n      sequence: 20..15\n",
            "  - protein:\n      id: C\n      sequence: 17\n",
        ));
        assert_eq!(report.codes(), ["bad-sequence"]);
        assert_eq!(document.entities.len(), 1);
        let EntityKind::Protein(kept) = &document.entities[0].kind else {
            panic!("expected a protein");
        };
        assert_eq!(kept.id.value, "C");
    }

    #[test]
    fn missing_entities_and_wrong_shapes_are_reported() {
        let (document, report) = read("constraints: []\n");
        assert_eq!(report.codes(), ["missing-entities"]);
        assert!(document.entities.is_empty());

        let (_, report) = read("entities:\n  - protein: B\n");
        assert_eq!(report.codes(), ["wrong-shape"]);

        let (_, report) = read("entities: 3\n");
        assert_eq!(report.codes(), ["wrong-shape"]);
    }

    #[test]
    fn malformed_atom_references_name_the_problem() {
        let bad = |text: &str| {
            let (_, report) = read(text);
            report.codes().first().copied().unwrap_or("none")
        };
        let spec =
            "entities:\n  - protein:\n      id: S\n      sequence: 17\nconstraints:\n  - bond:\n";
        assert_eq!(
            bad(&format!(
                "{spec}      atom1: [S, 11]\n      atom2: [S, 12, SG]\n"
            )),
            "wrong-shape"
        );
        assert_eq!(
            bad(&format!(
                "{spec}      atom1: [S, x, SG]\n      atom2: [S, 12, SG]\n"
            )),
            "bad-residue-index"
        );
        assert_eq!(
            bad(&format!(
                "{spec}      atom1: [S, 0, SG]\n      atom2: [S, 12, SG]\n"
            )),
            "bad-residue-index"
        );
        assert_eq!(
            bad(&format!("{spec}      atom2: [S, 12, SG]\n")),
            "missing-key"
        );
    }

    #[test]
    fn bad_ranges_are_refused() {
        let (_, report) = read(concat!(
            "entities:\n  - file:\n      path: a.cif\n      include:\n",
            "        - chain:\n            id: A\n            res_index: 50..2\n",
        ));
        assert_eq!(report.codes(), ["bad-range"]);

        let (_, report) = read(concat!(
            "entities:\n  - file:\n      path: a.cif\n      include:\n",
            "        - chain:\n            id: A\n            res_index: 0..4\n",
        ));
        assert_eq!(
            report.diagnostics()[0].labels[0].message,
            "residue indices start at 1"
        );
    }

    #[test]
    fn other_constraint_kinds_warn_and_are_skipped() {
        let (document, report) = read(concat!(
            "entities:\n  - protein:\n      id: B\n      sequence: 17\n",
            "constraints:\n  - pocket:\n      binder: B\n",
        ));
        assert_eq!(report.codes(), ["unsupported-constraint"]);
        assert!(document.bonds.is_empty());
    }

    #[test]
    fn edit_distance_only_suggests_near_misses() {
        assert_eq!(closest("sequance", &["sequence", "id"]), Some("sequence"));
        assert_eq!(closest("totally-different", &["sequence", "id"]), None);
        assert_eq!(distance("", "abc"), 3);
        assert_eq!(distance("abc", "abc"), 0);
    }
}
