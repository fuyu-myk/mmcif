use std::path::{Path, PathBuf};

use mmcif_ir::SegmentId;
use mmcif_spec::{Checked, check_text};
use proptest::prelude::*;

fn fixtures() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures")
}

fn check(yaml: &str) -> Checked {
    check_text("design.yaml", yaml, &fixtures())
}

#[derive(Debug, Clone)]
enum Piece {
    Run { min: u32, max: u32 },
    Residue(char),
}

impl Piece {
    fn text(&self) -> String {
        match self {
            Piece::Run { min, max } if min == max => min.to_string(),
            Piece::Run { min, max } => format!("{min}..{max}"),
            Piece::Residue(letter) => letter.to_string(),
        }
    }

    fn shortest(&self) -> u32 {
        match self {
            Piece::Run { min, .. } => *min,
            Piece::Residue(_) => 1,
        }
    }
}

#[derive(Debug, Clone)]
struct Chain {
    id: String,
    pieces: Vec<Piece>,
}

impl Chain {
    fn sequence(&self) -> String {
        self.pieces.iter().map(Piece::text).collect()
    }

    fn shortest(&self) -> u32 {
        self.pieces.iter().map(Piece::shortest).sum()
    }

    fn cysteines(&self) -> Vec<u32> {
        let mut positions = Vec::new();
        let mut seen = 0;
        for piece in &self.pieces {
            if matches!(piece, Piece::Residue('C')) {
                positions.push(seen + 1);
            }
            seen += piece.shortest();
        }
        positions
    }
}

fn piece() -> impl Strategy<Value = Piece> {
    prop_oneof![
        (1u32..12, 0u32..6).prop_map(|(min, span)| Piece::Run {
            min,
            max: min + span
        }),
        prop::sample::select(vec!['A', 'C', 'D', 'G', 'K', 'S', 'W', 'Y']).prop_map(Piece::Residue),
    ]
}

fn run() -> impl Strategy<Value = Piece> {
    (1u32..12, 0u32..6).prop_map(|(min, span)| Piece::Run {
        min,
        max: min + span,
    })
}

fn chain(id: &'static str) -> impl Strategy<Value = Chain> {
    prop::collection::vec(piece(), 1..6).prop_map(move |pieces| Chain {
        id: id.to_string(),
        pieces: writable(pieces),
    })
}

fn disulfide_chain(id: &'static str) -> impl Strategy<Value = Chain> {
    (run(), run(), run()).prop_map(move |(first, middle, last)| Chain {
        id: id.to_string(),
        pieces: vec![
            first,
            Piece::Residue('C'),
            middle,
            Piece::Residue('C'),
            last,
        ],
    })
}

fn writable(pieces: Vec<Piece>) -> Vec<Piece> {
    let mut kept: Vec<Piece> = Vec::new();
    for piece in pieces {
        let runs_together = matches!(
            (kept.last(), &piece),
            (Some(Piece::Run { .. }), Piece::Run { .. })
        );
        if !runs_together {
            kept.push(piece);
        }
    }
    kept
}

fn spec_of(chains: &[Chain], bonds: &[(String, u32, u32)]) -> String {
    let mut yaml = String::from("entities:\n");
    for chain in chains {
        yaml.push_str(&format!(
            "  - protein:\n      id: {}\n      sequence: {}\n",
            chain.id,
            chain.sequence()
        ));
    }
    if !bonds.is_empty() {
        yaml.push_str("constraints:\n");
        for (id, first, second) in bonds {
            yaml.push_str(&format!(
                "  - bond:\n      atom1: [{id}, {first}, SG]\n      atom2: [{id}, {second}, SG]\n"
            ));
        }
    }
    yaml
}

proptest! {
    #[test]
    fn arbitrary_text_is_never_a_panic(text in ".{0,400}") {
        let checked = check(&text);
        prop_assert!(checked.has_errors() || checked.shortest.is_some());
    }

    #[test]
    fn arbitrary_yaml_shaped_text_is_never_a_panic(
        keys in prop::collection::vec(
            prop::sample::select(vec![
                "entities", "constraints", "protein", "ligand", "file", "id", "sequence",
                "path", "include", "design", "res_index", "bond", "atom1", "chain",
            ]),
            1..8,
        ),
        values in prop::collection::vec(
            prop::sample::select(vec!["A", "17", "3..5C6C3", "[A, 1, SG]", "all", "2..50,55..", "-"]),
            1..8,
        ),
    ) {
        let mut yaml = String::new();
        for (depth, key) in keys.iter().enumerate() {
            let value = values[depth % values.len()];
            yaml.push_str(&format!("{}{key}: {value}\n", "  ".repeat(depth % 3)));
        }
        let checked = check(&yaml);
        prop_assert!(checked.has_errors() || checked.shortest.is_some());
    }

    #[test]
    fn every_diagnostic_points_inside_its_source(text in ".{0,200}") {
        let checked = check(&text);
        for diagnostic in checked.report.diagnostics() {
            for label in &diagnostic.labels {
                let Some(file) = checked.sources.file(label.source) else {
                    continue;
                };
                prop_assert!(
                    label.span.end() <= file.text().len(),
                    "{:?} ends past the {} characters of {}",
                    label.span,
                    file.text().len(),
                    file.name()
                );
            }
        }
    }

    #[test]
    fn a_generated_specification_checks_clean(chains in prop::collection::vec(chain("S"), 1..2)) {
        let chains: Vec<Chain> = chains
            .into_iter()
            .enumerate()
            .map(|(index, mut chain)| {
                chain.id = format!("S{index}");
                chain
            })
            .collect();
        let checked = check(&spec_of(&chains, &[]));
        prop_assert!(!checked.has_errors(), "{}", checked.render());
        let shortest = checked.shortest.expect("a clean spec instantiates");
        let expected: u32 = chains.iter().map(Chain::shortest).sum();
        prop_assert_eq!(shortest.slot_count(), expected as usize);
        for chain in &chains {
            let span = shortest
                .segment_slots(&SegmentId::new(chain.id.as_str()))
                .expect("every chain becomes a segment");
            prop_assert_eq!(span.len(), chain.shortest() as usize);
        }
    }

    #[test]
    fn a_bond_lands_on_the_position_the_pattern_puts_it(chain in disulfide_chain("S")) {
        let cysteines = chain.cysteines();
        let (first, second) = (cysteines[0], cysteines[1]);
        let checked = check(&spec_of(
            std::slice::from_ref(&chain),
            &[(chain.id.clone(), first, second)],
        ));
        prop_assert!(!checked.has_errors(), "{}", checked.render());
        let shortest = checked.shortest.expect("a clean spec instantiates");
        prop_assert!(
            shortest
                .is_constrained(mmcif_ir::SlotIndex::new(first - 1))
                .unwrap_or(false),
            "slot {} should carry the bond",
            first - 1
        );
        prop_assert!(
            shortest
                .is_constrained(mmcif_ir::SlotIndex::new(second - 1))
                .unwrap_or(false),
            "slot {} should carry the bond",
            second - 1
        );
    }

    #[test]
    fn an_index_past_the_shortest_length_is_always_refused(
        chain in chain("S"),
        beyond in 1u32..50,
    ) {
        let index = chain.shortest() + beyond;
        let checked = check(&spec_of(
            std::slice::from_ref(&chain),
            &[(chain.id.clone(), 1, index)],
        ));
        prop_assert!(
            checked.report.codes().contains(&"unknown-residue"),
            "index {index} past {} was accepted: {:?}",
            chain.shortest(),
            checked.report.codes()
        );
    }
}
