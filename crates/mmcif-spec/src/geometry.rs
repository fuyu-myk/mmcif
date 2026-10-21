use mmcif_ir::{Bounds, Constraint, Geometry, Instantiated, ResidueRef, SlotIndex};
use mmcif_model::LabelResidueKey;

use crate::check::span_of;
use crate::diagnostic::{Diagnostic, Report};
use crate::lower::{Addressing, Lowered};
use crate::resolve;
use crate::span::SourceId;

const BOND_MIN: f64 = 0.9;
const BOND_MAX: f64 = 2.5;
const RESIDUE_STEP: f64 = 3.8;
const SIDE_CHAIN_REACH: f64 = 10.0;
const TOLERANCE: f64 = 1e-6;
const NODE_BUDGET: usize = 256;

pub struct Endpoint {
    pub slot: SlotIndex,
    pub atom: String,
    pub anchor: Option<LabelResidueKey>,
    pub target: Option<usize>,
    pub position: Option<[f64; 3]>,
}

pub fn endpoints(lowered: &Lowered, instantiated: &Instantiated) -> Vec<Vec<Endpoint>> {
    instantiated
        .constraints()
        .iter()
        .map(|constraint| match constraint {
            Constraint::Bond { first, second } | Constraint::Distance { first, second, .. } => {
                vec![
                    endpoint(lowered, instantiated, first.slot, first.atom.as_str()),
                    endpoint(lowered, instantiated, second.slot, second.atom.as_str()),
                ]
            }
            _ => Vec::new(),
        })
        .collect()
}

fn endpoint(
    lowered: &Lowered,
    instantiated: &Instantiated,
    slot: SlotIndex,
    atom: &str,
) -> Endpoint {
    let anchor = instantiated
        .slot(slot)
        .ok()
        .and_then(|slot| match &slot.geometry {
            Geometry::Anchored(ResidueRef::Label(key)) => Some(key.clone()),
            _ => None,
        });
    let home = home(lowered, instantiated, slot);
    let position = anchor
        .as_ref()
        .zip(home.and_then(|index| lowered.targets.get(index)))
        .and_then(|(key, target)| resolve::position(target, key, atom));
    Endpoint {
        slot,
        atom: atom.to_string(),
        anchor,
        target: home,
        position,
    }
}

fn home(lowered: &Lowered, instantiated: &Instantiated, slot: SlotIndex) -> Option<usize> {
    let origin = instantiated.origin(slot).ok()?;
    let segment = lowered.segments.get(origin.segment.as_usize())?;
    segment.elements.get(origin.element.as_usize())?.target
}

pub fn target(lowered: &Lowered, index: Option<usize>) -> Option<&resolve::Target> {
    lowered.targets.get(index?)
}

struct Edge {
    from: usize,
    to: usize,
    constraint: usize,
}

struct Limits {
    count: usize,
    upper: Vec<f64>,
    lower: Vec<f64>,
}

impl Limits {
    fn new(count: usize) -> Limits {
        let mut limits = Limits {
            count,
            upper: vec![f64::INFINITY; count * count],
            lower: vec![0.0; count * count],
        };
        for node in 0..count {
            limits.upper[node * count + node] = 0.0;
        }
        limits
    }

    fn tighten(&mut self, first: usize, second: usize, min: f64, max: f64) {
        for index in [first * self.count + second, second * self.count + first] {
            self.upper[index] = self.upper[index].min(max);
            self.lower[index] = self.lower[index].max(min);
        }
    }

    fn close(&mut self) {
        for middle in 0..self.count {
            for first in 0..self.count {
                for second in 0..self.count {
                    let through = self.upper[first * self.count + middle]
                        + self.upper[middle * self.count + second];
                    if through < self.upper[first * self.count + second] {
                        self.upper[first * self.count + second] = through;
                        self.upper[second * self.count + first] = through;
                    }
                }
            }
        }
    }

    fn upper_at(&self, first: usize, second: usize) -> f64 {
        self.upper[first * self.count + second]
    }

    fn lower_at(&self, first: usize, second: usize) -> f64 {
        self.lower[first * self.count + second]
    }
}

pub fn feasibility(
    source: SourceId,
    lowered: &Lowered,
    instantiated: &Instantiated,
    ends: &[Vec<Endpoint>],
    report: &mut Report,
) {
    let nodes = nodes(ends);
    if nodes.is_empty() {
        return;
    }
    if nodes.len() > NODE_BUDGET {
        report.push(
            Diagnostic::warning(
                "geometry-not-checked",
                format!(
                    "{} constrained atoms is more than this check covers",
                    nodes.len()
                ),
            )
            .at(
                source,
                span_of(lowered, Some(0)),
                "geometry was left unchecked",
            )
            .help(format!("distances close over at most {NODE_BUDGET} atoms")),
        );
        return;
    }

    let count = nodes.len();
    let mut limits = Limits::new(count);
    for first in 0..count {
        for second in first + 1..count {
            match (nodes[first].position, nodes[second].position) {
                (Some((left, home)), Some((right, away))) if home == away => {
                    let measured = distance(left, right);
                    limits.tighten(first, second, measured, measured);
                }
                _ => {
                    if let Some(limit) =
                        reach(lowered, instantiated, nodes[first].slot, nodes[second].slot)
                    {
                        limits.tighten(first, second, 0.0, limit);
                    }
                }
            }
        }
    }

    let mut edges = Vec::new();
    for (constraint, pair) in instantiated.constraints().iter().zip(ends).enumerate() {
        let (window, pair) = match pair {
            (Constraint::Bond { .. }, ends) => (Bounds::new(BOND_MIN, BOND_MAX).ok(), ends),
            (Constraint::Distance { bounds, .. }, ends) => (Some(*bounds), ends),
            (_, ends) => (None, ends),
        };
        let (Some(window), [first, second]) = (window, pair.as_slice()) else {
            continue;
        };
        let (Some(from), Some(to)) = (place(&nodes, first), place(&nodes, second)) else {
            continue;
        };
        if from == to {
            report.push(
                Diagnostic::error(
                    "self-reference",
                    format!(
                        "{} cannot be constrained against itself",
                        name(lowered, instantiated, &nodes[from])
                    ),
                )
                .at(
                    source,
                    span_of(lowered, Some(constraint)),
                    "both ends name the same atom",
                ),
            );
            continue;
        }
        edges.push(Edge {
            from,
            to,
            constraint,
        });
        limits.tighten(from, to, window.min(), window.max());
    }

    limits.close();
    let graph = Graph {
        nodes,
        edges,
        limits,
    };
    for pair in graph.worst_per_group() {
        report.push(graph.contradiction(source, lowered, instantiated, pair));
    }
}

struct Graph {
    nodes: Vec<Node>,
    edges: Vec<Edge>,
    limits: Limits,
}

impl Graph {
    fn worst_per_group(&self) -> Vec<(usize, usize)> {
        let count = self.nodes.len();
        let groups = groups(count, &self.edges);
        let mut worst: Vec<Option<(f64, usize, usize)>> = vec![None; count];
        for (first, group) in groups.iter().enumerate() {
            for second in first + 1..count {
                let margin =
                    self.limits.lower_at(first, second) - self.limits.upper_at(first, second);
                if margin <= TOLERANCE {
                    continue;
                }
                if worst[*group]
                    .map(|(known, _, _)| margin > known)
                    .unwrap_or(true)
                {
                    worst[*group] = Some((margin, first, second));
                }
            }
        }
        worst
            .into_iter()
            .flatten()
            .map(|(_, first, second)| (first, second))
            .collect()
    }

    fn blame(&self, first: usize, second: usize) -> Option<usize> {
        self.edges
            .iter()
            .find(|edge| {
                (edge.from, edge.to) == (first, second) || (edge.from, edge.to) == (second, first)
            })
            .or_else(|| {
                self.edges.iter().find(|edge| {
                    [edge.from, edge.to]
                        .iter()
                        .any(|node| *node == first || *node == second)
                })
            })
            .map(|edge| edge.constraint)
    }

    fn contradiction(
        &self,
        source: SourceId,
        lowered: &Lowered,
        instantiated: &Instantiated,
        pair: (usize, usize),
    ) -> Diagnostic {
        let (first, second) = pair;
        let diagnostic = Diagnostic::error(
            "infeasible-geometry",
            format!(
                "{} and {} cannot satisfy every constraint",
                name(lowered, instantiated, &self.nodes[first]),
                name(lowered, instantiated, &self.nodes[second])
            ),
        )
        .at(
            source,
            span_of(lowered, self.blame(first, second)),
            format!(
                "needs at least {:.2} angstrom, allows at most {:.2}",
                self.limits.lower_at(first, second),
                self.limits.upper_at(first, second)
            ),
        );
        match self.nodes[first].position.zip(self.nodes[second].position) {
            Some(((left, home), (right, away))) if home == away => diagnostic.help(format!(
                "the structure file fixes them {:.2} angstrom apart",
                distance(left, right)
            )),
            _ => diagnostic.help("loosen or drop one of these constraints"),
        }
    }
}

struct Node {
    slot: SlotIndex,
    atom: String,
    position: Option<([f64; 3], usize)>,
}

fn nodes(ends: &[Vec<Endpoint>]) -> Vec<Node> {
    let mut nodes: Vec<Node> = Vec::new();
    for end in ends.iter().flatten() {
        if nodes
            .iter()
            .any(|node| node.slot == end.slot && node.atom == end.atom)
        {
            continue;
        }
        nodes.push(Node {
            slot: end.slot,
            atom: end.atom.clone(),
            position: end.position.zip(end.target),
        });
    }
    nodes
}

fn place(nodes: &[Node], end: &Endpoint) -> Option<usize> {
    nodes
        .iter()
        .position(|node| node.slot == end.slot && node.atom == end.atom)
}

fn groups(count: usize, edges: &[Edge]) -> Vec<usize> {
    let mut parent: Vec<usize> = (0..count).collect();
    for edge in edges {
        let (from, to) = (root(&mut parent, edge.from), root(&mut parent, edge.to));
        parent[from] = to;
    }
    (0..count).map(|node| root(&mut parent, node)).collect()
}

fn root(parent: &mut [usize], mut node: usize) -> usize {
    while parent[node] != node {
        let above = parent[node];
        parent[node] = parent[above];
        node = above;
    }
    node
}

fn name(lowered: &Lowered, instantiated: &Instantiated, node: &Node) -> String {
    let Ok(origin) = instantiated.origin(node.slot) else {
        return node.atom.clone();
    };
    let Some(segment) = lowered.segments.get(origin.segment.as_usize()) else {
        return node.atom.clone();
    };
    let index = segment
        .elements
        .get(origin.element.as_usize())
        .and_then(|element| element.label)
        .or_else(|| {
            instantiated
                .segment_slots(&segment.id)
                .ok()
                .map(|span| node.slot.get() - span.start().get() + 1)
        })
        .unwrap_or_else(|| node.slot.get());
    format!("{} {index} {}", segment.id, node.atom)
}

fn reach(
    lowered: &Lowered,
    instantiated: &Instantiated,
    first: SlotIndex,
    second: SlotIndex,
) -> Option<f64> {
    if first == second {
        return None;
    }
    let left = instantiated.origin(first).ok()?;
    let right = instantiated.origin(second).ok()?;
    if left.segment != right.segment {
        return None;
    }
    let segment = lowered.segments.get(left.segment.as_usize())?;
    let steps = match segment.kind.addressing() {
        Addressing::Label => {
            let from = segment.elements.get(left.element.as_usize())?.label?;
            let to = segment.elements.get(right.element.as_usize())?.label?;
            u64::from(from.abs_diff(to))
        }
        Addressing::Positional => {
            let definition = lowered.spec.segment(&segment.id).ok()?;
            let longest = |element: usize, offset: u32| -> u64 {
                definition
                    .elements()
                    .iter()
                    .take(element)
                    .map(|element| u64::from(element.length().max()))
                    .sum::<u64>()
                    + u64::from(offset)
            };
            longest(left.element.as_usize(), left.offset)
                .abs_diff(longest(right.element.as_usize(), right.offset))
        }
    };
    if steps == 0 {
        return None;
    }
    Some(steps as f64 * RESIDUE_STEP + 2.0 * SIDE_CHAIN_REACH)
}

fn distance(first: [f64; 3], second: [f64; 3]) -> f64 {
    let [x, y, z] = [
        first[0] - second[0],
        first[1] - second[1],
        first[2] - second[2],
    ];
    (x * x + y * y + z * z).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn distances_are_euclidean() {
        assert_eq!(distance([0.0, 0.0, 0.0], [3.0, 4.0, 0.0]), 5.0);
        assert_eq!(distance([1.0, 1.0, 1.0], [1.0, 1.0, 1.0]), 0.0);
    }

    #[test]
    fn a_bond_window_covers_real_bond_lengths() {
        let bounds = Bounds::new(BOND_MIN, BOND_MAX).unwrap();
        assert!(bounds.contains(1.33) && bounds.contains(2.05));
        assert!(!bounds.contains(0.5) && !bounds.contains(3.0));
    }

    #[test]
    fn closing_a_path_tightens_the_pair_it_spans() {
        let mut limits = Limits::new(3);
        limits.tighten(0, 1, 0.0, 2.5);
        limits.tighten(1, 2, 0.0, 2.5);
        limits.tighten(0, 2, 15.0, 15.0);
        limits.close();
        assert_eq!(limits.upper_at(0, 2), 5.0);
        assert_eq!(limits.lower_at(0, 2), 15.0);
        assert_eq!(limits.upper_at(0, 0), 0.0);
    }

    #[test]
    fn groups_follow_the_constraint_edges() {
        let edges = [
            Edge {
                from: 0,
                to: 1,
                constraint: 0,
            },
            Edge {
                from: 2,
                to: 3,
                constraint: 1,
            },
        ];
        let groups = groups(5, &edges);
        assert_eq!(groups[0], groups[1]);
        assert_eq!(groups[2], groups[3]);
        assert_ne!(groups[0], groups[2]);
        assert_ne!(groups[0], groups[4]);
    }

    #[test]
    fn only_the_worst_violation_of_each_group_is_kept() {
        let mut limits = Limits::new(4);
        limits.tighten(0, 1, 0.0, 2.0);
        limits.tighten(1, 2, 0.0, 2.0);
        limits.tighten(0, 2, 20.0, 20.0);
        limits.tighten(1, 3, 9.0, 1.0);
        limits.close();
        let graph = Graph {
            nodes: (0..4)
                .map(|slot| Node {
                    slot: SlotIndex::new(slot),
                    atom: "CA".to_string(),
                    position: None,
                })
                .collect(),
            edges: vec![
                Edge {
                    from: 0,
                    to: 1,
                    constraint: 0,
                },
                Edge {
                    from: 1,
                    to: 2,
                    constraint: 1,
                },
            ],
            limits,
        };
        assert_eq!(graph.worst_per_group(), [(0, 2)]);
        assert_eq!(graph.blame(0, 2), Some(0));
        assert_eq!(graph.blame(0, 1), Some(0));
        assert_eq!(graph.blame(1, 2), Some(1));
    }
}
