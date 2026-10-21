use mmcif_ir::{
    Bounds, Constraint, Geometry, Instantiated, IrError, Lengths, MAX_SLOTS, ResidueRef, Resolved,
    SlotIndex,
};
use mmcif_model::LabelResidueKey;

use crate::diagnostic::{Diagnostic, Report};
use crate::lower::Lowered;
use crate::resolve;
use crate::span::{SourceId, Span};

const BOND_MIN: f64 = 0.9;
const BOND_MAX: f64 = 2.5;
const RESIDUE_STEP: f64 = 3.8;
const SIDE_CHAIN_REACH: f64 = 10.0;
const TOLERANCE: f64 = 1e-6;

pub fn check(source: SourceId, lowered: &Lowered, report: &mut Report) -> Option<Instantiated> {
    longest_fits(source, lowered, report)?;
    let shortest = instantiate(source, lowered, &Lengths::minimum(&lowered.spec), report)?;
    atoms(source, lowered, &shortest, report);
    feasibility(source, lowered, &shortest, report);
    Some(shortest)
}

fn longest_fits(source: SourceId, lowered: &Lowered, report: &mut Report) -> Option<()> {
    let longest = lowered.spec.length_bounds().1;
    if longest <= u64::from(MAX_SLOTS) {
        return Some(());
    }
    let span = lowered
        .segments
        .first()
        .map(|segment| segment.span)
        .unwrap_or(Span::new(0, 0));
    report.push(
        Diagnostic::error(
            "unsatisfiable-lengths",
            format!("the longest sampled length is {longest} slots, over the limit of {MAX_SLOTS}"),
        )
        .at(source, span, "no sampled length satisfies this")
        .help("narrow the length ranges"),
    );
    None
}

fn instantiate(
    source: SourceId,
    lowered: &Lowered,
    lengths: &Lengths,
    report: &mut Report,
) -> Option<Instantiated> {
    match lowered.spec.instantiate(lengths) {
        Ok(instantiated) => Some(instantiated),
        Err(error) => {
            report.push(unsatisfiable(source, lowered, error));
            None
        }
    }
}

fn unsatisfiable(source: SourceId, lowered: &Lowered, error: IrError) -> Diagnostic {
    let span = match &error {
        IrError::OffsetOutOfRange {
            segment, element, ..
        }
        | IrError::UnknownElement { segment, element }
        | IrError::MixedDeterminacy { segment, element } => lowered
            .segment(segment.as_str())
            .and_then(|info| info.elements.get(element.as_usize()))
            .map(|element| element.span),
        IrError::UnknownSegment(segment) => lowered.segment(segment.as_str()).map(|info| info.span),
        _ => None,
    };
    let span = span
        .or_else(|| lowered.segments.first().map(|segment| segment.span))
        .unwrap_or(Span::new(0, 0));
    Diagnostic::error("unsatisfiable-lengths", error.to_string()).at(
        source,
        span,
        "no sampled length satisfies this",
    )
}

struct Endpoint {
    slot: SlotIndex,
    atom: String,
    anchor: Option<LabelResidueKey>,
    position: Option<[f64; 3]>,
}

fn endpoints(lowered: &Lowered, instantiated: &Instantiated) -> Vec<Vec<Endpoint>> {
    instantiated
        .constraints()
        .iter()
        .map(|constraint| match constraint {
            Constraint::Bond { first, second } => vec![
                endpoint(lowered, instantiated, first.slot, first.atom.as_str()),
                endpoint(lowered, instantiated, second.slot, second.atom.as_str()),
            ],
            Constraint::Distance { first, second, .. } => vec![
                endpoint(lowered, instantiated, first.slot, first.atom.as_str()),
                endpoint(lowered, instantiated, second.slot, second.atom.as_str()),
            ],
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
    let position = anchor
        .as_ref()
        .zip(target(lowered, instantiated, slot))
        .and_then(|(key, target)| resolve::position(target, key, atom));
    Endpoint {
        slot,
        atom: atom.to_string(),
        anchor,
        position,
    }
}

fn target<'a>(
    lowered: &'a Lowered,
    instantiated: &Instantiated,
    slot: SlotIndex,
) -> Option<&'a resolve::Target> {
    let origin = instantiated.origin(slot).ok()?;
    let segment = lowered.segments.get(origin.segment.as_usize())?;
    let index = segment.elements.get(origin.element.as_usize())?.target?;
    lowered.targets.get(index)
}

fn atoms(source: SourceId, lowered: &Lowered, instantiated: &Instantiated, report: &mut Report) {
    for (index, ends) in endpoints(lowered, instantiated).iter().enumerate() {
        for end in ends {
            let Some(key) = &end.anchor else { continue };
            let Some(target) = target(lowered, instantiated, end.slot) else {
                continue;
            };
            if resolve::atom_exists(target, key, &end.atom) == Some(false) {
                let component = resolve::component(target, key)
                    .map(|comp| comp.as_str().to_string())
                    .unwrap_or_default();
                let mut names = resolve::atom_names(target, key);
                names.truncate(12);
                report.push(
                    Diagnostic::error(
                        "unknown-atom",
                        format!("{} has no atom {}", component, end.atom),
                    )
                    .at(
                        source,
                        span_of(lowered, index),
                        format!("residue {key} is {component}"),
                    )
                    .help(format!("atoms in the file: {}", names.join(", "))),
                );
            }
        }
    }
}

fn feasibility(
    source: SourceId,
    lowered: &Lowered,
    instantiated: &Instantiated,
    report: &mut Report,
) {
    let ends = endpoints(lowered, instantiated);
    let mut nodes: Vec<(SlotIndex, String)> = Vec::new();
    for pair in &ends {
        for end in pair {
            let key = (end.slot, end.atom.clone());
            if !nodes.contains(&key) {
                nodes.push(key);
            }
        }
    }
    if nodes.len() < 2 {
        return;
    }
    let count = nodes.len();
    let index_of = |slot: SlotIndex, atom: &str| {
        nodes
            .iter()
            .position(|(known, name)| *known == slot && name == atom)
    };
    let mut upper = vec![f64::INFINITY; count * count];
    let mut lower = vec![0.0f64; count * count];
    for (position, _) in nodes.iter().enumerate() {
        upper[position * count + position] = 0.0;
    }

    let known: Vec<Option<[f64; 3]>> = nodes
        .iter()
        .map(|(slot, atom)| {
            ends.iter()
                .flatten()
                .find(|end| end.slot == *slot && end.atom == *atom)
                .and_then(|end| end.position)
        })
        .collect();
    for first in 0..count {
        for second in first + 1..count {
            if let (Some(left), Some(right)) = (known[first], known[second]) {
                let measured = distance(left, right);
                tighten(
                    &mut upper, &mut lower, count, first, second, measured, measured,
                );
            } else if let Some(limit) = reach(instantiated, nodes[first].0, nodes[second].0) {
                tighten(&mut upper, &mut lower, count, first, second, 0.0, limit);
            }
        }
    }

    let mut edges = Vec::new();
    for (constraint, pair) in instantiated.constraints().iter().zip(&ends) {
        let bounds = match constraint {
            Constraint::Bond { .. } => Bounds::new(BOND_MIN, BOND_MAX).ok(),
            Constraint::Distance { bounds, .. } => Some(*bounds),
            _ => None,
        };
        let (Some(bounds), [first, second]) = (bounds, pair.as_slice()) else {
            continue;
        };
        let (Some(from), Some(to)) = (
            index_of(first.slot, &first.atom),
            index_of(second.slot, &second.atom),
        ) else {
            continue;
        };
        edges.push((from, to, bounds));
        tighten(
            &mut upper,
            &mut lower,
            count,
            from,
            to,
            bounds.min(),
            bounds.max(),
        );
    }

    for middle in 0..count {
        for first in 0..count {
            for second in 0..count {
                let through = upper[first * count + middle] + upper[middle * count + second];
                if through < upper[first * count + second] {
                    upper[first * count + second] = through;
                    upper[second * count + first] = through;
                }
            }
        }
    }

    let mut owner = vec![None; count];
    for (index, (from, to, _)) in edges.iter().enumerate() {
        for node in [*from, *to] {
            owner[node] = owner[node].or(Some(index));
        }
    }

    for first in 0..count {
        for second in first + 1..count {
            let closed = upper[first * count + second];
            let floor = lower[first * count + second];
            if floor <= closed + TOLERANCE {
                continue;
            }
            let edge = owner[first].or(owner[second]).unwrap_or(0);
            let (slot, atom) = &nodes[first];
            let (other, other_atom) = &nodes[second];
            let diagnostic = Diagnostic::error(
                "infeasible-geometry",
                format!("the constraints on {atom} and {other_atom} cannot all hold"),
            )
            .at(
                source,
                span_of(lowered, constraint_at(lowered, instantiated, edge)),
                format!("needs at least {floor:.2} angstrom, allows at most {closed:.2}"),
            );
            report.push(match known[first].zip(known[second]) {
                Some((left, right)) => diagnostic.help(format!(
                    "the structure file fixes slots {} and {} {:.2} angstrom apart",
                    slot.get(),
                    other.get(),
                    distance(left, right)
                )),
                None => diagnostic.help("loosen or drop one of these constraints"),
            });
        }
    }
}

fn constraint_at(lowered: &Lowered, instantiated: &Instantiated, edge: usize) -> usize {
    instantiated
        .constraints()
        .iter()
        .enumerate()
        .filter(|(_, constraint)| {
            matches!(
                constraint,
                Constraint::<Resolved>::Bond { .. } | Constraint::<Resolved>::Distance { .. }
            )
        })
        .nth(edge)
        .map(|(index, _)| index)
        .unwrap_or_else(|| lowered.constraint_spans.len())
}

fn tighten(
    upper: &mut [f64],
    lower: &mut [f64],
    count: usize,
    first: usize,
    second: usize,
    min: f64,
    max: f64,
) {
    let forward = first * count + second;
    let backward = second * count + first;
    for index in [forward, backward] {
        upper[index] = upper[index].min(max);
        lower[index] = lower[index].max(min);
    }
}

fn reach(instantiated: &Instantiated, first: SlotIndex, second: SlotIndex) -> Option<f64> {
    let left = instantiated.origin(first).ok()?;
    let right = instantiated.origin(second).ok()?;
    if left.segment != right.segment {
        return None;
    }
    let steps = first.get().abs_diff(second.get());
    Some(f64::from(steps) * RESIDUE_STEP + 2.0 * SIDE_CHAIN_REACH)
}

fn distance(first: [f64; 3], second: [f64; 3]) -> f64 {
    let [x, y, z] = [
        first[0] - second[0],
        first[1] - second[1],
        first[2] - second[2],
    ];
    (x * x + y * y + z * z).sqrt()
}

fn span_of(lowered: &Lowered, constraint: usize) -> Span {
    lowered
        .constraint_spans
        .get(constraint)
        .copied()
        .or_else(|| lowered.segments.first().map(|segment| segment.span))
        .unwrap_or(Span::new(0, 0))
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
        assert!(bounds.contains(1.33));
        assert!(bounds.contains(2.05));
        assert!(!bounds.contains(0.5));
        assert!(!bounds.contains(3.0));
    }

    #[test]
    fn reach_grows_with_sequence_separation() {
        let step = |steps: u32| f64::from(steps) * RESIDUE_STEP + 2.0 * SIDE_CHAIN_REACH;
        assert_eq!(step(0), 20.0);
        assert!(step(3) > step(1));
    }
}
