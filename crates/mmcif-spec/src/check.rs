use mmcif_ir::{Instantiated, IrError, Lengths, MAX_SLOTS};

use crate::diagnostic::{Diagnostic, Report};
use crate::geometry::{self, Endpoint};
use crate::lower::Lowered;
use crate::resolve;
use crate::span::{SourceId, Span};

pub fn check(source: SourceId, lowered: &Lowered, report: &mut Report) -> Option<Instantiated> {
    longest_fits(source, lowered, report)?;
    let shortest = instantiate(source, lowered, &Lengths::minimum(&lowered.spec), report)?;
    let ends = geometry::endpoints(lowered, &shortest);
    atoms(source, lowered, &ends, report);
    geometry::feasibility(source, lowered, &shortest, &ends, report);
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

fn atoms(source: SourceId, lowered: &Lowered, ends: &[Vec<Endpoint>], report: &mut Report) {
    for (constraint, pair) in ends.iter().enumerate() {
        for end in pair {
            let Some(key) = &end.anchor else { continue };
            let Some(target) = geometry::target(lowered, end.target) else {
                continue;
            };
            if resolve::atom_exists(target, key, &end.atom) != Some(false) {
                continue;
            }
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
                    span_of(lowered, Some(constraint)),
                    format!("residue {key} is {component}"),
                )
                .help(format!("atoms in the file: {}", names.join(", "))),
            );
        }
    }
}

pub fn span_of(lowered: &Lowered, constraint: Option<usize>) -> Span {
    constraint
        .and_then(|index| lowered.constraint_spans.get(index).copied())
        .or_else(|| lowered.segments.first().map(|segment| segment.span))
        .unwrap_or(Span::new(0, 0))
}
