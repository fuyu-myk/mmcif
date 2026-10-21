use std::path::{Path, PathBuf};

use mmcif_model::{LabelResidueKey, LabelSeqId, ResidueRef, Structure};

use crate::diagnostic::{Diagnostic, Report};
use crate::document::{Ranges, Spanned};
use crate::span::SourceId;
#[cfg(test)]
use crate::span::Span;

pub struct Target {
    pub path: PathBuf,
    pub structure: Structure,
}

pub fn load(
    source: SourceId,
    base: &Path,
    path: &Spanned<String>,
    report: &mut Report,
) -> Option<Target> {
    let full = base.join(&path.value);
    match mmcif_model::read(&full) {
        Ok(structure) => Some(Target {
            path: full,
            structure,
        }),
        Err(error) => {
            report.push(
                Diagnostic::error(
                    "unreadable-target",
                    format!("{} could not be read", path.value),
                )
                .at(source, path.span, error.to_string())
                .help(format!("looked in {}", full.display())),
            );
            None
        }
    }
}

pub struct Residue {
    pub seq: LabelSeqId,
    pub key: LabelResidueKey,
    pub comp: mmcif_model::CompId,
}

pub fn chain(
    source: SourceId,
    target: &Target,
    chain: &Spanned<String>,
    report: &mut Report,
) -> Option<Vec<Residue>> {
    let model = target.structure.first_model().or_else(|| {
        report.push(
            Diagnostic::error("empty-target", "the target holds no model").at(
                source,
                chain.span,
                "nothing to resolve against",
            ),
        );
        None
    })?;
    let Some(found) = model.chain(&chain.value) else {
        let mut names: Vec<&str> = model
            .chains()
            .filter(|chain| chain.operators().is_empty())
            .map(|chain| chain.label_asym().as_str())
            .collect();
        names.dedup();
        report.push(
            Diagnostic::error(
                "unknown-chain",
                format!("{} has no chain {}", display(&target.path), chain.value),
            )
            .at(source, chain.span, "not in the target")
            .help(format!("chains in the file: {}", names.join(", "))),
        );
        return None;
    };
    Some(
        found
            .residues()
            .filter_map(|residue| polymer(&residue))
            .collect(),
    )
}

fn polymer(residue: &ResidueRef<'_>) -> Option<Residue> {
    let seq = residue.label_seq()?;
    Some(Residue {
        seq,
        key: LabelResidueKey {
            asym: residue.chain().label_asym().clone(),
            seq,
        },
        comp: residue.comp().clone(),
    })
}

pub fn select(
    source: SourceId,
    target: &Target,
    chain: &Spanned<String>,
    residues: Vec<Residue>,
    ranges: Option<&Spanned<Ranges>>,
    report: &mut Report,
) -> Vec<Residue> {
    let Some(ranges) = ranges else {
        return residues;
    };
    let selected: Vec<Residue> = residues
        .into_iter()
        .filter(|residue| ranges.value.contains(residue.seq.get()))
        .collect();
    if selected.is_empty() {
        report.push(
            Diagnostic::error(
                "empty-selection",
                format!("no residue of chain {} is selected", chain.value),
            )
            .at(source, ranges.span, "selects nothing")
            .help(available(source, target, chain).unwrap_or_else(|| "the chain is empty".into())),
        );
    }
    selected
}

fn available(source: SourceId, target: &Target, chain: &Spanned<String>) -> Option<String> {
    let mut report = Report::new();
    let residues = self::chain(source, target, chain, &mut report)?;
    let first = residues.first()?.seq.get();
    let last = residues.last()?.seq.get();
    Some(format!(
        "chain {} covers residues {first} to {last}",
        chain.value
    ))
}

pub fn atom_exists(target: &Target, key: &LabelResidueKey, atom: &str) -> Option<bool> {
    let residue = residue(target, key)?;
    Some(residue.atom(atom).is_some())
}

pub fn position(target: &Target, key: &LabelResidueKey, atom: &str) -> Option<[f64; 3]> {
    Some(residue(target, key)?.atom(atom)?.position())
}

pub fn atom_names(target: &Target, key: &LabelResidueKey) -> Vec<String> {
    match residue(target, key) {
        None => Vec::new(),
        Some(residue) => residue
            .atoms()
            .map(|atom| atom.name().as_str().to_string())
            .collect(),
    }
}

pub fn component(target: &Target, key: &LabelResidueKey) -> Option<mmcif_model::CompId> {
    residue(target, key).map(|residue| residue.comp().clone())
}

fn residue<'a>(target: &'a Target, key: &LabelResidueKey) -> Option<ResidueRef<'a>> {
    target
        .structure
        .first_model()?
        .chain(key.asym.as_str())?
        .residues()
        .find(|residue| residue.label_seq() == Some(key.seq))
}

pub fn display(path: &Path) -> String {
    path.file_name()
        .map(|name| name.to_string_lossy().into_owned())
        .unwrap_or_else(|| path.display().to_string())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::span::Sources;

    const CIF: &str = "data_test
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.auth_seq_id
_atom_site.auth_asym_id
ATOM 1 N N CYS A 1 1 0.000 0.000 0.000 1 A
ATOM 2 C CA CYS A 1 1 1.500 0.000 0.000 1 A
ATOM 3 S SG CYS A 1 1 2.000 1.000 0.000 1 A
ATOM 4 N N GLY A 1 2 3.800 0.000 0.000 2 A
ATOM 5 C CA GLY A 1 2 5.000 0.000 0.000 2 A
ATOM 6 N N CYS A 1 3 7.600 0.000 0.000 3 A
ATOM 7 C CA CYS A 1 3 8.000 0.000 0.000 3 A
ATOM 8 S SG CYS A 1 3 4.000 1.000 0.000 3 A
ATOM 9 N N ALA B 2 1 20.000 0.000 0.000 1 B
ATOM 10 C CA ALA B 2 1 21.000 0.000 0.000 1 B
";

    fn target() -> Target {
        Target {
            path: PathBuf::from("/tmp/test.cif"),
            structure: mmcif_model::parse_mmcif(CIF).expect("the fixture parses"),
        }
    }

    fn named(value: &str) -> Spanned<String> {
        Spanned::new(Span::new(0, value.len()), value.to_string())
    }

    fn source() -> SourceId {
        Sources::new().add("design.yaml", "")
    }

    fn key(seq: u32) -> LabelResidueKey {
        LabelResidueKey {
            asym: "A".into(),
            seq: LabelSeqId::new(seq).unwrap(),
        }
    }

    #[test]
    fn a_chain_resolves_to_its_polymer_residues() {
        let mut report = Report::new();
        let residues = chain(source(), &target(), &named("A"), &mut report).unwrap();
        assert!(report.is_empty());
        assert_eq!(residues.len(), 3);
        assert_eq!(residues[0].comp.as_str(), "CYS");
        assert_eq!(residues[2].seq.get(), 3);
        assert_eq!(residues[2].key, key(3));
    }

    #[test]
    fn an_unknown_chain_lists_the_ones_in_the_file() {
        let mut report = Report::new();
        assert!(chain(source(), &target(), &named("Z"), &mut report).is_none());
        assert_eq!(report.codes(), ["unknown-chain"]);
        assert_eq!(
            report.diagnostics()[0].help.as_deref(),
            Some("chains in the file: A, B")
        );
    }

    #[test]
    fn a_range_keeps_only_the_residues_it_covers() {
        let target = target();
        let mut report = Report::new();
        let residues = chain(source(), &target, &named("A"), &mut report).unwrap();
        let ranges = Spanned::new(
            Span::new(0, 4),
            Ranges::List(vec![crate::document::Range {
                min: 2,
                max: Some(3),
            }]),
        );
        let selected = select(
            source(),
            &target,
            &named("A"),
            residues,
            Some(&ranges),
            &mut report,
        );
        assert!(report.is_empty());
        assert_eq!(selected.len(), 2);
        assert_eq!(selected[0].seq.get(), 2);
    }

    #[test]
    fn a_range_that_selects_nothing_reports_the_available_span() {
        let target = target();
        let mut report = Report::new();
        let residues = chain(source(), &target, &named("A"), &mut report).unwrap();
        let ranges = Spanned::new(
            Span::new(0, 6),
            Ranges::List(vec![crate::document::Range {
                min: 90,
                max: Some(95),
            }]),
        );
        let selected = select(
            source(),
            &target,
            &named("A"),
            residues,
            Some(&ranges),
            &mut report,
        );
        assert!(selected.is_empty());
        assert_eq!(report.codes(), ["empty-selection"]);
        assert_eq!(
            report.diagnostics()[0].help.as_deref(),
            Some("chain A covers residues 1 to 3")
        );
    }

    #[test]
    fn atoms_and_positions_come_from_the_file() {
        let target = target();
        assert_eq!(atom_exists(&target, &key(1), "SG"), Some(true));
        assert_eq!(atom_exists(&target, &key(2), "SG"), Some(false));
        assert_eq!(atom_exists(&target, &key(9), "SG"), None);
        assert_eq!(position(&target, &key(1), "SG"), Some([2.0, 1.0, 0.0]));
        assert_eq!(position(&target, &key(1), "ZZ"), None);
        assert_eq!(atom_names(&target, &key(2)), ["N", "CA"]);
        assert!(atom_names(&target, &key(9)).is_empty());
        assert_eq!(component(&target, &key(3)).unwrap().as_str(), "CYS");
    }

    #[test]
    fn an_unreadable_target_names_the_path_it_tried() {
        let mut report = Report::new();
        assert!(
            load(
                source(),
                Path::new("/tmp/nowhere"),
                &named("missing.cif"),
                &mut report
            )
            .is_none()
        );
        assert_eq!(report.codes(), ["unreadable-target"]);
        assert_eq!(
            report.diagnostics()[0].help.as_deref(),
            Some("looked in /tmp/nowhere/missing.cif")
        );
    }

    #[test]
    fn file_names_are_shown_without_their_directory() {
        assert_eq!(display(Path::new("/tmp/designs/8r3a.cif")), "8r3a.cif");
    }
}
