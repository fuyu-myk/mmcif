use crate::ids::AltLoc;
use crate::model::{Model, ModelBuilder};
use crate::view::{ChainRef, ResidueRef};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AltlocPolicy {
    First,
    HighestOccupancy,
    Prefer(AltLoc),
}

impl Model {
    pub fn single_conformer(&self, policy: AltlocPolicy) -> Model {
        let mut builder = ModelBuilder::derived_from(self);
        builder.reserve_atoms(self.atom_count());
        for chain in self.chains() {
            copy_single_conformer_chain(self, &mut builder, chain, policy);
        }
        builder.finish()
    }
}

fn copy_single_conformer_chain(
    model: &Model,
    builder: &mut ModelBuilder,
    chain: ChainRef<'_>,
    policy: AltlocPolicy,
) {
    let slot = builder.new_chain(
        chain.label_asym().clone(),
        chain.auth_asym().clone(),
        chain.entity().cloned(),
        chain.operators().to_vec(),
    );
    let residues: Vec<ResidueRef<'_>> = chain.residues().collect();
    let mut rest = residues.as_slice();
    while let Some(first) = rest.first() {
        let group_len = rest
            .iter()
            .take_while(|residue| residue.same_position(first))
            .count();
        let (group, tail) = rest.split_at(group_len);
        let chosen = choose_altloc(group, policy);
        for residue in group {
            let kept: Vec<usize> = residue
                .atoms()
                .filter(|atom| atom.altloc().is_none() || Some(atom.altloc()) == chosen)
                .map(|atom| atom.index().get())
                .collect();
            if kept.is_empty() {
                continue;
            }
            let residue_slot = builder.new_residue(
                slot,
                residue.comp().clone(),
                residue.label_seq(),
                residue.auth_seq(),
                residue.ins_code(),
                residue.is_hetero(),
            );
            for atom in kept {
                builder.push_copied_atom(residue_slot, model, atom, model.positions()[atom]);
            }
        }
        rest = tail;
    }
}

fn choose_altloc(group: &[ResidueRef<'_>], policy: AltlocPolicy) -> Option<AltLoc> {
    let mut totals: Vec<(AltLoc, f64, usize)> = Vec::new();
    for atom in group.iter().flat_map(|residue| residue.atoms()) {
        let altloc = atom.altloc();
        if altloc.is_none() {
            continue;
        }
        match totals.iter_mut().find(|(seen, _, _)| *seen == altloc) {
            Some((_, occupancy, count)) => {
                *occupancy += atom.occupancy();
                *count += 1;
            }
            None => totals.push((altloc, atom.occupancy(), 1)),
        }
    }
    let candidates: Vec<(AltLoc, f64)> = totals
        .into_iter()
        .map(|(altloc, occupancy, count)| (altloc, occupancy / count as f64))
        .collect();
    let first = candidates.first()?.0;
    match policy {
        AltlocPolicy::First => Some(first),
        AltlocPolicy::Prefer(wanted) => {
            if candidates.iter().any(|(altloc, _)| *altloc == wanted) {
                Some(wanted)
            } else {
                Some(first)
            }
        }
        AltlocPolicy::HighestOccupancy => candidates
            .iter()
            .fold(
                None::<(AltLoc, f64)>,
                |best, &(altloc, occupancy)| match best {
                    Some((_, best_occupancy)) if best_occupancy >= occupancy => best,
                    _ => Some((altloc, occupancy)),
                },
            )
            .map(|(altloc, _)| altloc),
    }
}
