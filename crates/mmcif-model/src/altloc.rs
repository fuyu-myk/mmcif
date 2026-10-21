use crate::ids::AltLoc;
use crate::structure::{Chain, Model, Residue};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AltlocPolicy {
    First,
    HighestOccupancy,
    Prefer(AltLoc),
}

impl Model {
    pub fn single_conformer(&self, policy: AltlocPolicy) -> Model {
        let chains = self
            .chains
            .iter()
            .map(|chain| single_conformer_chain(chain, policy))
            .collect();
        Model::new(self.number, chains)
    }
}

fn single_conformer_chain(chain: &Chain, policy: AltlocPolicy) -> Chain {
    let mut residues = Vec::with_capacity(chain.residues.len());
    let mut rest = chain.residues.as_slice();
    while let Some(first) = rest.first() {
        let group_len = rest
            .iter()
            .take_while(|residue| residue.same_position(first))
            .count();
        let (group, tail) = rest.split_at(group_len);
        let chosen = choose_altloc(group, policy);
        residues.extend(
            group
                .iter()
                .filter_map(|residue| keep_altloc(residue, chosen)),
        );
        rest = tail;
    }
    Chain {
        label_asym: chain.label_asym.clone(),
        auth_asym: chain.auth_asym.clone(),
        entity: chain.entity.clone(),
        operators: chain.operators.clone(),
        residues,
    }
}

fn choose_altloc(group: &[Residue], policy: AltlocPolicy) -> Option<AltLoc> {
    let mut totals: Vec<(AltLoc, f64, usize)> = Vec::new();
    for atom in group.iter().flat_map(|residue| &residue.atoms) {
        if atom.altloc.is_none() {
            continue;
        }
        match totals
            .iter_mut()
            .find(|(altloc, _, _)| *altloc == atom.altloc)
        {
            Some((_, occupancy, count)) => {
                *occupancy += atom.occupancy;
                *count += 1;
            }
            None => totals.push((atom.altloc, atom.occupancy, 1)),
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

fn keep_altloc(residue: &Residue, chosen: Option<AltLoc>) -> Option<Residue> {
    let atoms: Vec<_> = residue
        .atoms
        .iter()
        .filter(|atom| atom.altloc.is_none() || Some(atom.altloc) == chosen)
        .cloned()
        .collect();
    (!atoms.is_empty()).then(|| Residue {
        comp: residue.comp.clone(),
        label_seq: residue.label_seq,
        auth_seq: residue.auth_seq,
        ins_code: residue.ins_code,
        hetero: residue.hetero,
        atoms,
    })
}
