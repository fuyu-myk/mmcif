use std::fmt;
use std::ops::Range;

use crate::ids::{
    AltLoc, AtomName, AuthAsymId, AuthSeqId, CompId, Element, EntityId, InsCode, LabelAsymId,
    LabelSeqId, OperatorId,
};
use crate::key::{AuthResidueKey, ResidueKey};
use crate::model::{AtomIndex, ChainIndex, ChainRecord, Model, ResidueIndex, ResidueRecord, span};

macro_rules! model_ref {
    ($name:ident) => {
        #[derive(Clone, Copy)]
        pub struct $name<'a> {
            model: &'a Model,
            index: u32,
        }

        impl<'a> $name<'a> {
            pub(crate) fn new(model: &'a Model, index: u32) -> Self {
                $name { model, index }
            }

            pub fn model(&self) -> &'a Model {
                self.model
            }
        }

        impl PartialEq for $name<'_> {
            fn eq(&self, other: &Self) -> bool {
                std::ptr::eq(self.model, other.model) && self.index == other.index
            }
        }

        impl Eq for $name<'_> {}
    };
}

model_ref!(ChainRef);
model_ref!(ResidueRef);
model_ref!(AtomRef);

impl<'a> ChainRef<'a> {
    fn record(&self) -> &'a ChainRecord {
        &self.model.chains[self.index as usize]
    }

    pub fn index(&self) -> ChainIndex {
        ChainIndex(self.index)
    }

    pub fn label_asym(&self) -> &'a LabelAsymId {
        &self.record().label_asym
    }

    pub fn auth_asym(&self) -> &'a AuthAsymId {
        &self.record().auth_asym
    }

    pub fn entity(&self) -> Option<&'a EntityId> {
        self.record().entity.as_ref()
    }

    pub fn operators(&self) -> &'a [OperatorId] {
        &self.record().operators
    }

    pub fn residue_count(&self) -> usize {
        self.record().residues.len()
    }

    pub fn residues(&self) -> impl ExactSizeIterator<Item = ResidueRef<'a>> + use<'a> {
        let model = self.model;
        self.record()
            .residues
            .clone()
            .map(move |index| ResidueRef::new(model, index))
    }

    pub fn atom_range(&self) -> Range<usize> {
        let residues = &self.record().residues;
        if residues.is_empty() {
            return 0..0;
        }
        let first = &self.model.residues[residues.start as usize];
        let last = &self.model.residues[residues.end as usize - 1];
        first.atoms.start as usize..last.atoms.end as usize
    }

    pub fn atom_count(&self) -> usize {
        self.atom_range().len()
    }

    pub fn atoms(&self) -> impl ExactSizeIterator<Item = AtomRef<'a>> + use<'a> {
        let model = self.model;
        self.atom_range()
            .map(move |index| AtomRef::new(model, index as u32))
    }

    pub fn positions(&self) -> &'a [[f64; 3]] {
        &self.model.atoms.positions[self.atom_range()]
    }

    pub fn residues_at_label(
        &self,
        seq: LabelSeqId,
    ) -> impl Iterator<Item = ResidueRef<'a>> + use<'a> {
        self.residues()
            .filter(move |residue| residue.label_seq() == Some(seq))
    }

    pub fn residues_at_auth(
        &self,
        seq: AuthSeqId,
        ins_code: InsCode,
    ) -> impl Iterator<Item = ResidueRef<'a>> + use<'a> {
        self.residues()
            .filter(move |residue| residue.auth_seq() == seq && residue.ins_code() == ins_code)
    }
}

impl fmt::Debug for ChainRef<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("Chain")
            .field("label_asym", self.label_asym())
            .field("auth_asym", self.auth_asym())
            .field("residues", &self.residue_count())
            .finish()
    }
}

impl<'a> ResidueRef<'a> {
    fn record(&self) -> &'a ResidueRecord {
        &self.model.residues[self.index as usize]
    }

    pub fn index(&self) -> ResidueIndex {
        ResidueIndex(self.index)
    }

    pub fn chain(&self) -> ChainRef<'a> {
        ChainRef::new(self.model, self.record().chain)
    }

    pub fn comp(&self) -> &'a CompId {
        &self.record().comp
    }

    pub fn label_seq(&self) -> Option<LabelSeqId> {
        self.record().label_seq
    }

    pub fn auth_seq(&self) -> AuthSeqId {
        self.record().auth_seq
    }

    pub fn ins_code(&self) -> InsCode {
        self.record().ins_code
    }

    pub fn is_hetero(&self) -> bool {
        self.record().hetero
    }

    pub fn atom_range(&self) -> Range<usize> {
        span(&self.record().atoms)
    }

    pub fn atom_count(&self) -> usize {
        self.record().atoms.len()
    }

    pub fn atoms(&self) -> impl ExactSizeIterator<Item = AtomRef<'a>> + use<'a> {
        let model = self.model;
        self.record()
            .atoms
            .clone()
            .map(move |index| AtomRef::new(model, index))
    }

    pub fn positions(&self) -> &'a [[f64; 3]] {
        &self.model.atoms.positions[self.atom_range()]
    }

    pub fn atom(&self, name: &str) -> Option<AtomRef<'a>> {
        self.atoms().find(|atom| atom.name().as_str() == name)
    }

    pub fn altlocs(&self) -> Vec<AltLoc> {
        let mut found = Vec::new();
        for &altloc in &self.model.atoms.altlocs[self.atom_range()] {
            if !altloc.is_none() && !found.contains(&altloc) {
                found.push(altloc);
            }
        }
        found
    }

    pub fn conformer(&self, altloc: AltLoc) -> impl Iterator<Item = AtomRef<'a>> + use<'a> {
        self.atoms()
            .filter(move |atom| atom.altloc().is_none() || atom.altloc() == altloc)
    }

    pub fn key(&self) -> ResidueKey {
        let record = self.record();
        ResidueKey {
            asym: self.chain().label_asym().clone(),
            label_seq: record.label_seq,
            auth_seq: record.auth_seq,
            ins_code: record.ins_code,
            comp: record.comp.clone(),
        }
    }

    pub fn auth_key(&self) -> AuthResidueKey {
        AuthResidueKey {
            chain: self.chain().auth_asym().clone(),
            seq: self.auth_seq(),
            ins_code: self.ins_code(),
        }
    }

    pub(crate) fn same_position(&self, other: &ResidueRef<'_>) -> bool {
        match (self.label_seq(), other.label_seq()) {
            (Some(a), Some(b)) => a == b,
            _ => self.auth_seq() == other.auth_seq() && self.ins_code() == other.ins_code(),
        }
    }
}

impl fmt::Debug for ResidueRef<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("Residue")
            .field("key", &self.key())
            .field("atoms", &self.atom_count())
            .finish()
    }
}

impl<'a> AtomRef<'a> {
    fn at(&self) -> usize {
        self.index as usize
    }

    pub fn index(&self) -> AtomIndex {
        AtomIndex(self.index)
    }

    pub fn residue(&self) -> ResidueRef<'a> {
        ResidueRef::new(self.model, self.model.atoms.residues[self.at()])
    }

    pub fn serial(&self) -> Option<i64> {
        self.model.atoms.serials[self.at()]
    }

    pub fn name(&self) -> &'a AtomName {
        self.model.atom_name(self.at())
    }

    pub fn element(&self) -> Element {
        self.model.atoms.elements[self.at()]
    }

    pub fn altloc(&self) -> AltLoc {
        self.model.atoms.altlocs[self.at()]
    }

    pub fn position(&self) -> [f64; 3] {
        self.model.atoms.positions[self.at()]
    }

    pub fn occupancy(&self) -> f64 {
        self.model.atoms.occupancies[self.at()]
    }

    pub fn b_factor(&self) -> f64 {
        self.model.atoms.b_factors[self.at()]
    }

    pub fn charge(&self) -> Option<i8> {
        self.model.atoms.charges[self.at()]
    }
}

impl fmt::Debug for AtomRef<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("Atom")
            .field("name", self.name())
            .field("altloc", &self.altloc())
            .field("position", &self.position())
            .finish()
    }
}
