use std::fmt;
use std::ops::Range;

use crate::ids::{
    AltLoc, AtomName, AuthAsymId, AuthSeqId, CompId, Element, EntityId, InsCode, LabelAsymId,
    LabelSeqId, OperatorId,
};
use crate::key::{AuthResidueKey, ResidueKey};
use crate::model::{AtomIndex, ChainIndex, ChainRecord, Model, ResidueIndex, ResidueRecord, span};

macro_rules! model_ref {
    ($(#[$meta:meta])* $name:ident) => {
        $(#[$meta])*
        #[derive(Clone, Copy)]
        pub struct $name<'a> {
            model: &'a Model,
            index: u32,
        }

        impl<'a> $name<'a> {
            pub(crate) fn new(model: &'a Model, index: u32) -> Self {
                $name { model, index }
            }

            /// Returns the model this view belongs to.
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

model_ref!(
    /// A view of one chain of a [`Model`].
    ///
    /// Two views are equal when they refer to the same chain of the same
    /// `Model` value.
    ///
    /// # Examples
    ///
    /// ```
    /// # let structure = mmcif_model::parse_mmcif(include_str!("../tests/fixtures/pathologies.cif"))?;
    /// let model = structure.first_model().unwrap();
    /// let chain = model.chain("AA").unwrap();
    ///
    /// assert_eq!(chain.auth_asym().as_str(), "H1");
    /// assert_eq!(chain.entity().unwrap().as_str(), "1");
    /// assert_eq!(chain.residue_count(), 5);
    /// assert_eq!(chain.atom_count(), 10);
    /// # Ok::<(), mmcif_model::ReadError>(())
    /// ```
    ChainRef
);
model_ref!(
    /// A view of one residue of a [`Model`].
    ///
    /// Two views are equal when they refer to the same residue of the same
    /// `Model` value.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_model::{InsCode, LabelSeqId};
    ///
    /// # let structure = mmcif_model::parse_mmcif(include_str!("../tests/fixtures/pathologies.cif"))?;
    /// let model = structure.first_model().unwrap();
    /// let alanine = model.chain("AA").unwrap().residues().nth(2).unwrap();
    ///
    /// assert_eq!(alanine.comp().as_str(), "ALA");
    /// assert_eq!(alanine.label_seq(), LabelSeqId::new(3));
    /// assert_eq!(alanine.auth_seq().get(), 100);
    /// assert_eq!(alanine.ins_code(), InsCode::new('A'));
    /// assert_eq!(alanine.atom_count(), 4);
    /// # Ok::<(), mmcif_model::ReadError>(())
    /// ```
    ResidueRef
);
model_ref!(
    /// A view of one atom of a [`Model`].
    ///
    /// Two views are equal when they refer to the same atom of the same
    /// `Model` value.
    ///
    /// # Examples
    ///
    /// ```
    /// # let structure = mmcif_model::parse_mmcif(include_str!("../tests/fixtures/pathologies.cif"))?;
    /// let model = structure.first_model().unwrap();
    /// let zinc = model.chain("C").unwrap().atoms().next().unwrap();
    ///
    /// assert_eq!(zinc.name().as_str(), "ZN");
    /// assert_eq!(zinc.element().symbol(), "Zn");
    /// assert_eq!(zinc.position(), [0.0, 0.0, 5.0]);
    /// assert_eq!(zinc.charge(), Some(2));
    /// assert_eq!(zinc.residue().auth_seq().get(), 201);
    /// # Ok::<(), mmcif_model::ReadError>(())
    /// ```
    AtomRef
);

impl<'a> ChainRef<'a> {
    fn record(&self) -> &'a ChainRecord {
        &self.model.chains[self.index as usize]
    }

    /// Returns the position of the chain in its model.
    pub fn index(&self) -> ChainIndex {
        ChainIndex(self.index)
    }

    /// Returns the label chain id.
    pub fn label_asym(&self) -> &'a LabelAsymId {
        &self.record().label_asym
    }

    /// Returns the author chain id.
    pub fn auth_asym(&self) -> &'a AuthAsymId {
        &self.record().auth_asym
    }

    /// Returns the id of the chain's entity, or `None` if the file does not
    /// give one.
    pub fn entity(&self) -> Option<&'a EntityId> {
        self.record().entity.as_ref()
    }

    /// Returns the ids of the operators applied to the chain, in the order
    /// they appear in the operator expression.
    ///
    /// The slice is empty for a chain as read from the file, and holds one id
    /// per group of the expression for a chain built by
    /// [`Structure::assemble`].
    ///
    /// [`Structure::assemble`]: crate::Structure::assemble
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_model::AssemblyId;
    ///
    /// # let structure = mmcif_model::parse_mmcif(include_str!("../tests/fixtures/pathologies.cif"))?;
    /// let model = structure.first_model().unwrap();
    /// assert!(model.chain("AA").unwrap().operators().is_empty());
    ///
    /// let assembled = structure.assemble(&AssemblyId::new("3"), model)?;
    /// let chain = assembled.chains().next().unwrap();
    /// let operators: Vec<&str> = chain.operators().iter().map(|id| id.as_str()).collect();
    /// assert_eq!(operators, ["2", "3"]);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn operators(&self) -> &'a [OperatorId] {
        &self.record().operators
    }

    /// Returns the number of residues in the chain.
    pub fn residue_count(&self) -> usize {
        self.record().residues.len()
    }

    /// Returns an iterator over the residues of the chain, in index order.
    pub fn residues(&self) -> impl ExactSizeIterator<Item = ResidueRef<'a>> + use<'a> {
        let model = self.model;
        self.record()
            .residues
            .clone()
            .map(move |index| ResidueRef::new(model, index))
    }

    /// Returns the range of atom indices the chain covers.
    ///
    /// The range is empty if the chain has no residues.
    pub fn atom_range(&self) -> Range<usize> {
        let residues = &self.record().residues;
        if residues.is_empty() {
            return 0..0;
        }
        let first = &self.model.residues[residues.start as usize];
        let last = &self.model.residues[residues.end as usize - 1];
        first.atoms.start as usize..last.atoms.end as usize
    }

    /// Returns the number of atoms in the chain.
    pub fn atom_count(&self) -> usize {
        self.atom_range().len()
    }

    /// Returns an iterator over the atoms of the chain, in index order.
    pub fn atoms(&self) -> impl ExactSizeIterator<Item = AtomRef<'a>> + use<'a> {
        let model = self.model;
        self.atom_range()
            .map(move |index| AtomRef::new(model, index as u32))
    }

    /// Returns the coordinates of the chain's atoms, in index order.
    pub fn positions(&self) -> &'a [[f64; 3]] {
        &self.model.atoms.positions[self.atom_range()]
    }

    /// Returns the residues of the chain with label residue number `seq`, in
    /// index order.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_model::LabelSeqId;
    ///
    /// # let structure = mmcif_model::parse_mmcif(include_str!("../tests/fixtures/pathologies.cif"))?;
    /// let chain = structure.first_model().unwrap().chain("AA").unwrap();
    ///
    /// let comps: Vec<&str> = chain
    ///     .residues_at_label(LabelSeqId::new(2).unwrap())
    ///     .map(|residue| residue.comp().as_str())
    ///     .collect();
    /// assert_eq!(comps, ["SER", "THR"]);
    /// # Ok::<(), mmcif_model::ReadError>(())
    /// ```
    pub fn residues_at_label(
        &self,
        seq: LabelSeqId,
    ) -> impl Iterator<Item = ResidueRef<'a>> + use<'a> {
        self.residues()
            .filter(move |residue| residue.label_seq() == Some(seq))
    }

    /// Returns the residues of the chain with author residue number `seq` and
    /// insertion code `ins_code`, in index order.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_model::{AuthSeqId, InsCode};
    ///
    /// # let structure = mmcif_model::parse_mmcif(include_str!("../tests/fixtures/pathologies.cif"))?;
    /// let chain = structure.first_model().unwrap().chain("AA").unwrap();
    ///
    /// let inserted = chain.residues_at_auth(AuthSeqId::new(100), InsCode::new('A')).next().unwrap();
    /// assert_eq!(inserted.comp().as_str(), "ALA");
    /// # Ok::<(), mmcif_model::ReadError>(())
    /// ```
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

    /// Returns the position of the residue in its model.
    pub fn index(&self) -> ResidueIndex {
        ResidueIndex(self.index)
    }

    /// Returns the chain the residue belongs to.
    pub fn chain(&self) -> ChainRef<'a> {
        ChainRef::new(self.model, self.record().chain)
    }

    /// Returns the chemical component id.
    pub fn comp(&self) -> &'a CompId {
        &self.record().comp
    }

    /// Returns the label residue number, or `None` if the residue has none.
    pub fn label_seq(&self) -> Option<LabelSeqId> {
        self.record().label_seq
    }

    /// Returns the author residue number.
    pub fn auth_seq(&self) -> AuthSeqId {
        self.record().auth_seq
    }

    /// Returns the insertion code.
    pub fn ins_code(&self) -> InsCode {
        self.record().ins_code
    }

    /// Returns `true` if the residue's first atom was read from a `HETATM`
    /// record.
    pub fn is_hetero(&self) -> bool {
        self.record().hetero
    }

    /// Returns the range of atom indices the residue covers.
    pub fn atom_range(&self) -> Range<usize> {
        span(&self.record().atoms)
    }

    /// Returns the number of atoms in the residue.
    pub fn atom_count(&self) -> usize {
        self.record().atoms.len()
    }

    /// Returns an iterator over the atoms of the residue, in index order.
    pub fn atoms(&self) -> impl ExactSizeIterator<Item = AtomRef<'a>> + use<'a> {
        let model = self.model;
        self.record()
            .atoms
            .clone()
            .map(move |index| AtomRef::new(model, index))
    }

    /// Returns the coordinates of the residue's atoms, in index order.
    pub fn positions(&self) -> &'a [[f64; 3]] {
        &self.model.atoms.positions[self.atom_range()]
    }

    /// Returns the first atom of the residue named `name`, or `None` if there
    /// is none.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_model::AltLoc;
    ///
    /// # let structure = mmcif_model::parse_mmcif(include_str!("../tests/fixtures/pathologies.cif"))?;
    /// let alanine = structure.first_model().unwrap().chain("AA").unwrap().residues().nth(2).unwrap();
    ///
    /// assert_eq!(alanine.atom("CB").unwrap().altloc(), AltLoc::new('A'));
    /// assert!(alanine.atom("OG").is_none());
    /// # Ok::<(), mmcif_model::ReadError>(())
    /// ```
    pub fn atom(&self, name: &str) -> Option<AtomRef<'a>> {
        self.atoms().find(|atom| atom.name().as_str() == name)
    }

    /// Returns the distinct alternate location indicators of the residue's
    /// atoms, in order of first appearance.
    ///
    /// [`AltLoc::NONE`] is not included.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_model::AltLoc;
    ///
    /// # let structure = mmcif_model::parse_mmcif(include_str!("../tests/fixtures/pathologies.cif"))?;
    /// let chain = structure.first_model().unwrap().chain("AA").unwrap();
    ///
    /// let alanine = chain.residues().nth(2).unwrap();
    /// assert_eq!(alanine.altlocs(), [AltLoc::new('A'), AltLoc::new('B')]);
    ///
    /// let lysine = chain.residues().nth(3).unwrap();
    /// assert!(lysine.altlocs().is_empty());
    /// # Ok::<(), mmcif_model::ReadError>(())
    /// ```
    pub fn altlocs(&self) -> Vec<AltLoc> {
        let mut found = Vec::new();
        for &altloc in &self.model.atoms.altlocs[self.atom_range()] {
            if !altloc.is_none() && !found.contains(&altloc) {
                found.push(altloc);
            }
        }
        found
    }

    /// Returns the atoms of conformer `altloc`: those with that indicator and
    /// those with no alternate location, in index order.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_model::AltLoc;
    ///
    /// # let structure = mmcif_model::parse_mmcif(include_str!("../tests/fixtures/pathologies.cif"))?;
    /// let alanine = structure.first_model().unwrap().chain("AA").unwrap().residues().nth(2).unwrap();
    ///
    /// let names: Vec<String> = alanine
    ///     .conformer(AltLoc::new('B'))
    ///     .map(|atom| format!("{}{}", atom.name(), atom.altloc()))
    ///     .collect();
    /// assert_eq!(names, ["N", "CA", "CBB"]);
    /// # Ok::<(), mmcif_model::ReadError>(())
    /// ```
    pub fn conformer(&self, altloc: AltLoc) -> impl Iterator<Item = AtomRef<'a>> + use<'a> {
        self.atoms()
            .filter(move |atom| atom.altloc().is_none() || atom.altloc() == altloc)
    }

    /// Returns the full address of the residue.
    ///
    /// # Examples
    ///
    /// ```
    /// # let structure = mmcif_model::parse_mmcif(include_str!("../tests/fixtures/pathologies.cif"))?;
    /// let model = structure.first_model().unwrap();
    /// let zinc = model.chain("C").unwrap().residues().next().unwrap();
    ///
    /// assert_eq!(zinc.key().to_string(), "C/./201/ZN");
    /// assert_eq!(model.residue(&zinc.key()), Some(zinc));
    /// # Ok::<(), mmcif_model::ReadError>(())
    /// ```
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

    /// Returns the author address of the residue.
    ///
    /// # Examples
    ///
    /// ```
    /// # let structure = mmcif_model::parse_mmcif(include_str!("../tests/fixtures/pathologies.cif"))?;
    /// let alanine = structure.first_model().unwrap().chain("AA").unwrap().residues().nth(2).unwrap();
    ///
    /// assert_eq!(alanine.auth_key().to_string(), "H1:100A");
    /// # Ok::<(), mmcif_model::ReadError>(())
    /// ```
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

    /// Returns the position of the atom in its model.
    pub fn index(&self) -> AtomIndex {
        AtomIndex(self.index)
    }

    /// Returns the residue the atom belongs to.
    pub fn residue(&self) -> ResidueRef<'a> {
        ResidueRef::new(self.model, self.model.atoms.residues[self.at()])
    }

    /// Returns the atom serial number, or `None` if the file gives none.
    pub fn serial(&self) -> Option<i64> {
        self.model.atoms.serials[self.at()]
    }

    /// Returns the atom name.
    pub fn name(&self) -> &'a AtomName {
        self.model.atom_name(self.at())
    }

    /// Returns the element.
    pub fn element(&self) -> Element {
        self.model.atoms.elements[self.at()]
    }

    /// Returns the alternate location indicator.
    pub fn altloc(&self) -> AltLoc {
        self.model.atoms.altlocs[self.at()]
    }

    /// Returns the Cartesian coordinates.
    pub fn position(&self) -> [f64; 3] {
        self.model.atoms.positions[self.at()]
    }

    /// Returns the occupancy, or 1.0 if the file gives none.
    pub fn occupancy(&self) -> f64 {
        self.model.atoms.occupancies[self.at()]
    }

    /// Returns the isotropic B-factor, or 0.0 if the file gives none.
    pub fn b_factor(&self) -> f64 {
        self.model.atoms.b_factors[self.at()]
    }

    /// Returns the formal charge, or `None` if the file gives none.
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
