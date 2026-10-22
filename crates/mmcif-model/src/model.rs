use std::collections::{HashMap, TryReserveError};
use std::ops::Range;

use crate::ids::{
    AltLoc, AtomName, AuthAsymId, AuthSeqId, CompId, Element, EntityId, InsCode, LabelAsymId,
    LabelSeqId, ModelNumber, OperatorId,
};
use crate::key::{AuthResidueKey, ResidueKey};
use crate::view::{AtomRef, ChainRef, ResidueRef};

/// The position of a chain within its [`Model`], from [`ChainRef::index`].
///
/// # Examples
///
/// ```
/// # let structure = mmcif_model::parse_mmcif(include_str!("../tests/fixtures/pathologies.cif"))?;
/// let model = structure.first_model().unwrap();
/// let index = model.chain("B").unwrap().index();
///
/// assert_eq!(index.get(), 1);
/// assert_eq!(model.chain_at(index).unwrap().label_asym().as_str(), "B");
/// # Ok::<(), mmcif_model::ReadError>(())
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct ChainIndex(pub(crate) u32);

/// The position of a residue within its [`Model`], from
/// [`ResidueRef::index`].
///
/// # Examples
///
/// ```
/// # let structure = mmcif_model::parse_mmcif(include_str!("../tests/fixtures/pathologies.cif"))?;
/// let model = structure.first_model().unwrap();
/// let index = model.chain("C").unwrap().residues().next().unwrap().index();
///
/// assert_eq!(index.get(), 7);
/// assert_eq!(model.residue_at(index).unwrap().comp().as_str(), "ZN");
/// # Ok::<(), mmcif_model::ReadError>(())
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct ResidueIndex(pub(crate) u32);

/// The position of an atom within its [`Model`], from [`AtomRef::index`].
///
/// Atom indices also index the per-atom columns of a model, such as
/// [`Model::positions`].
///
/// # Examples
///
/// ```
/// # let structure = mmcif_model::parse_mmcif(include_str!("../tests/fixtures/pathologies.cif"))?;
/// let model = structure.first_model().unwrap();
/// let zinc = model.chain("C").unwrap().atoms().next().unwrap();
///
/// assert_eq!(zinc.index().get(), 13);
/// assert_eq!(model.positions()[zinc.index().get()], [0.0, 0.0, 5.0]);
/// # Ok::<(), mmcif_model::ReadError>(())
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct AtomIndex(pub(crate) u32);

macro_rules! index_accessors {
    ($($name:ident),*) => {
        $(impl $name {
            /// Returns the index as a `usize`.
            pub fn get(self) -> usize {
                self.0 as usize
            }
        })*
    };
}

index_accessors!(ChainIndex, ResidueIndex, AtomIndex);

pub(crate) fn span(range: &Range<u32>) -> Range<usize> {
    range.start as usize..range.end as usize
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct ChainRecord {
    pub label_asym: LabelAsymId,
    pub auth_asym: AuthAsymId,
    pub entity: Option<EntityId>,
    pub operators: Vec<OperatorId>,
    pub residues: Range<u32>,
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct ResidueRecord {
    pub comp: CompId,
    pub label_seq: Option<LabelSeqId>,
    pub auth_seq: AuthSeqId,
    pub ins_code: InsCode,
    pub hetero: bool,
    pub chain: u32,
    pub atoms: Range<u32>,
}

#[derive(Debug, Clone, PartialEq, Default)]
pub(crate) struct AtomColumns {
    pub names: Vec<u32>,
    pub name_table: Vec<AtomName>,
    pub elements: Vec<Element>,
    pub altlocs: Vec<AltLoc>,
    pub serials: Vec<Option<i64>>,
    pub positions: Vec<[f64; 3]>,
    pub occupancies: Vec<f64>,
    pub b_factors: Vec<f64>,
    pub charges: Vec<Option<i8>>,
    pub residues: Vec<u32>,
}

/// One model of a structure, holding its chains, residues and atoms.
///
/// Residues are grouped by chain and atoms by residue: each chain and each
/// residue covers a contiguous range of atom indices. Per-atom values are
/// also available as columns indexed by atom index, such as
/// [`Model::positions`].
///
/// # Examples
///
/// ```
/// # let structure = mmcif_model::parse_mmcif(include_str!("../tests/fixtures/pathologies.cif"))?;
/// let model = structure.first_model().unwrap();
///
/// assert_eq!(model.number().get(), 1);
/// assert_eq!(model.chain_count(), 4);
/// assert_eq!(model.residue_count(), 10);
/// assert_eq!(model.atom_count(), 16);
///
/// let zinc = model.chain("C").unwrap();
/// assert_eq!(zinc.atom_range(), 13..14);
/// assert_eq!(model.b_factors()[13], 30.0);
/// # Ok::<(), mmcif_model::ReadError>(())
/// ```
#[derive(Debug, Clone, PartialEq)]
pub struct Model {
    pub(crate) number: ModelNumber,
    pub(crate) chains: Vec<ChainRecord>,
    pub(crate) residues: Vec<ResidueRecord>,
    pub(crate) atoms: AtomColumns,
}

impl Model {
    /// Returns the model number.
    pub fn number(&self) -> ModelNumber {
        self.number
    }

    /// Returns the number of chains.
    pub fn chain_count(&self) -> usize {
        self.chains.len()
    }

    /// Returns the number of residues.
    pub fn residue_count(&self) -> usize {
        self.residues.len()
    }

    /// Returns the number of atoms.
    pub fn atom_count(&self) -> usize {
        self.atoms.positions.len()
    }

    /// Returns an iterator over the chains, in index order.
    ///
    /// # Examples
    ///
    /// ```
    /// # let structure = mmcif_model::parse_mmcif(include_str!("../tests/fixtures/pathologies.cif"))?;
    /// let model = structure.first_model().unwrap();
    ///
    /// let labels: Vec<&str> = model.chains().map(|chain| chain.label_asym().as_str()).collect();
    /// assert_eq!(labels, ["AA", "B", "C", "D"]);
    /// # Ok::<(), mmcif_model::ReadError>(())
    /// ```
    pub fn chains(&self) -> impl ExactSizeIterator<Item = ChainRef<'_>> {
        (0..self.chains.len() as u32).map(move |index| ChainRef::new(self, index))
    }

    /// Returns an iterator over every residue, in index order.
    pub fn residues(&self) -> impl ExactSizeIterator<Item = ResidueRef<'_>> {
        (0..self.residues.len() as u32).map(move |index| ResidueRef::new(self, index))
    }

    /// Returns an iterator over every atom, in index order.
    pub fn atoms(&self) -> impl ExactSizeIterator<Item = AtomRef<'_>> {
        (0..self.atom_count() as u32).map(move |index| AtomRef::new(self, index))
    }

    /// Returns the chain at `index`, or `None` if `index` is out of range.
    pub fn chain_at(&self, index: ChainIndex) -> Option<ChainRef<'_>> {
        (index.get() < self.chains.len()).then(|| ChainRef::new(self, index.0))
    }

    /// Returns the residue at `index`, or `None` if `index` is out of range.
    pub fn residue_at(&self, index: ResidueIndex) -> Option<ResidueRef<'_>> {
        (index.get() < self.residues.len()).then(|| ResidueRef::new(self, index.0))
    }

    /// Returns the atom at `index`, or `None` if `index` is out of range.
    pub fn atom_at(&self, index: AtomIndex) -> Option<AtomRef<'_>> {
        (index.get() < self.atom_count()).then(|| AtomRef::new(self, index.0))
    }

    /// Returns the chain with label id `label_asym` that has no operators
    /// applied, or `None` if there is none.
    ///
    /// Every chain of a model built by [`Structure::assemble`] has operators
    /// applied.
    ///
    /// [`Structure::assemble`]: crate::Structure::assemble
    ///
    /// # Examples
    ///
    /// ```
    /// # let structure = mmcif_model::parse_mmcif(include_str!("../tests/fixtures/pathologies.cif"))?;
    /// let model = structure.first_model().unwrap();
    ///
    /// let waters = model.chain("D").unwrap();
    /// assert_eq!(waters.residue_count(), 2);
    /// assert!(model.chain("H1").is_none());
    /// # Ok::<(), mmcif_model::ReadError>(())
    /// ```
    pub fn chain(&self, label_asym: &str) -> Option<ChainRef<'_>> {
        self.chains()
            .find(|chain| chain.operators().is_empty() && chain.label_asym().as_str() == label_asym)
    }

    /// Returns the chains with author chain id `auth_asym`, in index order.
    ///
    /// # Examples
    ///
    /// ```
    /// # let structure = mmcif_model::parse_mmcif(include_str!("../tests/fixtures/pathologies.cif"))?;
    /// let model = structure.first_model().unwrap();
    ///
    /// let labels: Vec<&str> = model.auth_chain("H1").map(|chain| chain.label_asym().as_str()).collect();
    /// assert_eq!(labels, ["AA", "B", "C", "D"]);
    /// # Ok::<(), mmcif_model::ReadError>(())
    /// ```
    pub fn auth_chain<'a, 'k>(
        &'a self,
        auth_asym: &'k str,
    ) -> impl Iterator<Item = ChainRef<'a>> + use<'a, 'k> {
        self.chains()
            .filter(move |chain| chain.auth_asym().as_str() == auth_asym)
    }

    /// Returns the residues at author address `key`, in index order.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_model::{AuthResidueKey, AuthSeqId, InsCode};
    ///
    /// # let structure = mmcif_model::parse_mmcif(include_str!("../tests/fixtures/pathologies.cif"))?;
    /// let model = structure.first_model().unwrap();
    /// let key = AuthResidueKey {
    ///     chain: "H1".into(),
    ///     seq: AuthSeqId::new(100),
    ///     ins_code: InsCode::NONE,
    /// };
    ///
    /// let comps: Vec<&str> = model.residues_at(&key).iter().map(|residue| residue.comp().as_str()).collect();
    /// assert_eq!(comps, ["SER", "THR"]);
    /// # Ok::<(), mmcif_model::ReadError>(())
    /// ```
    pub fn residues_at(&self, key: &AuthResidueKey) -> Vec<ResidueRef<'_>> {
        self.auth_chain(key.chain.as_str())
            .flat_map(|chain| chain.residues_at_auth(key.seq, key.ins_code))
            .collect()
    }

    /// Returns the residue that matches every field of `key` in the chain
    /// [`Model::chain`] returns for `key.asym`, or `None` if there is none.
    ///
    /// # Examples
    ///
    /// ```
    /// # let structure = mmcif_model::parse_mmcif(include_str!("../tests/fixtures/pathologies.cif"))?;
    /// let model = structure.first_model().unwrap();
    ///
    /// for residue in model.residues() {
    ///     assert_eq!(model.residue(&residue.key()), Some(residue));
    /// }
    /// # Ok::<(), mmcif_model::ReadError>(())
    /// ```
    pub fn residue(&self, key: &ResidueKey) -> Option<ResidueRef<'_>> {
        self.chain(key.asym.as_str())?.residues().find(|residue| {
            residue.label_seq() == key.label_seq
                && residue.auth_seq() == key.auth_seq
                && residue.ins_code() == key.ins_code
                && residue.comp() == &key.comp
        })
    }

    /// Returns the Cartesian coordinates of every atom, by atom index.
    pub fn positions(&self) -> &[[f64; 3]] {
        &self.atoms.positions
    }

    /// Returns the element of every atom, by atom index.
    pub fn elements(&self) -> &[Element] {
        &self.atoms.elements
    }

    /// Returns the alternate location indicator of every atom, by atom index.
    pub fn altlocs(&self) -> &[AltLoc] {
        &self.atoms.altlocs
    }

    /// Returns the occupancy of every atom, by atom index.
    ///
    /// An atom whose file gives no occupancy has 1.0.
    pub fn occupancies(&self) -> &[f64] {
        &self.atoms.occupancies
    }

    /// Returns the isotropic B-factor of every atom, by atom index.
    ///
    /// An atom whose file gives no B-factor has 0.0.
    pub fn b_factors(&self) -> &[f64] {
        &self.atoms.b_factors
    }

    /// Returns the formal charge of every atom, by atom index, with `None`
    /// where the file gives none.
    pub fn charges(&self) -> &[Option<i8>] {
        &self.atoms.charges
    }

    /// Returns the serial number of every atom, by atom index, with `None`
    /// where the file gives none.
    pub fn serials(&self) -> &[Option<i64>] {
        &self.atoms.serials
    }

    /// Returns the residue index of every atom, by atom index.
    pub fn atom_residue_indices(&self) -> &[u32] {
        &self.atoms.residues
    }

    pub(crate) fn atom_name(&self, atom: usize) -> &AtomName {
        &self.atoms.name_table[self.atoms.names[atom] as usize]
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct ChainSlot(u32);

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct ResidueSlot(u32);

pub(crate) struct AtomInput<'a> {
    pub serial: Option<i64>,
    pub name: &'a str,
    pub element: Element,
    pub altloc: AltLoc,
    pub position: [f64; 3],
    pub occupancy: f64,
    pub b_factor: f64,
    pub charge: Option<i8>,
}

pub(crate) struct ModelBuilder {
    number: ModelNumber,
    chains: Vec<ChainRecord>,
    chain_residues: Vec<Vec<u32>>,
    residues: Vec<ResidueRecord>,
    atoms: AtomColumns,
    name_ids: HashMap<Box<str>, u32>,
}

fn to_u32(value: usize) -> u32 {
    u32::try_from(value).expect("models hold fewer than 2^32 chains, residues and atoms")
}

impl ModelBuilder {
    pub fn new(number: ModelNumber) -> Self {
        ModelBuilder {
            number,
            chains: Vec::new(),
            chain_residues: Vec::new(),
            residues: Vec::new(),
            atoms: AtomColumns::default(),
            name_ids: HashMap::new(),
        }
    }

    pub fn derived_from(source: &Model) -> Self {
        let mut builder = ModelBuilder::new(source.number);
        builder.atoms.name_table = source.atoms.name_table.clone();
        builder.name_ids = source
            .atoms
            .name_table
            .iter()
            .enumerate()
            .map(|(index, name)| (Box::from(name.as_str()), to_u32(index)))
            .collect();
        builder
    }

    pub fn try_reserve_atoms(&mut self, additional: usize) -> Result<(), TryReserveError> {
        let atoms = &mut self.atoms;
        atoms.names.try_reserve(additional)?;
        atoms.elements.try_reserve(additional)?;
        atoms.altlocs.try_reserve(additional)?;
        atoms.serials.try_reserve(additional)?;
        atoms.positions.try_reserve(additional)?;
        atoms.occupancies.try_reserve(additional)?;
        atoms.b_factors.try_reserve(additional)?;
        atoms.charges.try_reserve(additional)?;
        atoms.residues.try_reserve(additional)
    }

    pub fn reserve_atoms(&mut self, additional: usize) {
        let atoms = &mut self.atoms;
        atoms.names.reserve(additional);
        atoms.elements.reserve(additional);
        atoms.altlocs.reserve(additional);
        atoms.serials.reserve(additional);
        atoms.positions.reserve(additional);
        atoms.occupancies.reserve(additional);
        atoms.b_factors.reserve(additional);
        atoms.charges.reserve(additional);
        atoms.residues.reserve(additional);
    }

    pub fn chain_record(&mut self, slot: ChainSlot) -> &mut ChainRecord {
        &mut self.chains[slot.0 as usize]
    }

    pub fn new_chain(
        &mut self,
        label_asym: LabelAsymId,
        auth_asym: AuthAsymId,
        entity: Option<EntityId>,
        operators: Vec<OperatorId>,
    ) -> ChainSlot {
        self.chains.push(ChainRecord {
            label_asym,
            auth_asym,
            entity,
            operators,
            residues: 0..0,
        });
        self.chain_residues.push(Vec::new());
        ChainSlot(to_u32(self.chains.len() - 1))
    }

    pub fn new_residue(
        &mut self,
        chain: ChainSlot,
        comp: CompId,
        label_seq: Option<LabelSeqId>,
        auth_seq: AuthSeqId,
        ins_code: InsCode,
        hetero: bool,
    ) -> ResidueSlot {
        let slot = to_u32(self.residues.len());
        self.residues.push(ResidueRecord {
            comp,
            label_seq,
            auth_seq,
            ins_code,
            hetero,
            chain: chain.0,
            atoms: 0..0,
        });
        self.chain_residues[chain.0 as usize].push(slot);
        ResidueSlot(slot)
    }

    fn name_id(&mut self, name: &str) -> u32 {
        if let Some(&id) = self.name_ids.get(name) {
            return id;
        }
        let id = to_u32(self.atoms.name_table.len());
        self.atoms.name_table.push(AtomName::new(name));
        self.name_ids.insert(Box::from(name), id);
        id
    }

    pub fn push_atom(&mut self, residue: ResidueSlot, atom: AtomInput<'_>) {
        let name = self.name_id(atom.name);
        let atoms = &mut self.atoms;
        atoms.names.push(name);
        atoms.elements.push(atom.element);
        atoms.altlocs.push(atom.altloc);
        atoms.serials.push(atom.serial);
        atoms.positions.push(atom.position);
        atoms.occupancies.push(atom.occupancy);
        atoms.b_factors.push(atom.b_factor);
        atoms.charges.push(atom.charge);
        atoms.residues.push(residue.0);
    }

    pub fn push_copied_atom(
        &mut self,
        residue: ResidueSlot,
        source: &Model,
        atom: usize,
        position: [f64; 3],
    ) {
        let from = &source.atoms;
        let atoms = &mut self.atoms;
        atoms.names.push(from.names[atom]);
        atoms.elements.push(from.elements[atom]);
        atoms.altlocs.push(from.altlocs[atom]);
        atoms.serials.push(from.serials[atom]);
        atoms.positions.push(position);
        atoms.occupancies.push(from.occupancies[atom]);
        atoms.b_factors.push(from.b_factors[atom]);
        atoms.charges.push(from.charges[atom]);
        atoms.residues.push(residue.0);
    }

    pub fn finish(self) -> Model {
        let ModelBuilder {
            number,
            mut chains,
            chain_residues,
            residues,
            mut atoms,
            ..
        } = self;
        let mut final_index = vec![0u32; residues.len()];
        let mut order = Vec::with_capacity(residues.len());
        for (chain, slots) in chains.iter_mut().zip(&chain_residues) {
            let start = to_u32(order.len());
            for &slot in slots {
                final_index[slot as usize] = to_u32(order.len());
                order.push(slot);
            }
            chain.residues = start..to_u32(order.len());
        }
        for residue in &mut atoms.residues {
            *residue = final_index[*residue as usize];
        }
        sort_atoms_by_residue(&mut atoms, order.len());
        compact_names(&mut atoms);
        let mut slots: Vec<Option<ResidueRecord>> = residues.into_iter().map(Some).collect();
        let mut residues: Vec<ResidueRecord> = order
            .iter()
            .map(|&slot| {
                slots[slot as usize]
                    .take()
                    .expect("each residue slot belongs to exactly one chain")
            })
            .collect();
        let mut start = 0usize;
        for (index, residue) in residues.iter_mut().enumerate() {
            let end = start
                + atoms.residues[start..]
                    .iter()
                    .take_while(|&&r| r as usize == index)
                    .count();
            residue.atoms = to_u32(start)..to_u32(end);
            start = end;
        }
        Model {
            number,
            chains,
            residues,
            atoms,
        }
    }
}

fn compact_names(atoms: &mut AtomColumns) {
    let mut remap = vec![u32::MAX; atoms.name_table.len()];
    let mut table = Vec::new();
    for id in &mut atoms.names {
        let slot = &mut remap[*id as usize];
        if *slot == u32::MAX {
            *slot = to_u32(table.len());
            table.push(atoms.name_table[*id as usize].clone());
        }
        *id = *slot;
    }
    atoms.name_table = table;
}

fn sort_atoms_by_residue(atoms: &mut AtomColumns, residue_count: usize) {
    if atoms.residues.is_sorted() {
        return;
    }
    let mut starts = vec![0usize; residue_count + 1];
    for &residue in &atoms.residues {
        starts[residue as usize + 1] += 1;
    }
    for index in 1..starts.len() {
        starts[index] += starts[index - 1];
    }
    let mut destination = vec![0usize; atoms.residues.len()];
    for (atom, &residue) in atoms.residues.iter().enumerate() {
        destination[atom] = starts[residue as usize];
        starts[residue as usize] += 1;
    }
    fn permute<T: Copy>(column: &mut Vec<T>, destination: &[usize]) {
        let mut sorted = column.clone();
        for (from, &to) in destination.iter().enumerate() {
            sorted[to] = column[from];
        }
        *column = sorted;
    }
    permute(&mut atoms.names, &destination);
    permute(&mut atoms.elements, &destination);
    permute(&mut atoms.altlocs, &destination);
    permute(&mut atoms.serials, &destination);
    permute(&mut atoms.positions, &destination);
    permute(&mut atoms.occupancies, &destination);
    permute(&mut atoms.b_factors, &destination);
    permute(&mut atoms.charges, &destination);
    permute(&mut atoms.residues, &destination);
}
