use crate::assembly::{Assembly, Operator};
use crate::entity::Entity;
use crate::ids::{
    AltLoc, AssemblyId, AtomName, AuthAsymId, AuthSeqId, CompId, Element, EntityId, InsCode,
    LabelAsymId, LabelSeqId, ModelNumber, OperatorId,
};
use crate::key::{AuthResidueKey, ResidueKey};
use crate::numbering::Numbering;

#[derive(Debug, Clone, PartialEq)]
pub struct Atom {
    pub(crate) serial: Option<i64>,
    pub(crate) name: AtomName,
    pub(crate) element: Element,
    pub(crate) altloc: AltLoc,
    pub(crate) position: [f64; 3],
    pub(crate) occupancy: f64,
    pub(crate) b_factor: f64,
    pub(crate) charge: Option<i8>,
}

impl Atom {
    pub fn new(name: AtomName, element: Element, position: [f64; 3]) -> Self {
        Atom {
            serial: None,
            name,
            element,
            altloc: AltLoc::NONE,
            position,
            occupancy: 1.0,
            b_factor: 0.0,
            charge: None,
        }
    }

    pub fn with_altloc(self, altloc: AltLoc) -> Self {
        Atom { altloc, ..self }
    }

    pub fn with_occupancy(self, occupancy: f64) -> Self {
        Atom { occupancy, ..self }
    }

    pub fn serial(&self) -> Option<i64> {
        self.serial
    }

    pub fn name(&self) -> &AtomName {
        &self.name
    }

    pub fn element(&self) -> Element {
        self.element
    }

    pub fn altloc(&self) -> AltLoc {
        self.altloc
    }

    pub fn position(&self) -> [f64; 3] {
        self.position
    }

    pub fn occupancy(&self) -> f64 {
        self.occupancy
    }

    pub fn b_factor(&self) -> f64 {
        self.b_factor
    }

    pub fn charge(&self) -> Option<i8> {
        self.charge
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Residue {
    pub(crate) comp: CompId,
    pub(crate) label_seq: Option<LabelSeqId>,
    pub(crate) auth_seq: AuthSeqId,
    pub(crate) ins_code: InsCode,
    pub(crate) hetero: bool,
    pub(crate) atoms: Vec<Atom>,
}

impl Residue {
    pub fn new(
        comp: CompId,
        label_seq: Option<LabelSeqId>,
        auth_seq: AuthSeqId,
        ins_code: InsCode,
    ) -> Self {
        Residue {
            comp,
            label_seq,
            auth_seq,
            ins_code,
            hetero: false,
            atoms: Vec::new(),
        }
    }

    pub fn with_atoms(self, atoms: Vec<Atom>) -> Self {
        Residue { atoms, ..self }
    }

    pub fn comp(&self) -> &CompId {
        &self.comp
    }

    pub fn label_seq(&self) -> Option<LabelSeqId> {
        self.label_seq
    }

    pub fn auth_seq(&self) -> AuthSeqId {
        self.auth_seq
    }

    pub fn ins_code(&self) -> InsCode {
        self.ins_code
    }

    pub fn is_hetero(&self) -> bool {
        self.hetero
    }

    pub fn atoms(&self) -> &[Atom] {
        &self.atoms
    }

    pub fn atom(&self, name: &str) -> Option<&Atom> {
        self.atoms.iter().find(|atom| atom.name.as_str() == name)
    }

    pub fn altlocs(&self) -> Vec<AltLoc> {
        let mut found = Vec::new();
        for atom in &self.atoms {
            if !atom.altloc.is_none() && !found.contains(&atom.altloc) {
                found.push(atom.altloc);
            }
        }
        found
    }

    pub fn conformer(&self, altloc: AltLoc) -> impl Iterator<Item = &Atom> {
        self.atoms
            .iter()
            .filter(move |atom| atom.altloc.is_none() || atom.altloc == altloc)
    }

    pub(crate) fn same_position(&self, other: &Residue) -> bool {
        match (self.label_seq, other.label_seq) {
            (Some(a), Some(b)) => a == b,
            _ => self.auth_seq == other.auth_seq && self.ins_code == other.ins_code,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Chain {
    pub(crate) label_asym: LabelAsymId,
    pub(crate) auth_asym: AuthAsymId,
    pub(crate) entity: Option<EntityId>,
    pub(crate) operators: Vec<OperatorId>,
    pub(crate) residues: Vec<Residue>,
}

impl Chain {
    pub fn new(label_asym: LabelAsymId, auth_asym: AuthAsymId, entity: Option<EntityId>) -> Self {
        Chain {
            label_asym,
            auth_asym,
            entity,
            operators: Vec::new(),
            residues: Vec::new(),
        }
    }

    pub fn with_residues(self, residues: Vec<Residue>) -> Self {
        Chain { residues, ..self }
    }

    pub fn label_asym(&self) -> &LabelAsymId {
        &self.label_asym
    }

    pub fn auth_asym(&self) -> &AuthAsymId {
        &self.auth_asym
    }

    pub fn entity(&self) -> Option<&EntityId> {
        self.entity.as_ref()
    }

    pub fn operators(&self) -> &[OperatorId] {
        &self.operators
    }

    pub fn residues(&self) -> &[Residue] {
        &self.residues
    }

    pub fn atom_count(&self) -> usize {
        self.residues
            .iter()
            .map(|residue| residue.atoms.len())
            .sum()
    }

    pub fn residues_at_label(&self, seq: LabelSeqId) -> impl Iterator<Item = &Residue> {
        self.residues
            .iter()
            .filter(move |residue| residue.label_seq == Some(seq))
    }

    pub fn residues_at_auth(
        &self,
        seq: AuthSeqId,
        ins_code: InsCode,
    ) -> impl Iterator<Item = &Residue> {
        self.residues
            .iter()
            .filter(move |residue| residue.auth_seq == seq && residue.ins_code == ins_code)
    }

    pub fn residue_key(&self, residue: &Residue) -> ResidueKey {
        ResidueKey {
            asym: self.label_asym.clone(),
            label_seq: residue.label_seq,
            auth_seq: residue.auth_seq,
            ins_code: residue.ins_code,
            comp: residue.comp.clone(),
        }
    }

    pub fn auth_key(&self, residue: &Residue) -> AuthResidueKey {
        AuthResidueKey {
            chain: self.auth_asym.clone(),
            seq: residue.auth_seq,
            ins_code: residue.ins_code,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Model {
    pub(crate) number: ModelNumber,
    pub(crate) chains: Vec<Chain>,
}

impl Model {
    pub fn new(number: ModelNumber, chains: Vec<Chain>) -> Self {
        Model { number, chains }
    }

    pub fn number(&self) -> ModelNumber {
        self.number
    }

    pub fn chains(&self) -> &[Chain] {
        &self.chains
    }

    pub fn chain(&self, label_asym: &str) -> Option<&Chain> {
        self.chains
            .iter()
            .find(|chain| chain.operators.is_empty() && chain.label_asym.as_str() == label_asym)
    }

    pub fn auth_chain<'a, 'k>(
        &'a self,
        auth_asym: &'k str,
    ) -> impl Iterator<Item = &'a Chain> + use<'a, 'k> {
        self.chains
            .iter()
            .filter(move |chain| chain.auth_asym.as_str() == auth_asym)
    }

    pub fn residues_at(&self, key: &AuthResidueKey) -> Vec<(&Chain, &Residue)> {
        self.auth_chain(key.chain.as_str())
            .flat_map(|chain| {
                chain
                    .residues_at_auth(key.seq, key.ins_code)
                    .map(move |residue| (chain, residue))
            })
            .collect()
    }

    pub fn residue(&self, key: &ResidueKey) -> Option<&Residue> {
        self.chain(key.asym.as_str())?
            .residues
            .iter()
            .find(|residue| {
                residue.label_seq == key.label_seq
                    && residue.auth_seq == key.auth_seq
                    && residue.ins_code == key.ins_code
                    && residue.comp == key.comp
            })
    }

    pub fn atom_count(&self) -> usize {
        self.chains.iter().map(Chain::atom_count).sum()
    }

    pub fn atoms(&self) -> impl Iterator<Item = (&Chain, &Residue, &Atom)> {
        self.chains.iter().flat_map(|chain| {
            chain.residues.iter().flat_map(move |residue| {
                residue.atoms.iter().map(move |atom| (chain, residue, atom))
            })
        })
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Structure {
    pub(crate) entry_id: Option<String>,
    pub(crate) entities: Vec<Entity>,
    pub(crate) models: Vec<Model>,
    pub(crate) operators: Vec<Operator>,
    pub(crate) assemblies: Vec<Assembly>,
    pub(crate) numbering: Numbering,
}

impl Structure {
    pub fn entry_id(&self) -> Option<&str> {
        self.entry_id.as_deref()
    }

    pub fn entities(&self) -> &[Entity] {
        &self.entities
    }

    pub fn entity(&self, id: &str) -> Option<&Entity> {
        self.entities
            .iter()
            .find(|entity| entity.id().as_str() == id)
    }

    pub fn models(&self) -> &[Model] {
        &self.models
    }

    pub fn model(&self, number: ModelNumber) -> Option<&Model> {
        self.models.iter().find(|model| model.number == number)
    }

    pub fn first_model(&self) -> Option<&Model> {
        self.models.first()
    }

    pub fn operators(&self) -> &[Operator] {
        &self.operators
    }

    pub fn operator(&self, id: &str) -> Option<&Operator> {
        self.operators
            .iter()
            .find(|operator| operator.id().as_str() == id)
    }

    pub fn assemblies(&self) -> &[Assembly] {
        &self.assemblies
    }

    pub fn assembly(&self, id: &AssemblyId) -> Option<&Assembly> {
        self.assemblies.iter().find(|assembly| assembly.id() == id)
    }

    pub fn numbering(&self) -> &Numbering {
        &self.numbering
    }
}
