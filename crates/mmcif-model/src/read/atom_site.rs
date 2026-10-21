use std::collections::HashMap;

use mmcif_cif::Block;

use crate::error::ReadError;
use crate::ids::{
    AltLoc, AtomName, AuthAsymId, AuthSeqId, CompId, Element, EntityId, InsCode, LabelAsymId,
    LabelSeqId, ModelNumber,
};
use crate::read::table::{Field, Table};
use crate::structure::{Atom, Chain, Model, Residue};

struct AtomSiteFields<'a> {
    group: Field<'a>,
    serial: Field<'a>,
    element: Field<'a>,
    atom_name: Field<'a>,
    auth_atom_name: Field<'a>,
    altloc: Field<'a>,
    comp: Field<'a>,
    auth_comp: Field<'a>,
    label_asym: Field<'a>,
    auth_asym: Field<'a>,
    entity: Field<'a>,
    label_seq: Field<'a>,
    auth_seq: Field<'a>,
    ins_code: Field<'a>,
    x: Field<'a>,
    y: Field<'a>,
    z: Field<'a>,
    occupancy: Field<'a>,
    b_factor: Field<'a>,
    charge: Field<'a>,
    model: Field<'a>,
}

impl<'a> AtomSiteFields<'a> {
    fn new(table: &Table<'a>) -> Result<Self, ReadError> {
        Ok(AtomSiteFields {
            group: table.column("group_PDB"),
            serial: table.column("id"),
            element: table.required_column("type_symbol")?,
            atom_name: table.column("label_atom_id"),
            auth_atom_name: table.column("auth_atom_id"),
            altloc: table.column("label_alt_id"),
            comp: table.column("label_comp_id"),
            auth_comp: table.column("auth_comp_id"),
            label_asym: table.required_column("label_asym_id")?,
            auth_asym: table.column("auth_asym_id"),
            entity: table.column("label_entity_id"),
            label_seq: table.column("label_seq_id"),
            auth_seq: table.column("auth_seq_id"),
            ins_code: table.column("pdbx_PDB_ins_code"),
            x: table.required_column("Cartn_x")?,
            y: table.required_column("Cartn_y")?,
            z: table.required_column("Cartn_z")?,
            occupancy: table.column("occupancy"),
            b_factor: table.column("B_iso_or_equiv"),
            charge: table.column("pdbx_formal_charge"),
            model: table.column("pdbx_PDB_model_num"),
        })
    }

    fn preferred_text(
        &self,
        primary: &Field<'a>,
        fallback: &Field<'a>,
        row: usize,
    ) -> Result<&'a str, ReadError> {
        match primary.text(row).or_else(|| fallback.text(row)) {
            Some(text) => Ok(text),
            None => primary.required_text(row),
        }
    }

    fn label_seq(&self, row: usize) -> Result<Option<LabelSeqId>, ReadError> {
        let Some(value) = self.label_seq.parse::<u32>(row)? else {
            return Ok(None);
        };
        LabelSeqId::new(value)
            .map(Some)
            .ok_or_else(|| self.label_seq.invalid(row, "0", "label_seq_id starts at 1"))
    }

    fn auth_seq(&self, row: usize, label_seq: Option<LabelSeqId>) -> Result<AuthSeqId, ReadError> {
        if self.auth_seq.is_present() {
            return self.auth_seq.required_parse(row).map(AuthSeqId::new);
        }
        match label_seq {
            Some(label) => i32::try_from(label.get()).map(AuthSeqId::new).map_err(|_| {
                self.label_seq
                    .invalid(row, &label.to_string(), "out of range")
            }),
            None => Err(self.auth_seq.invalid(row, "", "missing value")),
        }
    }

    fn auth_asym(&self, row: usize, label_asym: &'a str) -> Result<&'a str, ReadError> {
        if self.auth_asym.is_present() {
            self.auth_asym.required_text(row)
        } else {
            Ok(label_asym)
        }
    }

    fn has_author_numbering(&self) -> bool {
        self.auth_asym.is_present() && self.auth_seq.is_present()
    }

    fn element(&self, row: usize) -> Result<Element, ReadError> {
        let text = self.element.required_text(row)?;
        Element::parse(text).ok_or_else(|| self.element.invalid(row, text, "not an element symbol"))
    }

    fn atom(&self, row: usize) -> Result<Atom, ReadError> {
        let position = [
            self.x.required_parse(row)?,
            self.y.required_parse(row)?,
            self.z.required_parse(row)?,
        ];
        Ok(Atom {
            serial: self.serial.parse(row)?,
            name: AtomName::new(self.preferred_text(&self.atom_name, &self.auth_atom_name, row)?),
            element: self.element(row)?,
            altloc: self
                .altloc
                .single_char(row)?
                .map_or(AltLoc::NONE, AltLoc::new),
            position,
            occupancy: self.occupancy.parse(row)?.unwrap_or(1.0),
            b_factor: self.b_factor.parse(row)?.unwrap_or(0.0),
            charge: self.charge.parse(row)?,
        })
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
struct ResidueSlot {
    label_seq: Option<LabelSeqId>,
    auth_seq: AuthSeqId,
    ins_code: InsCode,
    comp: CompId,
}

struct ChainBuilder {
    chain: Chain,
    slots: HashMap<ResidueSlot, usize>,
    last_slot: Option<(ResidueSlot, usize)>,
}

impl ResidueSlot {
    fn matches(
        &self,
        label_seq: Option<LabelSeqId>,
        auth_seq: AuthSeqId,
        ins_code: InsCode,
        comp: &str,
    ) -> bool {
        self.label_seq == label_seq
            && self.auth_seq == auth_seq
            && self.ins_code == ins_code
            && self.comp.as_str() == comp
    }
}

struct ModelBuilder {
    number: ModelNumber,
    chains: Vec<ChainBuilder>,
    chain_index: HashMap<Box<str>, usize>,
}

impl ModelBuilder {
    fn chain_for(
        &mut self,
        label_asym: &str,
        auth_asym: &str,
        entity: Option<&str>,
    ) -> Result<&mut ChainBuilder, ReadError> {
        let index = match self.chain_index.get(label_asym) {
            Some(&index) => index,
            None => {
                let chain = Chain::new(
                    LabelAsymId::new(label_asym),
                    AuthAsymId::new(auth_asym),
                    entity.map(EntityId::new),
                );
                self.chains.push(ChainBuilder {
                    chain,
                    slots: HashMap::new(),
                    last_slot: None,
                });
                self.chain_index
                    .insert(label_asym.into(), self.chains.len() - 1);
                self.chains.len() - 1
            }
        };
        let builder = &mut self.chains[index];
        let chain = &builder.chain;
        if chain.auth_asym.as_str() != auth_asym {
            return Err(ReadError::InconsistentAuthChain {
                label_asym: chain.label_asym.clone(),
                first: chain.auth_asym.clone(),
                second: AuthAsymId::new(auth_asym),
            });
        }
        Ok(builder)
    }

    fn finish(self) -> Model {
        let chains = self
            .chains
            .into_iter()
            .map(|builder| builder.chain)
            .collect();
        Model::new(self.number, chains)
    }
}

impl ChainBuilder {
    fn residue_for(
        &mut self,
        label_seq: Option<LabelSeqId>,
        auth_seq: AuthSeqId,
        ins_code: InsCode,
        comp: &str,
        hetero: bool,
    ) -> &mut Residue {
        let index = match &self.last_slot {
            Some((last, index)) if last.matches(label_seq, auth_seq, ins_code, comp) => *index,
            _ => {
                let slot = ResidueSlot {
                    label_seq,
                    auth_seq,
                    ins_code,
                    comp: CompId::new(comp),
                };
                let residues = &mut self.chain.residues;
                let index = *self.slots.entry(slot.clone()).or_insert_with(|| {
                    let mut residue =
                        Residue::new(slot.comp.clone(), label_seq, auth_seq, ins_code);
                    residue.hetero = hetero;
                    residues.push(residue);
                    residues.len() - 1
                });
                self.last_slot = Some((slot, index));
                index
            }
        };
        &mut self.chain.residues[index]
    }
}

pub(crate) struct AtomSite {
    pub models: Vec<Model>,
    pub has_author_numbering: bool,
}

pub(crate) fn read_models(block: &Block<'_>) -> Result<AtomSite, ReadError> {
    let table = Table::required(block, "atom_site")?;
    let fields = AtomSiteFields::new(&table)?;
    let mut models: Vec<ModelBuilder> = Vec::new();
    let mut model_index: HashMap<i32, usize> = HashMap::new();
    for row in 0..table.rows() {
        let number = fields.model.parse::<i32>(row)?.unwrap_or(1);
        let index = *model_index.entry(number).or_insert_with(|| {
            models.push(ModelBuilder {
                number: ModelNumber::new(number),
                chains: Vec::new(),
                chain_index: HashMap::new(),
            });
            models.len() - 1
        });
        let label_asym = fields.label_asym.required_text(row)?;
        let auth_asym = fields.auth_asym(row, label_asym)?;
        let chain = models[index].chain_for(label_asym, auth_asym, fields.entity.text(row))?;
        let label_seq = fields.label_seq(row)?;
        let auth_seq = fields.auth_seq(row, label_seq)?;
        let ins_code = fields
            .ins_code
            .single_char(row)?
            .map_or(InsCode::NONE, InsCode::new);
        let comp = fields.preferred_text(&fields.comp, &fields.auth_comp, row)?;
        let hetero = fields
            .group
            .text(row)
            .is_some_and(|group| group.eq_ignore_ascii_case("HETATM"));
        let atom = fields.atom(row)?;
        chain
            .residue_for(label_seq, auth_seq, ins_code, comp, hetero)
            .atoms
            .push(atom);
    }
    Ok(AtomSite {
        models: models.into_iter().map(ModelBuilder::finish).collect(),
        has_author_numbering: fields.has_author_numbering(),
    })
}
