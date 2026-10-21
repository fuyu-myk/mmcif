use std::collections::HashMap;

use mmcif_cif::Block;

use crate::error::ReadError;
use crate::ids::{
    AltLoc, AuthAsymId, AuthSeqId, CompId, Element, EntityId, InsCode, LabelAsymId, LabelSeqId,
    ModelNumber,
};
use crate::model::{AtomInput, ChainSlot, Model, ModelBuilder, ResidueSlot};
use crate::read::table::{Field, Table};

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

    fn atom(&self, row: usize) -> Result<AtomInput<'a>, ReadError> {
        let position = [
            self.x.required_parse(row)?,
            self.y.required_parse(row)?,
            self.z.required_parse(row)?,
        ];
        Ok(AtomInput {
            serial: self.serial.parse(row)?,
            name: self.preferred_text(&self.atom_name, &self.auth_atom_name, row)?,
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
struct ResidueGroup {
    label_seq: Option<LabelSeqId>,
    auth_seq: AuthSeqId,
    ins_code: InsCode,
    comp: CompId,
}

impl ResidueGroup {
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

struct ChainState {
    slot: ChainSlot,
    residues: HashMap<ResidueGroup, ResidueSlot>,
    last: Option<(ResidueGroup, ResidueSlot)>,
}

struct ModelState {
    builder: ModelBuilder,
    chains: Vec<ChainState>,
    chain_index: HashMap<Box<str>, usize>,
}

struct ResidueInput<'a> {
    label_seq: Option<LabelSeqId>,
    auth_seq: AuthSeqId,
    ins_code: InsCode,
    comp: &'a str,
    hetero: bool,
}

impl ModelState {
    fn new(number: i32) -> Self {
        ModelState {
            builder: ModelBuilder::new(ModelNumber::new(number)),
            chains: Vec::new(),
            chain_index: HashMap::new(),
        }
    }

    fn chain_for(
        &mut self,
        label_asym: &str,
        auth_asym: &str,
        entity: Option<&str>,
    ) -> Result<usize, ReadError> {
        let index = match self.chain_index.get(label_asym) {
            Some(&index) => index,
            None => {
                let slot = self.builder.new_chain(
                    LabelAsymId::new(label_asym),
                    AuthAsymId::new(auth_asym),
                    entity.map(EntityId::new),
                    Vec::new(),
                );
                self.chains.push(ChainState {
                    slot,
                    residues: HashMap::new(),
                    last: None,
                });
                self.chain_index
                    .insert(label_asym.into(), self.chains.len() - 1);
                self.chains.len() - 1
            }
        };
        let record = self.builder.chain_record(self.chains[index].slot);
        if record.auth_asym.as_str() != auth_asym {
            return Err(ReadError::InconsistentAuthChain {
                label_asym: record.label_asym.clone(),
                first: record.auth_asym.clone(),
                second: AuthAsymId::new(auth_asym),
            });
        }
        Ok(index)
    }

    fn residue_for(&mut self, chain: usize, input: ResidueInput<'_>) -> ResidueSlot {
        let state = &mut self.chains[chain];
        if let Some((last, slot)) = &state.last
            && last.matches(input.label_seq, input.auth_seq, input.ins_code, input.comp)
        {
            return *slot;
        }
        let group = ResidueGroup {
            label_seq: input.label_seq,
            auth_seq: input.auth_seq,
            ins_code: input.ins_code,
            comp: CompId::new(input.comp),
        };
        let builder = &mut self.builder;
        let slot = *state.residues.entry(group.clone()).or_insert_with(|| {
            builder.new_residue(
                state.slot,
                group.comp.clone(),
                group.label_seq,
                group.auth_seq,
                group.ins_code,
                input.hetero,
            )
        });
        state.last = Some((group, slot));
        slot
    }
}

pub(crate) struct AtomSite {
    pub models: Vec<Model>,
    pub has_author_numbering: bool,
}

pub(crate) fn read_models(block: &Block<'_>) -> Result<AtomSite, ReadError> {
    let table = Table::required(block, "atom_site")?;
    let fields = AtomSiteFields::new(&table)?;
    let rows = table.rows();
    let single_model =
        rows == 0 || fields.model.parse::<i32>(0)? == fields.model.parse::<i32>(rows - 1)?;
    let mut models: Vec<ModelState> = Vec::new();
    let mut model_index: HashMap<i32, usize> = HashMap::new();
    for row in 0..rows {
        let number = fields.model.parse::<i32>(row)?.unwrap_or(1);
        let index = *model_index.entry(number).or_insert_with(|| {
            let mut state = ModelState::new(number);
            if single_model && models.is_empty() {
                state.builder.reserve_atoms(rows);
            }
            models.push(state);
            models.len() - 1
        });
        let model = &mut models[index];
        let label_asym = fields.label_asym.required_text(row)?;
        let auth_asym = fields.auth_asym(row, label_asym)?;
        let chain = model.chain_for(label_asym, auth_asym, fields.entity.text(row))?;
        let label_seq = fields.label_seq(row)?;
        let residue = ResidueInput {
            label_seq,
            auth_seq: fields.auth_seq(row, label_seq)?,
            ins_code: fields
                .ins_code
                .single_char(row)?
                .map_or(InsCode::NONE, InsCode::new),
            comp: fields.preferred_text(&fields.comp, &fields.auth_comp, row)?,
            hetero: fields
                .group
                .text(row)
                .is_some_and(|group| group.eq_ignore_ascii_case("HETATM")),
        };
        let atom = fields.atom(row)?;
        let slot = model.residue_for(chain, residue);
        model.builder.push_atom(slot, atom);
    }
    Ok(AtomSite {
        models: models
            .into_iter()
            .map(|state| state.builder.finish())
            .collect(),
        has_author_numbering: fields.has_author_numbering(),
    })
}
