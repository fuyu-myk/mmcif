use std::collections::HashMap;

use mmcif_cif::Block;

use crate::entity::{
    Branch, BranchAtom, BranchLink, BranchMonomer, Entity, EntityKind, Polymer, PolymerType,
    SequencePosition,
};
use crate::error::ReadError;
use crate::ids::{AtomName, CompId, EntityId, LabelAsymId, LabelSeqId};
use crate::read::table::Table;
use crate::structure::Model;

pub(crate) fn read_entities(block: &Block<'_>, models: &[Model]) -> Result<Vec<Entity>, ReadError> {
    let Some(table) = Table::optional(block, "entity") else {
        return Ok(Vec::new());
    };
    let ids = table.required_column("id")?;
    let kinds = table.column("type");
    let descriptions = table.column("pdbx_description");
    let mut polymers = read_polymers(block)?;
    let mut branches = read_branches(block)?;
    let asym_ids = read_asym_ids(block, models)?;
    (0..table.rows())
        .map(|row| {
            let id = ids.required_text(row)?;
            Ok(Entity {
                id: EntityId::new(id),
                kind: kinds
                    .text(row)
                    .map_or(EntityKind::Other(String::new()), EntityKind::parse),
                description: descriptions.text(row).map(str::to_owned),
                polymer: polymers.remove(id),
                branch: branches.remove(id),
                asym_ids: asym_ids.get(id).cloned().unwrap_or_default(),
            })
        })
        .collect()
}

fn read_polymers(block: &Block<'_>) -> Result<HashMap<String, Polymer>, ReadError> {
    let mut polymers = HashMap::new();
    if let Some(table) = Table::optional(block, "entity_poly") {
        let ids = table.required_column("entity_id")?;
        let kinds = table.column("type");
        let canonical = table.column("pdbx_seq_one_letter_code_can");
        for row in 0..table.rows() {
            let polymer = Polymer {
                kind: kinds
                    .text(row)
                    .map_or(PolymerType::Other(String::new()), PolymerType::parse),
                canonical_sequence: canonical
                    .text(row)
                    .map(|text| text.chars().filter(|c| !c.is_whitespace()).collect()),
                sequence: Vec::new(),
            };
            polymers.insert(ids.required_text(row)?.to_owned(), polymer);
        }
    }
    if let Some(table) = Table::optional(block, "entity_poly_seq") {
        let ids = table.required_column("entity_id")?;
        let nums = table.required_column("num")?;
        let comps = table.required_column("mon_id")?;
        for row in 0..table.rows() {
            let num: u32 = nums.required_parse(row)?;
            let seq = LabelSeqId::new(num)
                .ok_or_else(|| nums.invalid(row, "0", "sequence numbers start at 1"))?;
            let comp = CompId::new(comps.required_text(row)?);
            let polymer = polymers
                .entry(ids.required_text(row)?.to_owned())
                .or_insert_with(|| Polymer {
                    kind: PolymerType::Other(String::new()),
                    canonical_sequence: None,
                    sequence: Vec::new(),
                });
            match polymer
                .sequence
                .iter_mut()
                .rev()
                .find(|position| position.seq == seq)
            {
                Some(position) if !position.comps.contains(&comp) => position.comps.push(comp),
                Some(_) => {}
                None => polymer.sequence.push(SequencePosition {
                    seq,
                    comps: vec![comp],
                }),
            }
        }
    }
    Ok(polymers)
}

fn read_branches(block: &Block<'_>) -> Result<HashMap<String, Branch>, ReadError> {
    let mut branches: HashMap<String, Branch> = HashMap::new();
    if let Some(table) = Table::optional(block, "pdbx_entity_branch") {
        let ids = table.required_column("entity_id")?;
        let kinds = table.column("type");
        for row in 0..table.rows() {
            branch_entry(&mut branches, ids.required_text(row)?).kind =
                kinds.text(row).map(str::to_owned);
        }
    }
    if let Some(table) = Table::optional(block, "pdbx_entity_branch_list") {
        let ids = table.required_column("entity_id")?;
        let nums = table.required_column("num")?;
        let comps = table.required_column("comp_id")?;
        for row in 0..table.rows() {
            let monomer = BranchMonomer {
                num: nums.required_parse(row)?,
                comp: CompId::new(comps.required_text(row)?),
            };
            branch_entry(&mut branches, ids.required_text(row)?)
                .monomers
                .push(monomer);
        }
    }
    if let Some(table) = Table::optional(block, "pdbx_entity_branch_link") {
        let ids = table.required_column("entity_id")?;
        let order = table.column("value_order");
        let side = |suffix: u8| -> Result<_, ReadError> {
            let (num, comp, atom, leaving) = match suffix {
                1 => (
                    "entity_branch_list_num_1",
                    "comp_id_1",
                    "atom_id_1",
                    "leaving_atom_id_1",
                ),
                _ => (
                    "entity_branch_list_num_2",
                    "comp_id_2",
                    "atom_id_2",
                    "leaving_atom_id_2",
                ),
            };
            Ok((
                table.required_column(num)?,
                table.required_column(comp)?,
                table.required_column(atom)?,
                table.column(leaving),
            ))
        };
        let sides = [side(1)?, side(2)?];
        for row in 0..table.rows() {
            let [first, second] = sides.each_ref().map(|(num, comp, atom, leaving)| {
                Ok::<_, ReadError>(BranchAtom {
                    num: num.required_parse(row)?,
                    comp: CompId::new(comp.required_text(row)?),
                    atom: AtomName::new(atom.required_text(row)?),
                    leaving_atom: leaving.text(row).map(AtomName::new),
                })
            });
            let link = BranchLink {
                first: first?,
                second: second?,
                order: order.text(row).map(str::to_owned),
            };
            branch_entry(&mut branches, ids.required_text(row)?)
                .links
                .push(link);
        }
    }
    Ok(branches)
}

fn branch_entry<'m>(branches: &'m mut HashMap<String, Branch>, id: &str) -> &'m mut Branch {
    branches.entry(id.to_owned()).or_insert_with(|| Branch {
        kind: None,
        monomers: Vec::new(),
        links: Vec::new(),
    })
}

fn read_asym_ids(
    block: &Block<'_>,
    models: &[Model],
) -> Result<HashMap<String, Vec<LabelAsymId>>, ReadError> {
    let mut asym_ids: HashMap<String, Vec<LabelAsymId>> = HashMap::new();
    let mut add = |entity: &str, asym: &LabelAsymId| {
        let list = asym_ids.entry(entity.to_owned()).or_default();
        if !list.contains(asym) {
            list.push(asym.clone());
        }
    };
    if let Some(table) = Table::optional(block, "struct_asym") {
        let ids = table.required_column("id")?;
        let entities = table.required_column("entity_id")?;
        for row in 0..table.rows() {
            if let Some(entity) = entities.text(row) {
                add(entity, &LabelAsymId::new(ids.required_text(row)?));
            }
        }
    }
    for chain in models.iter().flat_map(|model| &model.chains) {
        if let Some(entity) = &chain.entity {
            add(entity.as_str(), &chain.label_asym);
        }
    }
    for list in asym_ids.values_mut() {
        list.sort();
    }
    Ok(asym_ids)
}
