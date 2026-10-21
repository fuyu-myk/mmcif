mod assemblies;
mod atom_site;
mod entities;
mod numbering;
mod table;

use std::collections::HashMap;
use std::path::Path;

use mmcif_cif::{Block, Document};

use crate::entity::Entity;
use crate::error::ReadError;
use crate::ids::{EntityId, LabelAsymId};
use crate::pdb;
use crate::structure::{Model, Structure};

impl Structure {
    pub fn from_document(document: &Document) -> Result<Structure, ReadError> {
        Structure::from_block(&document.block(0).ok_or(ReadError::NoDataBlock)?)
    }

    pub fn from_block(block: &Block<'_>) -> Result<Structure, ReadError> {
        let atom_site::AtomSite {
            mut models,
            has_author_numbering,
        } = atom_site::read_models(block)?;
        let entry_id = block
            .category("entry")
            .and_then(|entry| entry.column("id"))
            .and_then(|column| column.get(0))
            .and_then(|value| value.as_str())
            .map(str::to_owned);
        let entities = entities::read_entities(block, &models)?;
        assign_missing_entities(&mut models, &entities);
        Ok(Structure {
            entry_id,
            numbering: numbering::read_numbering(block, &models, has_author_numbering)?,
            entities,
            models,
            operators: assemblies::read_operators(block)?,
            assemblies: assemblies::read_assemblies(block)?,
        })
    }
}

fn assign_missing_entities(models: &mut [Model], entities: &[Entity]) {
    let by_asym: HashMap<&LabelAsymId, &EntityId> = entities
        .iter()
        .flat_map(|entity| entity.asym_ids.iter().map(move |asym| (asym, &entity.id)))
        .collect();
    for chain in models.iter_mut().flat_map(|model| model.chains.iter_mut()) {
        if chain.entity.is_none() {
            chain.entity = by_asym.get(&chain.label_asym).map(|&id| id.clone());
        }
    }
}

pub fn parse_mmcif(text: &str) -> Result<Structure, ReadError> {
    Structure::from_document(&Document::parse(text)?)
}

pub fn read_mmcif(path: impl AsRef<Path>) -> Result<Structure, ReadError> {
    Structure::from_document(&Document::open(path)?)
}

pub fn parse_bcif(bytes: &[u8]) -> Result<Structure, ReadError> {
    Structure::from_document(&Document::from_bcif(bytes).map_err(mmcif_cif::Error::from)?)
}

pub fn read_bcif(path: impl AsRef<Path>) -> Result<Structure, ReadError> {
    Structure::from_document(&Document::open_bcif(path)?)
}

pub fn read(path: impl AsRef<Path>) -> Result<Structure, ReadError> {
    let path = path.as_ref();
    let extension = path
        .extension()
        .and_then(|extension| extension.to_str())
        .map(str::to_ascii_lowercase);
    match extension.as_deref() {
        Some("cif" | "mmcif") => read_mmcif(path),
        Some("bcif") => read_bcif(path),
        Some("pdb" | "ent") => pdb::read_pdb(path),
        _ => Err(ReadError::UnknownFormat(path.display().to_string())),
    }
}
