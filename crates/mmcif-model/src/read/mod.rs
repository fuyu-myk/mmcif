mod assemblies;
mod atom_site;
mod entities;
mod numbering;
mod table;

use std::path::Path;

use mmcif_cif::{Block, Document};

use crate::error::ReadError;
use crate::pdb;
use crate::structure::{Model, Structure};

impl Structure {
    pub fn from_document(document: &Document) -> Result<Structure, ReadError> {
        Structure::from_block(&document.block(0).ok_or(ReadError::NoDataBlock)?)
    }

    pub fn from_block(block: &Block<'_>) -> Result<Structure, ReadError> {
        let mut models = atom_site::read_models(block)?;
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
            numbering: numbering::read_numbering(block, &models)?,
            entities,
            models,
            operators: assemblies::read_operators(block)?,
            assemblies: assemblies::read_assemblies(block)?,
        })
    }
}

fn assign_missing_entities(models: &mut [Model], entities: &[crate::entity::Entity]) {
    for chain in models.iter_mut().flat_map(|model| model.chains.iter_mut()) {
        if chain.entity.is_none() {
            chain.entity = entities
                .iter()
                .find(|entity| entity.asym_ids.contains(&chain.label_asym))
                .map(|entity| entity.id.clone());
        }
    }
}

pub fn parse_mmcif(text: &str) -> Result<Structure, ReadError> {
    Structure::from_document(&Document::parse(text)?)
}

pub fn read_mmcif(path: impl AsRef<Path>) -> Result<Structure, ReadError> {
    Structure::from_document(&Document::open(path)?)
}

pub fn read(path: impl AsRef<Path>) -> Result<Structure, ReadError> {
    let path = path.as_ref();
    let extension = path
        .extension()
        .and_then(|extension| extension.to_str())
        .map(str::to_ascii_lowercase);
    match extension.as_deref() {
        Some("cif" | "mmcif") => read_mmcif(path),
        Some("pdb" | "ent") => pdb::read_pdb(path),
        _ => Err(ReadError::UnknownFormat(path.display().to_string())),
    }
}
