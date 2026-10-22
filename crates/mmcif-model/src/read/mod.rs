mod assemblies;
mod atom_site;
mod entities;
mod numbering;
mod table;

use std::collections::HashMap;
use std::path::Path;
use std::sync::Arc;

use mmcif_cif::{Block, Document};

use crate::entity::Entity;
use crate::error::ReadError;
use crate::ids::{EntityId, LabelAsymId};
use crate::model::Model;
use crate::pdb;
use crate::structure::Structure;

impl Structure {
    /// Reads a structure from the first data block of `document`.
    ///
    /// # Errors
    ///
    /// Returns [`ReadError::NoDataBlock`] if `document` has no data block, and
    /// otherwise the errors of [`Structure::from_block`].
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_cif::Document;
    /// use mmcif_model::Structure;
    ///
    /// # let text = include_str!("../../tests/fixtures/pathologies.cif");
    /// let document = Document::parse(text)?;
    /// let structure = Structure::from_document(&document)?;
    /// assert_eq!(structure.entry_id(), Some("TEST"));
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn from_document(document: &Document) -> Result<Structure, ReadError> {
        Structure::from_block(&document.block(0).ok_or(ReadError::NoDataBlock)?)
    }

    /// Reads a structure from one data block.
    ///
    /// `_atom_site` gives the models, chains, residues and atoms. Atoms are
    /// grouped into models by `pdbx_PDB_model_num`, which defaults to 1, into
    /// chains by `label_asym_id`, and into residues by label number, author
    /// number, insertion code and component, each in order of first
    /// appearance. `label_atom_id` and `label_comp_id` fall back to
    /// `auth_atom_id` and `auth_comp_id`. When the `auth_asym_id` or
    /// `auth_seq_id` column is absent, the label value is used in its place.
    /// A missing occupancy reads as 1 and a missing B-factor as 0.
    ///
    /// When present, `_entry`, `_entity`, `_entity_poly`, `_entity_poly_seq`,
    /// `_struct_asym`, `_pdbx_entity_branch`, `_pdbx_entity_branch_list`,
    /// `_pdbx_entity_branch_link`, `_pdbx_poly_seq_scheme`,
    /// `_pdbx_branch_scheme`, `_pdbx_struct_oper_list`,
    /// `_pdbx_struct_assembly` and `_pdbx_struct_assembly_gen` are read too.
    /// A chain whose atoms name no entity takes the entity that lists it in
    /// `_struct_asym`.
    ///
    /// # Errors
    ///
    /// Returns [`ReadError::MissingCategory`] if `_atom_site` is absent,
    /// [`ReadError::MissingItem`] if a category lacks an item it needs, such
    /// as `_atom_site.Cartn_x`, [`ReadError::InvalidValue`] if a required
    /// value is missing or a value cannot be read,
    /// [`ReadError::InconsistentAuthChain`] if the atoms of one label chain
    /// in a model have different author chain ids, and
    /// [`ReadError::OperatorExpression`] if an operator expression cannot be
    /// parsed.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_cif::Document;
    /// use mmcif_model::Structure;
    ///
    /// # let text = include_str!("../../tests/fixtures/pathologies.cif");
    /// let document = Document::parse(text)?;
    /// let block = document.block_named("TEST").unwrap();
    /// let structure = Structure::from_block(&block)?;
    /// assert_eq!(structure.first_model().unwrap().chain_count(), 4);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
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
            models: models.into_iter().map(Arc::new).collect(),
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

/// Parses a structure from mmCIF text, reading the first data block.
///
/// See [`Structure::from_block`] for how the block is read.
///
/// # Errors
///
/// Returns [`ReadError::Parse`] if `text` is not valid CIF,
/// [`ReadError::NoDataBlock`] if it has no data block, and otherwise the
/// errors of [`Structure::from_block`].
///
/// # Examples
///
/// ```
/// let text = "\
/// data_example
/// loop_
/// _atom_site.label_asym_id
/// _atom_site.label_comp_id
/// _atom_site.label_atom_id
/// _atom_site.type_symbol
/// _atom_site.label_seq_id
/// _atom_site.Cartn_x
/// _atom_site.Cartn_y
/// _atom_site.Cartn_z
/// A GLY N N 1 0.0 0.0 0.0
/// A GLY CA C 1 1.5 0.0 0.0
/// ";
/// let structure = mmcif_model::parse_mmcif(text)?;
/// let model = structure.first_model().unwrap();
///
/// assert_eq!(model.atom_count(), 2);
/// assert_eq!(model.chain("A").unwrap().residue_count(), 1);
/// # Ok::<(), mmcif_model::ReadError>(())
/// ```
pub fn parse_mmcif(text: &str) -> Result<Structure, ReadError> {
    Structure::from_document(&Document::parse(text)?)
}

/// Reads a structure from an mmCIF file, reading the first data block.
///
/// # Errors
///
/// Returns [`ReadError::Cif`] if the file cannot be read or is not valid CIF,
/// [`ReadError::NoDataBlock`] if it has no data block, and otherwise the
/// errors of [`Structure::from_block`].
///
/// # Examples
///
/// ```
/// let structure = mmcif_model::read_mmcif("tests/fixtures/pathologies.cif")?;
/// assert_eq!(structure.entry_id(), Some("TEST"));
/// # Ok::<(), mmcif_model::ReadError>(())
/// ```
pub fn read_mmcif(path: impl AsRef<Path>) -> Result<Structure, ReadError> {
    Structure::from_document(&Document::open(path)?)
}

/// Parses a structure from BinaryCIF data, reading the first data block.
///
/// # Errors
///
/// Returns [`ReadError::Cif`] if `bytes` are not valid BinaryCIF,
/// [`ReadError::NoDataBlock`] if they hold no data block, and otherwise the
/// errors of [`Structure::from_block`].
///
/// # Examples
///
/// ```no_run
/// let bytes = std::fs::read("1crn.bcif")?;
/// let structure = mmcif_model::parse_bcif(&bytes)?;
/// println!("{} atoms", structure.first_model().unwrap().atom_count());
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub fn parse_bcif(bytes: &[u8]) -> Result<Structure, ReadError> {
    Structure::from_document(&Document::from_bcif(bytes).map_err(mmcif_cif::Error::from)?)
}

/// Reads a structure from a BinaryCIF file, reading the first data block.
///
/// # Errors
///
/// Returns [`ReadError::Cif`] if the file cannot be read or is not valid
/// BinaryCIF, [`ReadError::NoDataBlock`] if it has no data block, and
/// otherwise the errors of [`Structure::from_block`].
///
/// # Examples
///
/// ```no_run
/// let structure = mmcif_model::read_bcif("1crn.bcif")?;
/// println!("{} atoms", structure.first_model().unwrap().atom_count());
/// # Ok::<(), mmcif_model::ReadError>(())
/// ```
pub fn read_bcif(path: impl AsRef<Path>) -> Result<Structure, ReadError> {
    Structure::from_document(&Document::open_bcif(path)?)
}

/// Reads a structure from a file, choosing the format by its extension.
///
/// Extensions are matched ignoring case: `cif` and `mmcif` are read with
/// [`read_mmcif`], `bcif` with [`read_bcif`], and `pdb` and `ent` with
/// [`read_pdb`](crate::read_pdb).
///
/// # Errors
///
/// Returns [`ReadError::UnknownFormat`] if the extension is missing or not
/// one of these, and otherwise the errors of the chosen reader.
///
/// # Examples
///
/// ```
/// use mmcif_model::ReadError;
///
/// let structure = mmcif_model::read("tests/fixtures/legacy.pdb")?;
/// assert_eq!(structure.entry_id(), Some("1ABC"));
///
/// assert!(matches!(
///     mmcif_model::read("structure.xyz"),
///     Err(ReadError::UnknownFormat(_))
/// ));
/// # Ok::<(), ReadError>(())
/// ```
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
