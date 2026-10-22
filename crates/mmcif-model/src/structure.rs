use std::sync::Arc;

use crate::assembly::{Assembly, Operator};
use crate::entity::Entity;
use crate::ids::{AssemblyId, ModelNumber};
use crate::model::Model;
use crate::numbering::Numbering;

/// A macromolecular structure read from one mmCIF, BinaryCIF or PDB entry.
///
/// It holds the models of the entry together with the entities, operators,
/// assemblies and residue numbering they share.
///
/// # Examples
///
/// ```
/// # let text = include_str!("../tests/fixtures/pathologies.cif");
/// let structure = mmcif_model::parse_mmcif(text)?;
///
/// assert_eq!(structure.entry_id(), Some("TEST"));
/// assert_eq!(structure.models().len(), 2);
/// assert_eq!(structure.entities().len(), 4);
/// assert_eq!(structure.first_model().unwrap().atom_count(), 16);
/// # Ok::<(), mmcif_model::ReadError>(())
/// ```
#[derive(Debug, Clone, PartialEq)]
pub struct Structure {
    pub(crate) entry_id: Option<String>,
    pub(crate) entities: Vec<Entity>,
    pub(crate) models: Vec<Arc<Model>>,
    pub(crate) operators: Vec<Operator>,
    pub(crate) assemblies: Vec<Assembly>,
    pub(crate) numbering: Numbering,
}

impl Structure {
    /// Returns the entry id from `_entry.id` or the PDB `HEADER` record, if
    /// the file gives one.
    pub fn entry_id(&self) -> Option<&str> {
        self.entry_id.as_deref()
    }

    /// Returns the entities, in file order.
    pub fn entities(&self) -> &[Entity] {
        &self.entities
    }

    /// Returns the entity with id `id`, if there is one.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_model::EntityKind;
    ///
    /// # let structure = mmcif_model::parse_mmcif(include_str!("../tests/fixtures/pathologies.cif"))?;
    /// let zinc = structure.entity("3").unwrap();
    /// assert_eq!(zinc.kind(), &EntityKind::NonPolymer);
    /// assert!(structure.entity("9").is_none());
    /// # Ok::<(), mmcif_model::ReadError>(())
    /// ```
    pub fn entity(&self, id: &str) -> Option<&Entity> {
        self.entities
            .iter()
            .find(|entity| entity.id().as_str() == id)
    }

    /// Returns the models, in file order.
    pub fn models(&self) -> &[Arc<Model>] {
        &self.models
    }

    /// Returns the model numbered `number`, if there is one.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_model::ModelNumber;
    ///
    /// # let structure = mmcif_model::parse_mmcif(include_str!("../tests/fixtures/pathologies.cif"))?;
    /// let second = structure.model(ModelNumber::new(2)).unwrap();
    /// assert_eq!(second.atom_count(), 1);
    /// assert!(structure.model(ModelNumber::new(3)).is_none());
    /// # Ok::<(), mmcif_model::ReadError>(())
    /// ```
    pub fn model(&self, number: ModelNumber) -> Option<&Arc<Model>> {
        self.models.iter().find(|model| model.number == number)
    }

    /// Returns the first model, or `None` if the structure has no models.
    pub fn first_model(&self) -> Option<&Arc<Model>> {
        self.models.first()
    }

    /// Returns the transformation operators from `_pdbx_struct_oper_list`, in
    /// file order.
    pub fn operators(&self) -> &[Operator] {
        &self.operators
    }

    /// Returns the operator with id `id`, if there is one.
    ///
    /// # Examples
    ///
    /// ```
    /// # let structure = mmcif_model::parse_mmcif(include_str!("../tests/fixtures/pathologies.cif"))?;
    /// let shift = structure.operator("3").unwrap();
    /// assert_eq!(shift.transform().translation, [10.0, 0.0, 0.0]);
    /// # Ok::<(), mmcif_model::ReadError>(())
    /// ```
    pub fn operator(&self, id: &str) -> Option<&Operator> {
        self.operators
            .iter()
            .find(|operator| operator.id().as_str() == id)
    }

    /// Returns the assemblies, in file order.
    pub fn assemblies(&self) -> &[Assembly] {
        &self.assemblies
    }

    /// Returns the assembly with id `id`, if there is one.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_model::AssemblyId;
    ///
    /// # let structure = mmcif_model::parse_mmcif(include_str!("../tests/fixtures/pathologies.cif"))?;
    /// let assembly = structure.assembly(&AssemblyId::new("1")).unwrap();
    /// assert_eq!(assembly.details(), Some("author_defined_assembly"));
    /// # Ok::<(), mmcif_model::ReadError>(())
    /// ```
    pub fn assembly(&self, id: &AssemblyId) -> Option<&Assembly> {
        self.assemblies.iter().find(|assembly| assembly.id() == id)
    }

    /// Returns the mapping between label and author residue numbering.
    pub fn numbering(&self) -> &Numbering {
        &self.numbering
    }
}
