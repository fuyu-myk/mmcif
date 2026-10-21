use std::sync::Arc;

use crate::assembly::{Assembly, Operator};
use crate::entity::Entity;
use crate::ids::{AssemblyId, ModelNumber};
use crate::model::Model;
use crate::numbering::Numbering;

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

    pub fn models(&self) -> &[Arc<Model>] {
        &self.models
    }

    pub fn model(&self, number: ModelNumber) -> Option<&Arc<Model>> {
        self.models.iter().find(|model| model.number == number)
    }

    pub fn first_model(&self) -> Option<&Arc<Model>> {
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
