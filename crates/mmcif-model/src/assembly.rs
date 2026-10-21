use std::collections::HashSet;

use crate::expression::{InvalidExpression, OperatorExpression};
use crate::ids::{AssemblyId, LabelAsymId, OperatorId};
use crate::model::{Model, ModelBuilder};
use crate::structure::Structure;
use crate::view::ChainRef;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Transform {
    pub rotation: [[f64; 3]; 3],
    pub translation: [f64; 3],
}

impl Transform {
    pub const IDENTITY: Transform = Transform {
        rotation: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        translation: [0.0; 3],
    };

    pub fn apply(&self, point: [f64; 3]) -> [f64; 3] {
        let r = &self.rotation;
        let t = &self.translation;
        [0, 1, 2].map(|i| r[i][0] * point[0] + r[i][1] * point[1] + r[i][2] * point[2] + t[i])
    }

    pub fn compose(&self, inner: &Transform) -> Transform {
        let a = &self.rotation;
        let b = &inner.rotation;
        let rotation =
            [0, 1, 2].map(|i| [0, 1, 2].map(|j| (0..3).map(|k| a[i][k] * b[k][j]).sum()));
        let moved = self.apply(inner.translation);
        Transform {
            rotation,
            translation: moved,
        }
    }

    pub fn is_identity(&self, tolerance: f64) -> bool {
        let close = |a: f64, b: f64| (a - b).abs() <= tolerance;
        (0..3).all(|i| {
            close(self.translation[i], 0.0)
                && (0..3).all(|j| close(self.rotation[i][j], if i == j { 1.0 } else { 0.0 }))
        })
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Operator {
    pub(crate) id: OperatorId,
    pub(crate) kind: Option<String>,
    pub(crate) transform: Transform,
}

impl Operator {
    pub fn new(id: OperatorId, transform: Transform) -> Self {
        Operator {
            id,
            kind: None,
            transform,
        }
    }

    pub fn id(&self) -> &OperatorId {
        &self.id
    }

    pub fn kind(&self) -> Option<&str> {
        self.kind.as_deref()
    }

    pub fn transform(&self) -> &Transform {
        &self.transform
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Generator {
    pub(crate) expression: String,
    pub(crate) operators: OperatorExpression,
    pub(crate) asym_ids: Vec<LabelAsymId>,
}

impl Generator {
    pub fn new(expression: &str, asym_ids: Vec<LabelAsymId>) -> Result<Self, InvalidExpression> {
        Ok(Generator {
            expression: expression.to_owned(),
            operators: OperatorExpression::parse(expression)?,
            asym_ids,
        })
    }

    pub fn expression(&self) -> &str {
        &self.expression
    }

    pub fn operators(&self) -> &OperatorExpression {
        &self.operators
    }

    pub fn asym_ids(&self) -> &[LabelAsymId] {
        &self.asym_ids
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Assembly {
    pub(crate) id: AssemblyId,
    pub(crate) details: Option<String>,
    pub(crate) generators: Vec<Generator>,
}

impl Assembly {
    pub fn new(id: AssemblyId, generators: Vec<Generator>) -> Self {
        Assembly {
            id,
            details: None,
            generators,
        }
    }

    pub fn id(&self) -> &AssemblyId {
        &self.id
    }

    pub fn details(&self) -> Option<&str> {
        self.details.as_deref()
    }

    pub fn generators(&self) -> &[Generator] {
        &self.generators
    }

    pub fn operator_count(&self) -> usize {
        self.generators.iter().map(|g| g.operators.len()).sum()
    }
}

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
#[non_exhaustive]
pub enum AssemblyError {
    #[error("no assembly {0}")]
    UnknownAssembly(AssemblyId),
    #[error("assembly {assembly} uses undefined operator {operator}")]
    UnknownOperator {
        assembly: AssemblyId,
        operator: OperatorId,
    },
    #[error("assembly {assembly} would hold {atoms} atoms")]
    TooLarge { assembly: AssemblyId, atoms: u64 },
}

impl Structure {
    pub fn assemble(&self, id: &AssemblyId, model: &Model) -> Result<Model, AssemblyError> {
        let assembly = self
            .assembly(id)
            .ok_or_else(|| AssemblyError::UnknownAssembly(id.clone()))?;
        let mut plans = Vec::with_capacity(assembly.generators.len());
        let mut total_atoms: u64 = 0;
        for generator in &assembly.generators {
            let wanted: HashSet<&LabelAsymId> = generator.asym_ids.iter().collect();
            let selected: Vec<ChainRef<'_>> = model
                .chains()
                .filter(|chain| chain.operators().is_empty() && wanted.contains(chain.label_asym()))
                .collect();
            let atoms: u64 = selected.iter().map(|chain| chain.atom_count() as u64).sum();
            total_atoms =
                total_atoms.saturating_add(atoms.saturating_mul(generator.operators.len() as u64));
            plans.push((generator, selected));
        }
        if total_atoms > u64::from(u32::MAX) {
            return Err(AssemblyError::TooLarge {
                assembly: id.clone(),
                atoms: total_atoms,
            });
        }
        let mut builder = ModelBuilder::derived_from(model);
        builder.reserve_atoms(total_atoms as usize);
        for (generator, selected) in plans {
            for combination in generator.operators.combinations() {
                let transform = self.combined_transform(id, &combination)?;
                for chain in &selected {
                    copy_transformed_chain(&mut builder, model, *chain, &transform, &combination);
                }
            }
        }
        Ok(builder.finish())
    }

    fn combined_transform(
        &self,
        assembly: &AssemblyId,
        combination: &[OperatorId],
    ) -> Result<Transform, AssemblyError> {
        combination
            .iter()
            .try_fold(Transform::IDENTITY, |acc, operator_id| {
                let operator = self.operator(operator_id.as_str()).ok_or_else(|| {
                    AssemblyError::UnknownOperator {
                        assembly: assembly.clone(),
                        operator: operator_id.clone(),
                    }
                })?;
                Ok(acc.compose(&operator.transform))
            })
    }
}

fn copy_transformed_chain(
    builder: &mut ModelBuilder,
    model: &Model,
    chain: ChainRef<'_>,
    transform: &Transform,
    combination: &[OperatorId],
) {
    let slot = builder.new_chain(
        chain.label_asym().clone(),
        chain.auth_asym().clone(),
        chain.entity().cloned(),
        combination.to_vec(),
    );
    for residue in chain.residues() {
        let residue_slot = builder.new_residue(
            slot,
            residue.comp().clone(),
            residue.label_seq(),
            residue.auth_seq(),
            residue.ins_code(),
            residue.is_hetero(),
        );
        for atom in residue.atom_range() {
            builder.push_copied_atom(
                residue_slot,
                model,
                atom,
                transform.apply(model.positions()[atom]),
            );
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn composition_applies_the_inner_transform_first() {
        let shift = Transform {
            translation: [1.0, 0.0, 0.0],
            ..Transform::IDENTITY
        };
        let quarter_turn = Transform {
            rotation: [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]],
            translation: [0.0; 3],
        };
        let point = [1.0, 0.0, 0.0];
        assert_eq!(quarter_turn.compose(&shift).apply(point), [0.0, 2.0, 0.0]);
        assert_eq!(shift.compose(&quarter_turn).apply(point), [1.0, 1.0, 0.0]);
        assert!(
            Transform::IDENTITY
                .compose(&Transform::IDENTITY)
                .is_identity(0.0)
        );
    }
}
