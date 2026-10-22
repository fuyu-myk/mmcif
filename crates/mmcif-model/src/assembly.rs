use std::collections::{HashMap, HashSet};

use crate::expression::{InvalidExpression, OperatorExpression};
use crate::ids::{AssemblyId, LabelAsymId, OperatorId};
use crate::model::{Model, ModelBuilder};
use crate::structure::Structure;
use crate::view::ChainRef;

const MAX_ASSEMBLY_ATOMS: u64 = 1 << 28;

/// A rotation matrix followed by a translation.
///
/// A point `p` maps to `rotation * p + translation`, with `rotation` stored
/// row by row.
///
/// # Examples
///
/// ```
/// use mmcif_model::Transform;
///
/// let shift = Transform {
///     translation: [10.0, 0.0, 0.0],
///     ..Transform::IDENTITY
/// };
/// assert_eq!(shift.apply([1.0, 2.0, 3.0]), [11.0, 2.0, 3.0]);
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Transform {
    /// The rotation matrix, row by row.
    pub rotation: [[f64; 3]; 3],
    /// The translation vector.
    pub translation: [f64; 3],
}

impl Transform {
    /// The transformation that leaves every point where it is.
    pub const IDENTITY: Transform = Transform {
        rotation: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        translation: [0.0; 3],
    };

    /// Returns `point` transformed.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_model::Transform;
    ///
    /// let quarter_turn = Transform {
    ///     rotation: [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]],
    ///     translation: [0.0, 0.0, 1.0],
    /// };
    /// assert_eq!(quarter_turn.apply([1.0, 0.0, 0.0]), [0.0, 1.0, 1.0]);
    /// ```
    pub fn apply(&self, point: [f64; 3]) -> [f64; 3] {
        let r = &self.rotation;
        let t = &self.translation;
        [0, 1, 2].map(|i| r[i][0] * point[0] + r[i][1] * point[1] + r[i][2] * point[2] + t[i])
    }

    /// Returns the transformation that applies `inner` and then `self`.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_model::Transform;
    ///
    /// let shift = Transform {
    ///     translation: [1.0, 0.0, 0.0],
    ///     ..Transform::IDENTITY
    /// };
    /// let quarter_turn = Transform {
    ///     rotation: [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]],
    ///     translation: [0.0; 3],
    /// };
    ///
    /// assert_eq!(quarter_turn.compose(&shift).apply([1.0, 0.0, 0.0]), [0.0, 2.0, 0.0]);
    /// assert_eq!(shift.compose(&quarter_turn).apply([1.0, 0.0, 0.0]), [1.0, 1.0, 0.0]);
    /// ```
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

    /// Returns `true` if every element of the rotation and translation is
    /// within `tolerance` of the identity, inclusive.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_model::Transform;
    ///
    /// let nearly = Transform {
    ///     translation: [0.0, 0.0, 1e-6],
    ///     ..Transform::IDENTITY
    /// };
    /// assert!(nearly.is_identity(1e-5));
    /// assert!(!nearly.is_identity(0.0));
    /// ```
    pub fn is_identity(&self, tolerance: f64) -> bool {
        let close = |a: f64, b: f64| (a - b).abs() <= tolerance;
        (0..3).all(|i| {
            close(self.translation[i], 0.0)
                && (0..3).all(|j| close(self.rotation[i][j], if i == j { 1.0 } else { 0.0 }))
        })
    }
}

/// A named transformation, from `_pdbx_struct_oper_list`.
///
/// # Examples
///
/// ```
/// # let structure = mmcif_model::parse_mmcif(include_str!("../tests/fixtures/pathologies.cif"))?;
/// let operator = &structure.operators()[2];
///
/// assert_eq!(operator.id().as_str(), "3");
/// assert_eq!(operator.kind(), Some("translation"));
/// assert_eq!(operator.transform().apply([1.0, 2.0, 3.0]), [11.0, 2.0, 3.0]);
/// # Ok::<(), mmcif_model::ReadError>(())
/// ```
#[derive(Debug, Clone, PartialEq)]
pub struct Operator {
    pub(crate) id: OperatorId,
    pub(crate) kind: Option<String>,
    pub(crate) transform: Transform,
}

impl Operator {
    /// Creates an operator with no type.
    pub fn new(id: OperatorId, transform: Transform) -> Self {
        Operator {
            id,
            kind: None,
            transform,
        }
    }

    /// Returns the operator id.
    pub fn id(&self) -> &OperatorId {
        &self.id
    }

    /// Returns the operator type from `_pdbx_struct_oper_list.type`, if the
    /// file gives one.
    pub fn kind(&self) -> Option<&str> {
        self.kind.as_deref()
    }

    /// Returns the transformation.
    pub fn transform(&self) -> &Transform {
        &self.transform
    }
}

/// An operator expression and the label chains it applies to, from one row
/// of `_pdbx_struct_assembly_gen`.
///
/// # Examples
///
/// ```
/// use mmcif_model::Generator;
///
/// let generator = Generator::new("(1-60)", vec!["A".into(), "B".into()])?;
///
/// assert_eq!(generator.expression(), "(1-60)");
/// assert_eq!(generator.operators().len(), 60);
/// assert_eq!(generator.asym_ids().len(), 2);
/// # Ok::<(), mmcif_model::InvalidExpression>(())
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Generator {
    pub(crate) expression: String,
    pub(crate) operators: OperatorExpression,
    pub(crate) asym_ids: Vec<LabelAsymId>,
}

impl Generator {
    /// Creates a generator that applies `expression` to the chains
    /// `asym_ids`.
    ///
    /// # Errors
    ///
    /// Returns [`InvalidExpression`] if `expression` cannot be parsed by
    /// [`OperatorExpression::parse`].
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_model::Generator;
    ///
    /// assert!(Generator::new("1,2", vec!["A".into()]).is_ok());
    /// assert!(Generator::new("(1", vec!["A".into()]).is_err());
    /// ```
    pub fn new(expression: &str, asym_ids: Vec<LabelAsymId>) -> Result<Self, InvalidExpression> {
        Ok(Generator {
            expression: expression.to_owned(),
            operators: OperatorExpression::parse(expression)?,
            asym_ids,
        })
    }

    /// Returns the operator expression as given.
    pub fn expression(&self) -> &str {
        &self.expression
    }

    /// Returns the parsed operator expression.
    pub fn operators(&self) -> &OperatorExpression {
        &self.operators
    }

    /// Returns the label ids of the chains the generator applies to.
    pub fn asym_ids(&self) -> &[LabelAsymId] {
        &self.asym_ids
    }
}

/// A biological assembly, from `_pdbx_struct_assembly` and
/// `_pdbx_struct_assembly_gen`.
///
/// [`Structure::assemble`] builds its model.
///
/// # Examples
///
/// ```
/// use mmcif_model::AssemblyId;
///
/// # let structure = mmcif_model::parse_mmcif(include_str!("../tests/fixtures/pathologies.cif"))?;
/// let assembly = structure.assembly(&AssemblyId::new("2")).unwrap();
///
/// assert_eq!(assembly.details(), Some("software_defined_assembly"));
/// assert_eq!(assembly.generators()[0].expression(), "(1,2)");
/// assert_eq!(assembly.operator_count(), 2);
/// # Ok::<(), mmcif_model::ReadError>(())
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Assembly {
    pub(crate) id: AssemblyId,
    pub(crate) details: Option<String>,
    pub(crate) generators: Vec<Generator>,
}

impl Assembly {
    /// Creates an assembly with no details.
    pub fn new(id: AssemblyId, generators: Vec<Generator>) -> Self {
        Assembly {
            id,
            details: None,
            generators,
        }
    }

    /// Returns the assembly id.
    pub fn id(&self) -> &AssemblyId {
        &self.id
    }

    /// Returns the description from `_pdbx_struct_assembly.details`, if the
    /// file gives one.
    pub fn details(&self) -> Option<&str> {
        self.details.as_deref()
    }

    /// Returns the generators, in file order.
    pub fn generators(&self) -> &[Generator] {
        &self.generators
    }

    /// Returns the number of operator combinations across all generators.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_model::{Assembly, AssemblyId, Generator};
    ///
    /// let assembly = Assembly::new(
    ///     AssemblyId::new("1"),
    ///     vec![
    ///         Generator::new("(1-2)(5,6)", vec!["A".into()])?,
    ///         Generator::new("1", vec!["B".into()])?,
    ///     ],
    /// );
    /// assert_eq!(assembly.operator_count(), 5);
    /// # Ok::<(), mmcif_model::InvalidExpression>(())
    /// ```
    pub fn operator_count(&self) -> usize {
        self.generators.iter().map(|g| g.operators.len()).sum()
    }
}

/// An error from [`Structure::assemble`].
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
#[non_exhaustive]
pub enum AssemblyError {
    /// The structure has no assembly with this id.
    #[error("no assembly {0}")]
    UnknownAssembly(AssemblyId),
    /// A generator of the assembly names an operator that the structure does
    /// not define.
    #[error("assembly {assembly} uses undefined operator {operator}")]
    UnknownOperator {
        /// The assembly id.
        assembly: AssemblyId,
        /// The undefined operator id.
        operator: OperatorId,
    },
    /// The assembly would hold more than 2^28 atoms, or memory for its atoms
    /// could not be reserved.
    #[error("assembly {assembly} would hold {atoms} atoms")]
    TooLarge {
        /// The assembly id.
        assembly: AssemblyId,
        /// The number of atoms the assembly would hold.
        atoms: u64,
    },
}

impl Structure {
    /// Builds the model of assembly `id` from the chains of `model`.
    ///
    /// For each generator in order, each operator combination of its
    /// expression is applied, in order, to the chains of `model` that the
    /// generator lists and that have no operators applied, taken in index
    /// order. Listed chains that `model` does not have are skipped. Within a
    /// combination the last operator is applied first: `(X0)(1)` applies `1`
    /// and then `X0`. Each new chain keeps the ids and entity of its
    /// source and records the combination in [`ChainRef::operators`]. The new
    /// model has the number of `model`.
    ///
    /// [`ChainRef::operators`]: crate::ChainRef::operators
    ///
    /// # Errors
    ///
    /// Returns [`AssemblyError::UnknownAssembly`] if the structure has no
    /// assembly `id`, [`AssemblyError::UnknownOperator`] if a generator names
    /// an operator that is not in [`Structure::operators`], and
    /// [`AssemblyError::TooLarge`] if the new model would hold more than 2^28
    /// atoms or memory for them cannot be reserved.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_model::{AssemblyError, AssemblyId};
    ///
    /// # let structure = mmcif_model::parse_mmcif(include_str!("../tests/fixtures/pathologies.cif"))?;
    /// let model = structure.first_model().unwrap();
    /// let dimer = structure.assemble(&AssemblyId::new("2"), model)?;
    ///
    /// let copies: Vec<String> = dimer
    ///     .chains()
    ///     .map(|chain| format!("{}x{}", chain.label_asym(), chain.operators()[0]))
    ///     .collect();
    /// assert_eq!(copies, ["AAx1", "Cx1", "AAx2", "Cx2"]);
    ///
    /// assert_eq!(
    ///     structure.assemble(&AssemblyId::new("9"), model),
    ///     Err(AssemblyError::UnknownAssembly(AssemblyId::new("9")))
    /// );
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn assemble(&self, id: &AssemblyId, model: &Model) -> Result<Model, AssemblyError> {
        let assembly = self
            .assembly(id)
            .ok_or_else(|| AssemblyError::UnknownAssembly(id.clone()))?;
        let operators: HashMap<&str, &Operator> = self
            .operators
            .iter()
            .map(|operator| (operator.id.as_str(), operator))
            .collect();
        let mut plans = Vec::with_capacity(assembly.generators.len());
        let mut total_atoms: u64 = 0;
        for generator in &assembly.generators {
            if let Some(missing) = generator
                .operators
                .operator_ids()
                .find(|operator| !operators.contains_key(operator.as_str()))
            {
                return Err(AssemblyError::UnknownOperator {
                    assembly: id.clone(),
                    operator: missing,
                });
            }
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
        let too_large = || AssemblyError::TooLarge {
            assembly: id.clone(),
            atoms: total_atoms,
        };
        if total_atoms > MAX_ASSEMBLY_ATOMS {
            return Err(too_large());
        }
        let mut builder = ModelBuilder::derived_from(model);
        builder
            .try_reserve_atoms(total_atoms as usize)
            .map_err(|_| too_large())?;
        for (generator, selected) in plans {
            for combination in generator.operators.combinations() {
                let transform = combined_transform(&operators, &combination);
                for chain in &selected {
                    copy_transformed_chain(&mut builder, model, *chain, &transform, &combination);
                }
            }
        }
        Ok(builder.finish())
    }
}

fn combined_transform(
    operators: &HashMap<&str, &Operator>,
    combination: &[OperatorId],
) -> Transform {
    combination
        .iter()
        .fold(Transform::IDENTITY, |acc, operator_id| {
            let operator = operators[operator_id.as_str()];
            acc.compose(&operator.transform)
        })
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
