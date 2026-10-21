use crate::ids::{AssemblyId, LabelAsymId, OperatorId};
use crate::structure::{Chain, Model, Structure};

const MAX_COMBINATIONS: usize = 1 << 20;

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
    pub(crate) combinations: Vec<Vec<OperatorId>>,
    pub(crate) asym_ids: Vec<LabelAsymId>,
}

impl Generator {
    pub fn new(expression: &str, asym_ids: Vec<LabelAsymId>) -> Result<Self, InvalidExpression> {
        Ok(Generator {
            expression: expression.to_owned(),
            combinations: parse_operator_expression(expression)?,
            asym_ids,
        })
    }

    pub fn expression(&self) -> &str {
        &self.expression
    }

    pub fn combinations(&self) -> &[Vec<OperatorId>] {
        &self.combinations
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
        self.generators.iter().map(|g| g.combinations.len()).sum()
    }
}

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
#[error("invalid operator expression {expression:?}: {reason}")]
pub struct InvalidExpression {
    pub expression: String,
    pub reason: &'static str,
}

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum AssemblyError {
    #[error("no assembly {0}")]
    UnknownAssembly(AssemblyId),
    #[error("assembly {assembly} uses undefined operator {operator}")]
    UnknownOperator {
        assembly: AssemblyId,
        operator: OperatorId,
    },
}

pub fn parse_operator_expression(
    expression: &str,
) -> Result<Vec<Vec<OperatorId>>, InvalidExpression> {
    let invalid = |reason| InvalidExpression {
        expression: expression.to_owned(),
        reason,
    };
    let compact: String = expression.chars().filter(|c| !c.is_whitespace()).collect();
    let groups: Vec<&str> = if compact.starts_with('(') {
        let mut groups = Vec::new();
        let mut rest = compact.as_str();
        while !rest.is_empty() {
            let inner = rest
                .strip_prefix('(')
                .ok_or_else(|| invalid("text outside parentheses"))?;
            let close = inner
                .find(')')
                .ok_or_else(|| invalid("unbalanced parentheses"))?;
            groups.push(&inner[..close]);
            rest = &inner[close + 1..];
        }
        groups
    } else {
        vec![compact.as_str()]
    };
    let mut combinations: Vec<Vec<OperatorId>> = vec![Vec::new()];
    for group in groups {
        let members = expand_group(group).map_err(invalid)?;
        let total = combinations.len().saturating_mul(members.len());
        if total > MAX_COMBINATIONS {
            return Err(invalid("too many operator combinations"));
        }
        combinations = combinations
            .iter()
            .flat_map(|prefix| {
                members.iter().map(move |member| {
                    let mut combination = prefix.clone();
                    combination.push(member.clone());
                    combination
                })
            })
            .collect();
    }
    Ok(combinations)
}

fn expand_group(group: &str) -> Result<Vec<OperatorId>, &'static str> {
    let mut members = Vec::new();
    for item in group.split(',') {
        if item.is_empty() || item.contains(['(', ')']) {
            return Err("empty or nested operator item");
        }
        let range = item
            .split_once('-')
            .and_then(|(low, high)| Some((low.parse::<u64>().ok()?, high.parse::<u64>().ok()?)));
        match range {
            Some((low, high)) if low > high => return Err("descending operator range"),
            Some((low, high)) => {
                if (high - low) as usize >= MAX_COMBINATIONS.saturating_sub(members.len()) {
                    return Err("operator range too large");
                }
                members.extend((low..=high).map(|n| OperatorId::new(n.to_string())));
            }
            None => members.push(OperatorId::new(item)),
        }
    }
    Ok(members)
}

impl Structure {
    pub fn assemble(&self, id: &AssemblyId, model: &Model) -> Result<Model, AssemblyError> {
        let assembly = self
            .assembly(id)
            .ok_or_else(|| AssemblyError::UnknownAssembly(id.clone()))?;
        let mut chains = Vec::new();
        for generator in &assembly.generators {
            for combination in &generator.combinations {
                let transform = self.combined_transform(id, combination)?;
                let selected = model.chains.iter().filter(|chain| {
                    chain.operators.is_empty() && generator.asym_ids.contains(&chain.label_asym)
                });
                chains.extend(
                    selected.map(|chain| transformed_chain(chain, &transform, combination)),
                );
            }
        }
        Ok(Model::new(model.number, chains))
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

fn transformed_chain(chain: &Chain, transform: &Transform, combination: &[OperatorId]) -> Chain {
    let mut copy = chain.clone();
    copy.operators = combination.to_vec();
    for atom in copy
        .residues
        .iter_mut()
        .flat_map(|residue| residue.atoms.iter_mut())
    {
        atom.position = transform.apply(atom.position);
    }
    copy
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ids(combinations: &[Vec<OperatorId>]) -> Vec<Vec<&str>> {
        combinations
            .iter()
            .map(|c| c.iter().map(OperatorId::as_str).collect())
            .collect()
    }

    #[test]
    fn parses_lists_ranges_and_products() {
        assert_eq!(ids(&parse_operator_expression("1").unwrap()), [["1"]]);
        assert_eq!(
            ids(&parse_operator_expression("1,2,5").unwrap()),
            [["1"], ["2"], ["5"]]
        );
        assert_eq!(
            ids(&parse_operator_expression("(1-3)").unwrap()),
            [["1"], ["2"], ["3"]]
        );
        assert_eq!(
            ids(&parse_operator_expression("(X0)(1-2,P)").unwrap()),
            [["X0", "1"], ["X0", "2"], ["X0", "P"]]
        );
        assert_eq!(
            parse_operator_expression("(1-60)(61-88)").unwrap().len(),
            60 * 28
        );
    }

    #[test]
    fn rejects_malformed_expressions() {
        for bad in [
            "",
            "1,,2",
            "(1",
            "1)(2",
            "(3-1)",
            "((1))",
            "(1-99999999999)",
        ] {
            assert!(parse_operator_expression(bad).is_err(), "{bad}");
        }
    }

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
