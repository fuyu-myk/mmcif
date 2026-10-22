use crate::ids::OperatorId;

const MAX_COMBINATIONS: u64 = 1 << 20;
const MAX_GROUPS: usize = 8;

/// An error from parsing an operator expression with
/// [`OperatorExpression::parse`].
///
/// # Examples
///
/// ```
/// use mmcif_model::OperatorExpression;
///
/// let error = OperatorExpression::parse("(3-1)").unwrap_err();
/// assert_eq!(error.expression, "(3-1)");
/// assert_eq!(error.reason, "descending operator range");
/// ```
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
#[error("invalid operator expression {expression:?}: {reason}")]
pub struct InvalidExpression {
    /// The expression as given.
    pub expression: String,
    /// A short description of the problem.
    pub reason: &'static str,
}

#[derive(Debug, Clone, PartialEq, Eq)]
enum Member {
    Single(OperatorId),
    Range { low: u64, high: u64 },
}

impl Member {
    fn len(&self) -> u64 {
        match self {
            Member::Single(_) => 1,
            Member::Range { low, high } => (high - low).saturating_add(1),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct Group {
    members: Vec<Member>,
    len: u64,
}

impl Group {
    fn parse(text: &str) -> Result<Group, &'static str> {
        let mut members = Vec::new();
        let mut len: u64 = 0;
        for item in text.split(',') {
            if item.is_empty() || item.contains(['(', ')']) {
                return Err("empty or nested operator item");
            }
            let range = item.split_once('-').and_then(|(low, high)| {
                Some((low.parse::<u64>().ok()?, high.parse::<u64>().ok()?))
            });
            let member = match range {
                Some((low, high)) if low > high => return Err("descending operator range"),
                Some((low, high)) => Member::Range { low, high },
                None => Member::Single(OperatorId::new(item)),
            };
            len = len
                .checked_add(member.len())
                .filter(|&len| len <= MAX_COMBINATIONS)
                .ok_or("too many operator combinations")?;
            members.push(member);
        }
        Ok(Group { members, len })
    }

    fn get(&self, mut index: u64) -> OperatorId {
        for member in &self.members {
            let len = member.len();
            if index < len {
                return match member {
                    Member::Single(id) => id.clone(),
                    Member::Range { low, .. } => OperatorId::new((low + index).to_string()),
                };
            }
            index -= len;
        }
        unreachable!("operator index is below the group length")
    }
}

/// A parsed `_pdbx_struct_assembly_gen.oper_expression`.
///
/// An expression is either one group, such as `1,2,5` or `1-60`, or a
/// sequence of parenthesised groups, such as `(X0)(1-60)`. Each combination
/// takes one operator from every group.
///
/// # Examples
///
/// ```
/// use mmcif_model::OperatorExpression;
///
/// let expression = OperatorExpression::parse("(1-2)(5,6)")?;
/// assert_eq!(expression.len(), 4);
///
/// let combinations: Vec<String> = expression
///     .combinations()
///     .map(|ids| ids.iter().map(|id| id.as_str()).collect::<Vec<_>>().join("x"))
///     .collect();
/// assert_eq!(combinations, ["1x5", "1x6", "2x5", "2x6"]);
/// # Ok::<(), mmcif_model::InvalidExpression>(())
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct OperatorExpression {
    groups: Vec<Group>,
    len: u64,
}

impl OperatorExpression {
    /// Parses an operator expression.
    ///
    /// Whitespace is ignored. Items within a group are separated by commas.
    /// An item of two unsigned integers joined by `-` is an inclusive range
    /// of operator ids, and any other item is a single operator id.
    ///
    /// # Errors
    ///
    /// Returns [`InvalidExpression`] if the parentheses are unbalanced or
    /// nested, text lies outside the parentheses of a parenthesised
    /// expression, an item is empty, a range is descending, there are more
    /// than 8 groups, or the expression has more than 2^20 combinations.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_model::OperatorExpression;
    ///
    /// assert_eq!(OperatorExpression::parse("1, 2, 5")?.len(), 3);
    /// assert_eq!(OperatorExpression::parse("(1-60)(61-88)")?.len(), 60 * 28);
    ///
    /// assert!(OperatorExpression::parse("").is_err());
    /// assert!(OperatorExpression::parse("1,,2").is_err());
    /// assert!(OperatorExpression::parse("(1)2").is_err());
    /// # Ok::<(), mmcif_model::InvalidExpression>(())
    /// ```
    pub fn parse(expression: &str) -> Result<OperatorExpression, InvalidExpression> {
        let invalid = |reason| InvalidExpression {
            expression: expression.to_owned(),
            reason,
        };
        let compact: String = expression.chars().filter(|c| !c.is_whitespace()).collect();
        let mut texts: Vec<&str> = Vec::new();
        if compact.starts_with('(') {
            let mut rest = compact.as_str();
            while !rest.is_empty() {
                let inner = rest
                    .strip_prefix('(')
                    .ok_or_else(|| invalid("text outside parentheses"))?;
                let close = inner
                    .find(')')
                    .ok_or_else(|| invalid("unbalanced parentheses"))?;
                if texts.len() == MAX_GROUPS {
                    return Err(invalid("too many operator groups"));
                }
                texts.push(&inner[..close]);
                rest = &inner[close + 1..];
            }
        } else {
            texts.push(&compact);
        }
        let mut len: u64 = 1;
        let mut groups = Vec::with_capacity(texts.len());
        for text in texts {
            let group = Group::parse(text).map_err(invalid)?;
            len = len
                .checked_mul(group.len)
                .filter(|&len| len <= MAX_COMBINATIONS)
                .ok_or_else(|| invalid("too many operator combinations"))?;
            groups.push(group);
        }
        Ok(OperatorExpression { groups, len })
    }

    /// Returns the number of operator combinations.
    pub fn len(&self) -> usize {
        self.len as usize
    }

    /// Returns `true` if the expression has no operator combinations.
    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    /// Returns every operator id the expression names, group by group, with
    /// ranges expanded.
    ///
    /// An id named in more than one group is returned once per group.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_model::OperatorExpression;
    ///
    /// let expression = OperatorExpression::parse("(X0)(1-2,P)")?;
    /// let ids: Vec<String> = expression.operator_ids().map(|id| id.to_string()).collect();
    /// assert_eq!(ids, ["X0", "1", "2", "P"]);
    /// # Ok::<(), mmcif_model::InvalidExpression>(())
    /// ```
    pub fn operator_ids(&self) -> impl Iterator<Item = OperatorId> + '_ {
        self.groups
            .iter()
            .flat_map(|group| (0..group.len).map(move |index| group.get(index)))
    }

    /// Returns every combination of one operator from each group, as ids in
    /// group order.
    ///
    /// The combinations vary the last group fastest.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_model::OperatorExpression;
    ///
    /// let expression = OperatorExpression::parse("(X0)(1-2)")?;
    /// let combinations: Vec<Vec<String>> = expression
    ///     .combinations()
    ///     .map(|ids| ids.iter().map(|id| id.to_string()).collect())
    ///     .collect();
    /// assert_eq!(combinations, [["X0", "1"], ["X0", "2"]]);
    /// # Ok::<(), mmcif_model::InvalidExpression>(())
    /// ```
    pub fn combinations(&self) -> impl ExactSizeIterator<Item = Vec<OperatorId>> + '_ {
        (0..self.len as usize).map(move |index| {
            let mut rest = index as u64;
            let mut combination: Vec<OperatorId> = self
                .groups
                .iter()
                .rev()
                .map(|group| {
                    let operator = group.get(rest % group.len);
                    rest /= group.len;
                    operator
                })
                .collect();
            combination.reverse();
            combination
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ids(expression: &str) -> Vec<Vec<String>> {
        OperatorExpression::parse(expression)
            .unwrap()
            .combinations()
            .map(|c| c.iter().map(|id| id.to_string()).collect())
            .collect()
    }

    #[test]
    fn parses_lists_ranges_and_products() {
        assert_eq!(ids("1"), [["1"]]);
        assert_eq!(ids("1,2,5"), [["1"], ["2"], ["5"]]);
        assert_eq!(ids("(1-3)"), [["1"], ["2"], ["3"]]);
        assert_eq!(ids("(X0)(1-2,P)"), [["X0", "1"], ["X0", "2"], ["X0", "P"]]);
        assert_eq!(
            ids("(1-2)(5,6)"),
            [["1", "5"], ["1", "6"], ["2", "5"], ["2", "6"]]
        );
        assert_eq!(
            OperatorExpression::parse("(1-60)(61-88)").unwrap().len(),
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
            assert!(OperatorExpression::parse(bad).is_err(), "{bad}");
        }
    }

    #[test]
    fn bounds_expressions_without_expanding_them() {
        assert!(OperatorExpression::parse("(1-1024)(1-1024)(1)(1)").is_ok());
        assert!(OperatorExpression::parse("(1-1024)(1-1025)").is_err());
        assert!(OperatorExpression::parse(&"(1)".repeat(MAX_GROUPS + 1)).is_err());
        assert!(OperatorExpression::parse("1-18446744073709551615").is_err());
        assert!(OperatorExpression::parse("0-18446744073709551615").is_err());
    }
}
