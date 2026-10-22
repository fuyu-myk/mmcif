use std::cmp::Ordering;
use std::collections::BTreeSet;

use mmcif_model::CompId;

use crate::error::IrError;
use crate::ids::ResidueRef;

/// A non-empty set of chemical component ids.
///
/// A repeated component is stored once, and the set iterates in ascending
/// order.
///
/// # Examples
///
/// ```
/// use mmcif_ir::AlphabetSet;
/// use mmcif_model::CompId;
///
/// let acidic = AlphabetSet::new([CompId::new("GLU"), CompId::new("ASP"), CompId::new("GLU")])?;
///
/// assert_eq!(acidic.len(), 2);
/// assert!(acidic.contains(&CompId::new("ASP")));
/// # Ok::<(), mmcif_ir::IrError>(())
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
#[allow(clippy::len_without_is_empty)]
pub struct AlphabetSet(BTreeSet<CompId>);

impl AlphabetSet {
    /// Creates a set from `components`.
    ///
    /// # Errors
    ///
    /// Returns [`IrError::EmptyAlphabet`] if `components` is empty.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_ir::{AlphabetSet, IrError};
    /// use mmcif_model::CompId;
    ///
    /// let acidic = AlphabetSet::new([CompId::new("ASP"), CompId::new("GLU")])?;
    /// assert_eq!(acidic.len(), 2);
    ///
    /// assert_eq!(AlphabetSet::new([]), Err(IrError::EmptyAlphabet));
    /// # Ok::<(), IrError>(())
    /// ```
    pub fn new(components: impl IntoIterator<Item = CompId>) -> Result<AlphabetSet, IrError> {
        let components: BTreeSet<CompId> = components.into_iter().collect();
        if components.is_empty() {
            return Err(IrError::EmptyAlphabet);
        }
        Ok(AlphabetSet(components))
    }

    /// Returns the number of components in the set.
    pub fn len(&self) -> usize {
        self.0.len()
    }

    /// Returns `true` if the set holds `component`.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_ir::AlphabetSet;
    /// use mmcif_model::CompId;
    ///
    /// let acidic = AlphabetSet::new([CompId::new("ASP"), CompId::new("GLU")])?;
    ///
    /// assert!(acidic.contains(&CompId::new("GLU")));
    /// assert!(!acidic.contains(&CompId::new("LYS")));
    /// # Ok::<(), mmcif_ir::IrError>(())
    /// ```
    pub fn contains(&self, component: &CompId) -> bool {
        self.0.contains(component)
    }

    /// Returns an iterator over the components, in ascending order.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_ir::AlphabetSet;
    /// use mmcif_model::CompId;
    ///
    /// let acidic = AlphabetSet::new([CompId::new("GLU"), CompId::new("ASP")])?;
    /// let names: Vec<&str> = acidic.iter().map(CompId::as_str).collect();
    ///
    /// assert_eq!(names, ["ASP", "GLU"]);
    /// # Ok::<(), mmcif_ir::IrError>(())
    /// ```
    pub fn iter(&self) -> impl ExactSizeIterator<Item = &CompId> {
        self.0.iter()
    }
}

/// The chemical components a slot may hold.
///
/// # Examples
///
/// ```
/// use mmcif_ir::{Identity, Level};
/// use mmcif_model::CompId;
///
/// let acidic = Identity::one_of([CompId::new("ASP"), CompId::new("GLU")])?;
///
/// assert!(acidic.allows(&CompId::new("ASP")));
/// assert!(!acidic.allows(&CompId::new("LYS")));
/// assert_eq!(acidic.level(), Level::Restricted);
/// # Ok::<(), mmcif_ir::IrError>(())
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum Identity {
    /// Any component.
    Any,
    /// Any component in the set.
    OneOf(AlphabetSet),
    /// Exactly this component.
    Exactly(CompId),
}

impl Identity {
    /// Creates an identity that allows each of `components`.
    ///
    /// Returns [`Identity::Exactly`] when `components` name a single distinct
    /// component, and [`Identity::OneOf`] otherwise.
    ///
    /// # Errors
    ///
    /// Returns [`IrError::EmptyAlphabet`] if `components` is empty.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_ir::{Identity, IrError};
    /// use mmcif_model::CompId;
    ///
    /// let alanine = Identity::one_of([CompId::new("ALA"), CompId::new("ALA")])?;
    /// assert_eq!(alanine, Identity::Exactly(CompId::new("ALA")));
    ///
    /// let acidic = Identity::one_of([CompId::new("ASP"), CompId::new("GLU")])?;
    /// assert!(matches!(acidic, Identity::OneOf(_)));
    ///
    /// assert_eq!(Identity::one_of([]), Err(IrError::EmptyAlphabet));
    /// # Ok::<(), IrError>(())
    /// ```
    pub fn one_of(components: impl IntoIterator<Item = CompId>) -> Result<Identity, IrError> {
        let allowed = AlphabetSet::new(components)?;
        match allowed.len() {
            1 => Ok(Identity::Exactly(
                allowed
                    .iter()
                    .next()
                    .expect("the set holds one component")
                    .clone(),
            )),
            _ => Ok(Identity::OneOf(allowed)),
        }
    }

    /// Returns `true` if the identity allows `component`.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_ir::Identity;
    /// use mmcif_model::CompId;
    ///
    /// let cysteine = Identity::Exactly(CompId::new("CYS"));
    ///
    /// assert!(cysteine.allows(&CompId::new("CYS")));
    /// assert!(!cysteine.allows(&CompId::new("SER")));
    /// assert!(Identity::Any.allows(&CompId::new("SER")));
    /// ```
    pub fn allows(&self, component: &CompId) -> bool {
        match self {
            Identity::Any => true,
            Identity::OneOf(allowed) => allowed.contains(component),
            Identity::Exactly(fixed) => fixed == component,
        }
    }

    /// Returns the determinacy level of the identity.
    ///
    /// [`Identity::Any`] is [`Level::Free`], a set of several components is
    /// [`Level::Restricted`], and a set of one component or
    /// [`Identity::Exactly`] is [`Level::Fixed`].
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_ir::{AlphabetSet, Identity, Level};
    /// use mmcif_model::CompId;
    ///
    /// let acidic = Identity::one_of([CompId::new("ASP"), CompId::new("GLU")])?;
    /// let alanine = Identity::OneOf(AlphabetSet::new([CompId::new("ALA")])?);
    ///
    /// assert_eq!(Identity::Any.level(), Level::Free);
    /// assert_eq!(acidic.level(), Level::Restricted);
    /// assert_eq!(alanine.level(), Level::Fixed);
    /// # Ok::<(), mmcif_ir::IrError>(())
    /// ```
    pub fn level(&self) -> Level {
        match self {
            Identity::Any => Level::Free,
            Identity::OneOf(allowed) if allowed.len() == 1 => Level::Fixed,
            Identity::OneOf(_) => Level::Restricted,
            Identity::Exactly(_) => Level::Fixed,
        }
    }
}

/// Where the position of a slot comes from.
///
/// # Examples
///
/// ```
/// use mmcif_ir::{Geometry, ResidueRef};
/// use mmcif_model::{LabelResidueKey, LabelSeqId};
///
/// let residue = ResidueRef::Label(LabelResidueKey {
///     asym: "A".into(),
///     seq: LabelSeqId::new(57).unwrap(),
/// });
/// let anchored = Geometry::Anchored(residue.clone());
///
/// assert_eq!(anchored.anchor(), Some(&residue));
/// assert_eq!(Geometry::Free.anchor(), None);
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum Geometry {
    /// The position is not taken from a structure.
    Free,
    /// The position of a residue in an existing structure.
    Anchored(ResidueRef),
}

impl Geometry {
    /// Returns the residue the geometry is anchored to, if any.
    pub fn anchor(&self) -> Option<&ResidueRef> {
        match self {
            Geometry::Free => None,
            Geometry::Anchored(reference) => Some(reference),
        }
    }
}

/// A position in a design, described by its identity and its geometry.
///
/// # Examples
///
/// ```
/// use mmcif_ir::{Geometry, Identity, Slot};
/// use mmcif_model::CompId;
///
/// let cysteine = Slot::new(Identity::Exactly(CompId::new("CYS")), Geometry::Free);
///
/// assert!(cysteine.identity.allows(&CompId::new("CYS")));
/// assert_eq!(cysteine.geometry, Geometry::Free);
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub struct Slot {
    /// The components the slot may hold.
    pub identity: Identity,
    /// Where the position of the slot comes from.
    pub geometry: Geometry,
}

impl Slot {
    /// Creates a slot from an identity and a geometry.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_ir::{Geometry, Identity, Slot};
    ///
    /// let slot = Slot::new(Identity::Any, Geometry::Free);
    ///
    /// assert_eq!(slot.identity, Identity::Any);
    /// ```
    pub fn new(identity: Identity, geometry: Geometry) -> Slot {
        Slot { identity, geometry }
    }

    /// Creates a slot that may hold any component at a free position.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_ir::{Geometry, Identity, Slot};
    ///
    /// assert_eq!(Slot::de_novo(), Slot::new(Identity::Any, Geometry::Free));
    /// ```
    pub fn de_novo() -> Slot {
        Slot::new(Identity::Any, Geometry::Free)
    }

    /// Creates a slot that holds exactly `component` at the position of
    /// `anchor`.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_ir::{Geometry, Identity, ResidueRef, Slot};
    /// use mmcif_model::{CompId, LabelResidueKey, LabelSeqId};
    ///
    /// let residue = ResidueRef::Label(LabelResidueKey {
    ///     asym: "A".into(),
    ///     seq: LabelSeqId::new(57).unwrap(),
    /// });
    /// let histidine = Slot::context(CompId::new("HIS"), residue.clone());
    ///
    /// assert_eq!(histidine.identity, Identity::Exactly(CompId::new("HIS")));
    /// assert_eq!(histidine.geometry, Geometry::Anchored(residue));
    /// ```
    pub fn context(component: CompId, anchor: ResidueRef) -> Slot {
        Slot::new(Identity::Exactly(component), Geometry::Anchored(anchor))
    }

    pub(crate) fn determinacy(&self, constrained: bool) -> Determinacy {
        let geometry = match (&self.geometry, constrained) {
            (Geometry::Anchored(_), _) => Level::Fixed,
            (Geometry::Free, true) => Level::Restricted,
            (Geometry::Free, false) => Level::Free,
        };
        Determinacy {
            identity: self.identity.level(),
            geometry,
        }
    }
}

/// How far one axis of a slot is determined.
///
/// Levels are ordered from [`Free`](Level::Free) to [`Fixed`](Level::Fixed).
///
/// # Examples
///
/// ```
/// use mmcif_ir::Level;
///
/// assert!(Level::Free < Level::Restricted);
/// assert!(Level::Restricted < Level::Fixed);
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Level {
    /// Not determined.
    Free,
    /// Partly determined.
    Restricted,
    /// Fully determined.
    Fixed,
}

/// The determinacy of a slot: the [`Level`] of its identity and of its
/// geometry.
///
/// The identity level is [`Identity::level`]. The geometry level is
/// [`Level::Fixed`] for an anchored slot, [`Level::Restricted`] for a free
/// slot that a constraint names, and [`Level::Free`] for any other free
/// slot. Determinacies are compared level by level: one is less than another
/// when neither of its levels is greater and at least one is less, and two
/// determinacies where each has one greater level are incomparable.
///
/// # Examples
///
/// ```
/// use mmcif_ir::{Corner, Determinacy, Level};
///
/// let de_novo = Determinacy::new(Level::Free, Level::Free);
/// let fixed_backbone = Determinacy::new(Level::Free, Level::Fixed);
/// let flexible = Determinacy::new(Level::Fixed, Level::Free);
///
/// assert_eq!(fixed_backbone.corner(), Some(Corner::FixedBackbone));
/// assert!(de_novo < fixed_backbone);
/// assert_eq!(fixed_backbone.partial_cmp(&flexible), None);
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub struct Determinacy {
    /// The level of the identity.
    pub identity: Level,
    /// The level of the geometry.
    pub geometry: Level,
}

/// A determinacy whose identity and geometry are each free or fixed.
///
/// # Examples
///
/// ```
/// use mmcif_ir::{Corner, Determinacy, Level};
///
/// let context = Determinacy::new(Level::Fixed, Level::Fixed);
///
/// assert_eq!(context.corner(), Some(Corner::Context));
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum Corner {
    /// A free identity and a free geometry.
    DeNovo,
    /// A free identity and a fixed geometry.
    FixedBackbone,
    /// A fixed identity and a free geometry.
    Flexible,
    /// A fixed identity and a fixed geometry.
    Context,
}

impl Determinacy {
    /// Creates a determinacy from the level of the identity and the level of
    /// the geometry.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_ir::{Determinacy, Level};
    ///
    /// let determinacy = Determinacy::new(Level::Fixed, Level::Restricted);
    ///
    /// assert_eq!(determinacy.identity, Level::Fixed);
    /// assert_eq!(determinacy.geometry, Level::Restricted);
    /// ```
    pub fn new(identity: Level, geometry: Level) -> Determinacy {
        Determinacy { identity, geometry }
    }

    /// Returns the corner the determinacy sits at, or `None` if either level
    /// is [`Level::Restricted`].
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_ir::{Corner, Determinacy, Level};
    ///
    /// assert_eq!(
    ///     Determinacy::new(Level::Fixed, Level::Free).corner(),
    ///     Some(Corner::Flexible)
    /// );
    /// assert_eq!(Determinacy::new(Level::Fixed, Level::Restricted).corner(), None);
    /// ```
    pub fn corner(self) -> Option<Corner> {
        match (self.identity, self.geometry) {
            (Level::Free, Level::Free) => Some(Corner::DeNovo),
            (Level::Free, Level::Fixed) => Some(Corner::FixedBackbone),
            (Level::Fixed, Level::Free) => Some(Corner::Flexible),
            (Level::Fixed, Level::Fixed) => Some(Corner::Context),
            _ => None,
        }
    }

    /// Returns `true` if either level is [`Level::Restricted`].
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_ir::{Determinacy, Level};
    ///
    /// assert!(Determinacy::new(Level::Fixed, Level::Restricted).is_partial());
    /// assert!(!Determinacy::new(Level::Fixed, Level::Free).is_partial());
    /// ```
    pub fn is_partial(self) -> bool {
        self.identity == Level::Restricted || self.geometry == Level::Restricted
    }
}

impl PartialOrd for Determinacy {
    fn partial_cmp(&self, other: &Determinacy) -> Option<Ordering> {
        let identity = self.identity.cmp(&other.identity);
        let geometry = self.geometry.cmp(&other.geometry);
        match (identity, geometry) {
            (Ordering::Equal, order) | (order, Ordering::Equal) => Some(order),
            (first, second) if first == second => Some(first),
            _ => None,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use mmcif_model::{LabelResidueKey, LabelSeqId};

    fn comp(id: &str) -> CompId {
        CompId::new(id)
    }

    fn anchor() -> ResidueRef {
        ResidueRef::Label(LabelResidueKey {
            asym: "A".into(),
            seq: LabelSeqId::new(1).unwrap(),
        })
    }

    #[test]
    fn a_single_allowed_component_is_a_fixed_identity() {
        let one = Identity::one_of([comp("ALA")]).unwrap();
        assert_eq!(one, Identity::Exactly(comp("ALA")));
        assert_eq!(one.level(), Level::Fixed);
        let acidic = Identity::one_of([comp("ASP"), comp("GLU")]).unwrap();
        assert_eq!(acidic.level(), Level::Restricted);
        assert!(acidic.allows(&comp("ASP")) && !acidic.allows(&comp("LYS")));
        assert_eq!(Identity::one_of([]), Err(IrError::EmptyAlphabet));
    }

    #[test]
    fn a_hand_built_singleton_set_is_also_fixed() {
        let singleton = Identity::OneOf(AlphabetSet::new([comp("ALA")]).unwrap());
        assert_eq!(singleton.level(), Level::Fixed);
        assert_eq!(singleton.level(), Identity::Exactly(comp("ALA")).level());
    }

    #[test]
    fn repeated_components_collapse() {
        let allowed = AlphabetSet::new([comp("ASP"), comp("GLU"), comp("ASP")]).unwrap();
        assert_eq!(allowed.len(), 2);
        assert_eq!(
            Identity::one_of([comp("ALA"), comp("ALA")])
                .unwrap()
                .level(),
            Level::Fixed
        );
    }

    #[test]
    fn the_four_corners_come_from_the_extremes() {
        let corner = |slot: &Slot, constrained| slot.determinacy(constrained).corner();
        assert_eq!(corner(&Slot::de_novo(), false), Some(Corner::DeNovo));
        assert_eq!(
            corner(&Slot::context(comp("ALA"), anchor()), false),
            Some(Corner::Context)
        );
        let fixed_backbone = Slot::new(Identity::Any, Geometry::Anchored(anchor()));
        assert_eq!(corner(&fixed_backbone, false), Some(Corner::FixedBackbone));
        let flexible = Slot::new(Identity::Exactly(comp("ALA")), Geometry::Free);
        assert_eq!(corner(&flexible, false), Some(Corner::Flexible));
    }

    #[test]
    fn partial_states_sit_off_the_corners() {
        let restricted = Slot::new(
            Identity::one_of([comp("ASP"), comp("GLU")]).unwrap(),
            Geometry::Free,
        );
        let determinacy = restricted.determinacy(false);
        assert_eq!(determinacy.corner(), None);
        assert!(determinacy.is_partial());
        let bounded = Slot::de_novo().determinacy(true);
        assert_eq!(bounded.corner(), None);
        assert_eq!(bounded.geometry, Level::Restricted);
    }

    #[test]
    fn determinacy_orders_componentwise_and_leaves_mixed_pairs_incomparable() {
        let free = Determinacy {
            identity: Level::Free,
            geometry: Level::Free,
        };
        let context = Determinacy {
            identity: Level::Fixed,
            geometry: Level::Fixed,
        };
        let fixed_backbone = Determinacy {
            identity: Level::Free,
            geometry: Level::Fixed,
        };
        let flexible = Determinacy {
            identity: Level::Fixed,
            geometry: Level::Free,
        };
        assert!(free < context);
        assert!(free < fixed_backbone);
        assert!(fixed_backbone < context);
        assert_eq!(fixed_backbone.partial_cmp(&flexible), None);
        assert_eq!(free.partial_cmp(&free), Some(Ordering::Equal));
    }

    #[test]
    fn the_ordering_matches_the_componentwise_order_for_every_pair() {
        let levels = [Level::Free, Level::Restricted, Level::Fixed];
        for identity_left in levels {
            for geometry_left in levels {
                for identity_right in levels {
                    for geometry_right in levels {
                        let left = Determinacy::new(identity_left, geometry_left);
                        let right = Determinacy::new(identity_right, geometry_right);
                        let identity = identity_left.cmp(&identity_right);
                        let geometry = geometry_left.cmp(&geometry_right);
                        let expected = match (identity, geometry) {
                            (Ordering::Equal, order) | (order, Ordering::Equal) => Some(order),
                            (first, second) if first == second => Some(first),
                            _ => None,
                        };
                        assert_eq!(left.partial_cmp(&right), expected, "{left:?} {right:?}");
                        assert_eq!(
                            left == right,
                            left.partial_cmp(&right) == Some(Ordering::Equal)
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn anchoring_exposes_the_context_reference() {
        let slot = Slot::context(comp("ALA"), anchor());
        assert_eq!(slot.geometry.anchor(), Some(&anchor()));
        assert_eq!(Slot::de_novo().geometry.anchor(), None);
    }
}
