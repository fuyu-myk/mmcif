use std::cmp::Ordering;
use std::collections::BTreeSet;

use mmcif_model::CompId;

use crate::error::IrError;
use crate::ids::ResidueRef;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AlphabetSet(BTreeSet<CompId>);

impl AlphabetSet {
    pub fn new(components: impl IntoIterator<Item = CompId>) -> Result<AlphabetSet, IrError> {
        let components: BTreeSet<CompId> = components.into_iter().collect();
        if components.is_empty() {
            return Err(IrError::EmptyAlphabet);
        }
        Ok(AlphabetSet(components))
    }

    pub fn len(&self) -> usize {
        self.0.len()
    }

    pub fn is_empty(&self) -> bool {
        false
    }

    pub fn contains(&self, component: &CompId) -> bool {
        self.0.contains(component)
    }

    pub fn iter(&self) -> impl ExactSizeIterator<Item = &CompId> {
        self.0.iter()
    }

    pub fn is_subset(&self, other: &AlphabetSet) -> bool {
        self.0.is_subset(&other.0)
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum Identity {
    Any,
    OneOf(AlphabetSet),
    Exactly(CompId),
}

impl Identity {
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

    pub fn allows(&self, component: &CompId) -> bool {
        match self {
            Identity::Any => true,
            Identity::OneOf(allowed) => allowed.contains(component),
            Identity::Exactly(fixed) => fixed == component,
        }
    }

    pub fn level(&self) -> Level {
        match self {
            Identity::Any => Level::Free,
            Identity::OneOf(_) => Level::Restricted,
            Identity::Exactly(_) => Level::Fixed,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum Geometry {
    Free,
    Anchored(ResidueRef),
}

impl Geometry {
    pub fn anchor(&self) -> Option<&ResidueRef> {
        match self {
            Geometry::Free => None,
            Geometry::Anchored(reference) => Some(reference),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub struct Slot {
    pub identity: Identity,
    pub geometry: Geometry,
}

impl Slot {
    pub fn new(identity: Identity, geometry: Geometry) -> Slot {
        Slot { identity, geometry }
    }

    pub fn de_novo() -> Slot {
        Slot::new(Identity::Any, Geometry::Free)
    }

    pub fn context(component: CompId, anchor: ResidueRef) -> Slot {
        Slot::new(Identity::Exactly(component), Geometry::Anchored(anchor))
    }

    pub fn determinacy(&self, constrained: bool) -> Determinacy {
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

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Level {
    Free,
    Restricted,
    Fixed,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Determinacy {
    pub identity: Level,
    pub geometry: Level,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Corner {
    DeNovo,
    FixedBackbone,
    Flexible,
    Context,
}

impl Determinacy {
    pub fn corner(self) -> Option<Corner> {
        match (self.identity, self.geometry) {
            (Level::Free, Level::Free) => Some(Corner::DeNovo),
            (Level::Free, Level::Fixed) => Some(Corner::FixedBackbone),
            (Level::Fixed, Level::Free) => Some(Corner::Flexible),
            (Level::Fixed, Level::Fixed) => Some(Corner::Context),
            _ => None,
        }
    }

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
    fn anchoring_exposes_the_context_reference() {
        let slot = Slot::context(comp("ALA"), anchor());
        assert_eq!(slot.geometry.anchor(), Some(&anchor()));
        assert_eq!(Slot::de_novo().geometry.anchor(), None);
    }
}
