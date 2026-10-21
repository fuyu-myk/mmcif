use mmcif_ir::{Identity, Slot};
use mmcif_spec::{Design, SegmentSource};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Role {
    Motif,
    MotifFreeSequence,
    Designed,
    Restricted,
}

pub fn role(slot: &Slot) -> Role {
    match (slot.geometry.anchor().is_some(), &slot.identity) {
        (true, Identity::Any) => Role::MotifFreeSequence,
        (true, _) => Role::Motif,
        (false, Identity::Any) => Role::Designed,
        (false, _) => Role::Restricted,
    }
}

pub fn roles(design: &Design, segment: &SegmentSource) -> Vec<(Role, Option<u32>)> {
    let Ok(definition) = design.spec.segment(&segment.id) else {
        return Vec::new();
    };
    definition
        .elements()
        .iter()
        .enumerate()
        .map(|(position, element)| {
            (
                role(element.slot()),
                segment.elements.get(position).and_then(|info| info.label),
            )
        })
        .collect()
}

pub fn redesigned(design: &Design, segment: &SegmentSource) -> Vec<u32> {
    roles(design, segment)
        .into_iter()
        .filter(|(role, _)| *role == Role::Designed)
        .filter_map(|(_, label)| label)
        .collect()
}

pub fn motif_labels(design: &Design, segment: &SegmentSource) -> Vec<u32> {
    roles(design, segment)
        .into_iter()
        .filter(|(role, _)| matches!(role, Role::Motif | Role::MotifFreeSequence))
        .filter_map(|(_, label)| label)
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use mmcif_ir::{Geometry, ResidueRef};
    use mmcif_model::{CompId, LabelResidueKey, LabelSeqId};

    fn anchor() -> ResidueRef {
        ResidueRef::Label(LabelResidueKey {
            asym: "A".into(),
            seq: LabelSeqId::new(1).unwrap(),
        })
    }

    #[test]
    fn each_corner_maps_to_the_role_a_backend_needs() {
        assert_eq!(role(&Slot::de_novo()), Role::Designed);
        assert_eq!(
            role(&Slot::context(CompId::new("CYS"), anchor())),
            Role::Motif
        );
        assert_eq!(
            role(&Slot::new(Identity::Any, Geometry::Anchored(anchor()))),
            Role::MotifFreeSequence
        );
        assert_eq!(
            role(&Slot::new(
                Identity::Exactly(CompId::new("CYS")),
                Geometry::Free
            )),
            Role::Restricted
        );
    }
}
