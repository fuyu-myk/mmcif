use mmcif_ir::{
    AtomAddress, Bounds, Constraint, Determinacy, Element, ElementIndex, Geometry, GroupId,
    Identity, Instantiated, IrError, Length, Lengths, Level, ResidueRef, Scope, SecondaryStructure,
    Segment, SegmentId, Slot, SlotIndex, Spec, Symbolic, SymbolicAddress,
};
use mmcif_model::{AtomName, CompId, LabelResidueKey, LabelSeqId};
use proptest::prelude::*;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum IdentityShape {
    Any,
    Restricted,
    Fixed,
}

#[derive(Debug, Clone, Copy)]
struct ElementShape {
    identity: IdentityShape,
    anchored: bool,
    min: u32,
    max: u32,
}

#[derive(Debug, Clone)]
struct Shape {
    segments: Vec<Vec<ElementShape>>,
    picks: Vec<(usize, usize, usize)>,
}

fn anchor() -> ResidueRef {
    ResidueRef::Label(LabelResidueKey {
        asym: "A".into(),
        seq: LabelSeqId::new(1).unwrap(),
    })
}

fn element_shape() -> impl Strategy<Value = ElementShape> {
    (
        prop_oneof![
            Just(IdentityShape::Any),
            Just(IdentityShape::Restricted),
            Just(IdentityShape::Fixed)
        ],
        any::<bool>(),
        0u32..4,
        1u32..5,
    )
        .prop_map(|(identity, anchored, low, high)| ElementShape {
            identity,
            anchored,
            min: low.min(high),
            max: high.max(low.min(high)),
        })
}

fn shape() -> impl Strategy<Value = Shape> {
    (
        prop::collection::vec(prop::collection::vec(element_shape(), 0..4), 1..4),
        prop::collection::vec((0usize..8, 0usize..8, 0usize..8), 0..5),
    )
        .prop_map(|(segments, picks)| Shape { segments, picks })
}

fn identity_of(shape: IdentityShape) -> Identity {
    match shape {
        IdentityShape::Any => Identity::Any,
        IdentityShape::Restricted => {
            Identity::one_of([CompId::new("ASP"), CompId::new("GLU")]).unwrap()
        }
        IdentityShape::Fixed => Identity::Exactly(CompId::new("ALA")),
    }
}

fn segment_id(index: usize) -> SegmentId {
    SegmentId::new(format!("s{index}"))
}

fn singular_address(shape: &Shape, pick: (usize, usize, usize)) -> Option<SymbolicAddress> {
    let occupied: Vec<usize> = (0..shape.segments.len())
        .filter(|&index| !shape.segments[index].is_empty())
        .collect();
    let segment = *occupied.get(pick.0 % occupied.len().max(1))?;
    let elements = &shape.segments[segment];
    let element = pick.1 % elements.len();
    let offset = pick.2 % elements[element].max as usize;
    Some(SymbolicAddress::at(
        segment_id(segment),
        ElementIndex::new(element as u32),
        offset as u32,
    ))
}

fn whole_address(shape: &Shape, pick: (usize, usize, usize)) -> Option<SymbolicAddress> {
    let address = singular_address(shape, pick)?;
    Some(SymbolicAddress::whole(address.segment, address.element))
}

#[derive(Debug)]
struct Built {
    spec: Spec,
    required: Vec<(usize, usize, u32)>,
}

fn build(shape: &Shape) -> Built {
    let mut spec = Spec::new();
    let mut required = Vec::new();
    for (index, elements) in shape.segments.iter().enumerate() {
        let elements = elements.iter().map(|element| {
            let geometry = if element.anchored {
                Geometry::Anchored(anchor())
            } else {
                Geometry::Free
            };
            Element::run(
                Slot::new(identity_of(element.identity), geometry),
                Length::range(element.min, element.max).unwrap(),
            )
        });
        spec.push_segment(Segment::new(segment_id(index), elements))
            .unwrap();
    }

    let mut members = Vec::new();
    for (position, pick) in shape.picks.iter().enumerate() {
        match position % 3 {
            0 => {
                if let (Some(first), Some(second)) = (
                    singular_address(shape, *pick),
                    singular_address(shape, (pick.2, pick.0, pick.1)),
                ) {
                    for address in [&first, &second] {
                        let segment = spec.segment_index(&address.segment).unwrap();
                        let offset = address.offset.expect("distance addresses carry an offset");
                        required.push((segment.as_usize(), address.element.as_usize(), offset + 1));
                    }
                    spec.push_constraint(Constraint::<Symbolic>::Distance {
                        first: AtomAddress::new(first, AtomName::new("CA")),
                        second: AtomAddress::new(second, AtomName::new("CB")),
                        bounds: Bounds::new(2.0, 8.0).unwrap(),
                    })
                    .unwrap();
                }
            }
            1 => {
                if let Some(address) = whole_address(shape, *pick) {
                    spec.push_constraint(Constraint::<Symbolic>::SecondaryStructure {
                        scope: Scope::Address(address),
                        kind: SecondaryStructure::Helix,
                    })
                    .unwrap();
                }
            }
            _ => members.extend(whole_address(shape, *pick)),
        }
    }
    if !members.is_empty() {
        spec.define_group(GroupId::new("group"), members).unwrap();
        spec.push_constraint(Constraint::<Symbolic>::SecondaryStructure {
            scope: Scope::Group(GroupId::new("group")),
            kind: SecondaryStructure::Strand,
        })
        .unwrap();
    }
    Built { spec, required }
}

fn without_constraints(shape: &Shape) -> Spec {
    build(&Shape {
        segments: shape.segments.clone(),
        picks: Vec::new(),
    })
    .spec
}

fn chosen_lengths(built: &Built, picks: &[u32]) -> Lengths {
    let counts: Vec<u32> = built
        .spec
        .symbolic_elements()
        .enumerate()
        .map(|(position, (segment, element, length))| {
            let span = length.max() - length.min() + 1;
            let pick = picks
                .get(position % picks.len().max(1))
                .copied()
                .unwrap_or(0);
            let chosen = length.min() + pick % span;
            let needed = built
                .required
                .iter()
                .filter(|(addressed_segment, addressed_element, _)| {
                    *addressed_segment == segment.as_usize()
                        && *addressed_element == element.as_usize()
                })
                .map(|(_, _, slots)| *slots)
                .max()
                .unwrap_or(0);
            chosen.max(needed)
        })
        .collect();
    Lengths::new(&built.spec, counts).expect("counts are inside their ranges")
}

fn expected_constrained(spec: &Spec, instantiated: &Instantiated) -> Vec<bool> {
    let mut touched: Vec<SlotIndex> = Vec::new();
    for constraint in spec.constraints() {
        match constraint {
            Constraint::Distance { first, second, .. } => {
                touched.extend(instantiated.resolve(&first.slot).unwrap());
                touched.extend(instantiated.resolve(&second.slot).unwrap());
            }
            Constraint::SecondaryStructure { scope, .. } => match scope {
                Scope::Address(address) => touched.extend(instantiated.resolve(address).unwrap()),
                Scope::Group(id) => {
                    for member in spec.group(id).unwrap() {
                        touched.extend(instantiated.resolve(member).unwrap());
                    }
                }
                _ => {}
            },
            _ => {}
        }
    }
    let mut constrained = vec![false; instantiated.slot_count()];
    for slot in touched {
        constrained[slot.as_usize()] = true;
    }
    constrained
}

fn expected_determinacy(
    shape: &Shape,
    instantiated: &Instantiated,
    index: SlotIndex,
    constrained: bool,
) -> Determinacy {
    let origin = instantiated.origin(index).unwrap();
    let element = shape.segments[origin.segment.as_usize()][origin.element.as_usize()];
    let identity = match element.identity {
        IdentityShape::Any => Level::Free,
        IdentityShape::Restricted => Level::Restricted,
        IdentityShape::Fixed => Level::Fixed,
    };
    let geometry = match (element.anchored, constrained) {
        (true, _) => Level::Fixed,
        (false, true) => Level::Restricted,
        (false, false) => Level::Free,
    };
    Determinacy::new(identity, geometry)
}

fn small_element() -> impl Strategy<Value = ElementShape> {
    (
        prop_oneof![
            Just(IdentityShape::Any),
            Just(IdentityShape::Restricted),
            Just(IdentityShape::Fixed)
        ],
        any::<bool>(),
        0u32..3,
        0u32..3,
    )
        .prop_map(|(identity, anchored, min, span)| ElementShape {
            identity,
            anchored,
            min,
            max: (min + span).max(1),
        })
}

fn small_shape() -> impl Strategy<Value = Shape> {
    prop::collection::vec(prop::collection::vec(small_element(), 0..3), 1..3).prop_map(|segments| {
        Shape {
            segments,
            picks: Vec::new(),
        }
    })
}

fn shape_with_a_symbolic_element() -> impl Strategy<Value = Shape> {
    shape().prop_map(|mut shape| {
        let symbolic = ElementShape {
            identity: IdentityShape::Any,
            anchored: false,
            min: 1,
            max: 3,
        };
        match shape.segments.first_mut() {
            Some(segment) => segment.push(symbolic),
            None => shape.segments.push(vec![symbolic]),
        }
        shape
    })
}

proptest! {
    #[test]
    fn minimum_and_maximum_lengths_hit_the_declared_bounds(shape in shape()) {
        let spec = without_constraints(&shape);
        let (min, max) = spec.length_bounds();
        prop_assert_eq!(spec.instantiate(&Lengths::minimum(&spec)).unwrap().slot_count() as u64, min);
        prop_assert_eq!(spec.instantiate(&Lengths::maximum(&spec)).unwrap().slot_count() as u64, max);
    }

    #[test]
    fn every_slot_resolves_back_to_itself(shape in shape(), picks in prop::collection::vec(0u32..6, 0..8)) {
        let built = build(&shape);
        let spec = &built.spec;
        let instantiated = spec.instantiate(&chosen_lengths(&built, &picks)).unwrap();
        for index in 0..instantiated.slot_count() as u32 {
            let index = SlotIndex::new(index);
            let source = instantiated.source(index).unwrap();
            prop_assert_eq!(instantiated.resolve(&source).unwrap(), vec![index]);
            let origin = instantiated.origin(index).unwrap();
            let segment = &shape.segments[origin.segment.as_usize()];
            prop_assert!(origin.element.as_usize() < segment.len());
            prop_assert!(origin.offset < segment[origin.element.as_usize()].max);
        }
    }

    #[test]
    fn segments_partition_the_slots_in_order(shape in shape(), picks in prop::collection::vec(0u32..6, 0..8)) {
        let built = build(&shape);
        let spec = &built.spec;
        let instantiated = spec.instantiate(&chosen_lengths(&built, &picks)).unwrap();
        let mut next = 0;
        for id in instantiated.segment_ids() {
            let span = instantiated.segment_slots(id).unwrap();
            prop_assert_eq!(span.start().get(), next);
            next = span.end().get();
        }
        prop_assert_eq!(next as usize, instantiated.slot_count());
    }

    #[test]
    fn instantiation_is_deterministic(shape in shape(), picks in prop::collection::vec(0u32..6, 0..8)) {
        let built = build(&shape);
        let lengths = chosen_lengths(&built, &picks);
        prop_assert_eq!(
            built.spec.instantiate(&lengths).unwrap(),
            built.spec.instantiate(&lengths).unwrap()
        );
    }

    #[test]
    fn constrained_slots_are_exactly_the_ones_constraints_name(
        shape in shape(),
        picks in prop::collection::vec(0u32..6, 0..8),
    ) {
        let built = build(&shape);
        let spec = &built.spec;
        let instantiated = spec.instantiate(&chosen_lengths(&built, &picks)).unwrap();
        let expected = expected_constrained(spec, &instantiated);
        for index in 0..instantiated.slot_count() as u32 {
            let index = SlotIndex::new(index);
            prop_assert_eq!(instantiated.is_constrained(index).unwrap(), expected[index.as_usize()]);
        }
        for slot in instantiated.constraints().iter().flat_map(Constraint::slots) {
            prop_assert!(slot.as_usize() < instantiated.slot_count());
        }
    }

    #[test]
    fn resolved_scopes_and_groups_are_sorted_without_repeats(
        shape in shape(),
        picks in prop::collection::vec(0u32..6, 0..8),
    ) {
        let built = build(&shape);
        let instantiated = built.spec.instantiate(&chosen_lengths(&built, &picks)).unwrap();
        let tidy = |slots: &[SlotIndex]| {
            let mut expected = slots.to_vec();
            expected.sort_unstable();
            expected.dedup();
            expected
        };
        for constraint in instantiated.constraints() {
            if let Constraint::SecondaryStructure { scope, .. } = constraint {
                prop_assert_eq!(scope, &tidy(scope));
            }
        }
        for (_, slots) in instantiated.groups() {
            prop_assert_eq!(slots, tidy(slots));
        }
    }

    #[test]
    fn determinacy_follows_the_slot_and_its_constraints(
        shape in shape(),
        picks in prop::collection::vec(0u32..6, 0..8),
    ) {
        let built = build(&shape);
        let spec = &built.spec;
        let instantiated = spec.instantiate(&chosen_lengths(&built, &picks)).unwrap();
        let constrained = expected_constrained(spec, &instantiated);
        for index in 0..instantiated.slot_count() as u32 {
            let index = SlotIndex::new(index);
            let determinacy = instantiated.determinacy(index).unwrap();
            prop_assert_eq!(
                determinacy,
                expected_determinacy(&shape, &instantiated, index, constrained[index.as_usize()])
            );
            prop_assert_eq!(determinacy.corner().is_some(), !determinacy.is_partial());
            let source = instantiated.source(index).unwrap();
            prop_assert_eq!(spec.determinacy(&source).unwrap(), determinacy);
        }
    }

    #[test]
    fn every_valid_length_choice_lands_inside_the_bounds(shape in small_shape()) {
        let spec = without_constraints(&shape);
        let ranges: Vec<(u32, u32)> = spec
            .symbolic_elements()
            .map(|(_, _, length)| (length.min(), length.max()))
            .collect();
        let combinations: u64 = ranges
            .iter()
            .map(|(min, max)| u64::from(max - min + 1))
            .product();
        let (min, max) = spec.length_bounds();
        let mut seen = Vec::new();
        for combination in 0..combinations {
            let mut rest = combination;
            let counts: Vec<u32> = ranges
                .iter()
                .map(|(low, high)| {
                    let span = u64::from(high - low + 1);
                    let value = low + (rest % span) as u32;
                    rest /= span;
                    value
                })
                .collect();
            let lengths = Lengths::new(&spec, counts).unwrap();
            let total = spec.instantiate(&lengths).unwrap().slot_count() as u64;
            prop_assert!((min..=max).contains(&total));
            seen.push(total);
        }
        prop_assert_eq!(seen.iter().copied().min(), Some(min));
        prop_assert_eq!(seen.iter().copied().max(), Some(max));
    }

    #[test]
    fn counts_outside_their_range_are_rejected(
        shape in shape_with_a_symbolic_element(),
        excess in 1u32..5,
    ) {
        let spec = without_constraints(&shape);
        let counts: Vec<u32> = spec
            .symbolic_elements()
            .map(|(_, _, length)| length.max() + excess)
            .collect();
        let rejected = Lengths::new(&spec, counts);
        prop_assert!(
            matches!(rejected, Err(IrError::LengthOutOfRange { .. })),
            "{:?}",
            rejected
        );
    }

    #[test]
    fn narrowing_either_axis_never_lowers_determinacy(
        identity in prop::sample::select(vec![Level::Free, Level::Restricted, Level::Fixed]),
        geometry in prop::sample::select(vec![Level::Free, Level::Restricted, Level::Fixed]),
        narrower in prop::sample::select(vec![Level::Free, Level::Restricted, Level::Fixed]),
    ) {
        let start = Determinacy::new(identity, geometry);
        let narrowed = Determinacy::new(identity.max(narrower), geometry);
        prop_assert!(start <= narrowed);
        let both = Determinacy::new(identity.max(narrower), geometry.max(narrower));
        prop_assert!(start <= both);
        prop_assert_eq!(start.corner().is_none(), start.is_partial());
    }
}
