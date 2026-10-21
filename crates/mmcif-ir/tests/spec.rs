use mmcif_ir::{
    AtomAddress, Bounds, Constraint, Corner, Element, Geometry, GroupId, Identity, IrError, Length,
    Lengths, Level, ResidueRef, Scope, SecondaryStructure, Segment, SegmentId, Slot, SlotIndex,
    Spec, Symbolic, SymbolicAddress, TokenLayout,
};
use mmcif_model::{AtomName, CompId, LabelResidueKey, LabelSeqId};

fn comp(id: &str) -> CompId {
    CompId::new(id)
}

fn anchor(seq: u32) -> ResidueRef {
    ResidueRef::Label(LabelResidueKey {
        asym: "A".into(),
        seq: LabelSeqId::new(seq).unwrap(),
    })
}

fn context_run(count: u32, first_seq: u32) -> Vec<Element> {
    (0..count)
        .map(|offset| Element::Slot(Slot::context(comp("GLY"), anchor(first_seq + offset))))
        .collect()
}

fn address(segment: &str, element: u32) -> SymbolicAddress {
    SymbolicAddress::whole(
        SegmentId::new(segment),
        mmcif_ir::ElementIndex::new(element),
    )
}

fn at(segment: &str, element: u32, offset: u32) -> SymbolicAddress {
    SymbolicAddress::at(
        SegmentId::new(segment),
        mmcif_ir::ElementIndex::new(element),
        offset,
    )
}

fn binder_against_target() -> Spec {
    let mut spec = Spec::new();
    spec.push_segment(Segment::new(SegmentId::new("target"), context_run(4, 10)))
        .unwrap();
    spec.push_segment(Segment::new(
        SegmentId::new("binder"),
        [Element::run(Slot::de_novo(), Length::range(8, 12).unwrap())],
    ))
    .unwrap();
    spec
}

#[test]
fn a_de_novo_binder_sits_beside_a_fixed_target() {
    let spec = binder_against_target();
    assert_eq!(spec.length_bounds(), (12, 16));
    assert!(spec.is_symbolic());
    assert_eq!(
        spec.determinacy(&address("target", 0)).unwrap().corner(),
        Some(Corner::Context)
    );
    assert_eq!(
        spec.determinacy(&address("binder", 0)).unwrap().corner(),
        Some(Corner::DeNovo)
    );

    let instantiated = spec.instantiate(&Lengths::minimum(&spec)).unwrap();
    assert_eq!(instantiated.slot_count(), 12);
    let target = instantiated
        .segment_slots(&SegmentId::new("target"))
        .unwrap();
    let binder = instantiated
        .segment_slots(&SegmentId::new("binder"))
        .unwrap();
    assert_eq!((target.start().get(), target.len()), (0, 4));
    assert_eq!((binder.start().get(), binder.len()), (4, 8));
    assert!(binder.contains(SlotIndex::new(11)) && !binder.contains(SlotIndex::new(12)));
    let corners: Vec<Option<Corner>> = (0..instantiated.slot_count() as u32)
        .map(|index| {
            instantiated
                .determinacy(SlotIndex::new(index))
                .unwrap()
                .corner()
        })
        .collect();
    assert_eq!(corners[..4], [Some(Corner::Context); 4]);
    assert_eq!(corners[4..], [Some(Corner::DeNovo); 8]);
}

#[test]
fn cdr_slots_sit_inside_a_context_framework() {
    let mut spec = Spec::new();
    let mut elements = context_run(3, 1);
    elements.push(Element::run(Slot::de_novo(), Length::range(5, 7).unwrap()));
    elements.extend(context_run(2, 20));
    spec.push_segment(Segment::new(SegmentId::new("heavy"), elements))
        .unwrap();

    let instantiated = spec.instantiate(&Lengths::maximum(&spec)).unwrap();
    assert_eq!(instantiated.slot_count(), 12);
    let de_novo: Vec<u32> = (0..instantiated.slot_count() as u32)
        .filter(|&index| {
            instantiated
                .determinacy(SlotIndex::new(index))
                .unwrap()
                .corner()
                == Some(Corner::DeNovo)
        })
        .collect();
    assert_eq!(de_novo, [3, 4, 5, 6, 7, 8, 9]);
}

#[test]
fn a_motif_lands_at_an_index_the_spec_does_not_fix() {
    let mut spec = Spec::new();
    spec.push_segment(Segment::new(
        SegmentId::new("scaffold"),
        [
            Element::run(Slot::de_novo(), Length::range(2, 5).unwrap()),
            Element::Slot(Slot::context(comp("HIS"), anchor(57))),
            Element::run(Slot::de_novo(), Length::range(1, 3).unwrap()),
        ],
    ))
    .unwrap();

    let motif_index = |lengths: &Lengths| {
        let instantiated = spec.instantiate(lengths).unwrap();
        let found: Vec<u32> = (0..instantiated.slot_count() as u32)
            .filter(|&index| {
                instantiated
                    .slot(SlotIndex::new(index))
                    .unwrap()
                    .geometry
                    .anchor()
                    .is_some()
            })
            .collect();
        assert_eq!(found.len(), 1);
        (found[0], instantiated)
    };

    let (shortest, short_spec) = motif_index(&Lengths::minimum(&spec));
    let (longest, _) = motif_index(&Lengths::maximum(&spec));
    assert_eq!((shortest, longest), (2, 5));
    assert_eq!(
        short_spec.source(SlotIndex::new(shortest)).unwrap(),
        at("scaffold", 1, 0)
    );
}

fn disulfide_spec() -> Spec {
    let mut spec = Spec::new();
    let cysteine = Slot::new(Identity::Exactly(comp("CYS")), Geometry::Free);
    spec.push_segment(Segment::new(
        SegmentId::new("loop"),
        [
            Element::run(Slot::de_novo(), Length::range(3, 5).unwrap()),
            Element::Slot(cysteine.clone()),
            Element::run(Slot::de_novo(), Length::exact(2).unwrap()),
            Element::Slot(cysteine),
            Element::run(Slot::de_novo(), Length::range(1, 4).unwrap()),
        ],
    ))
    .unwrap();
    spec.push_constraint(Constraint::<Symbolic>::Distance {
        first: AtomAddress::new(address("loop", 1), AtomName::new("SG")),
        second: AtomAddress::new(address("loop", 3), AtomName::new("SG")),
        bounds: Bounds::new(1.8, 2.2).unwrap(),
    })
    .unwrap();
    spec.define_group(
        GroupId::new("stapled"),
        [address("loop", 1), address("loop", 3)],
    )
    .unwrap();
    spec.push_constraint(Constraint::<Symbolic>::SecondaryStructure {
        scope: Scope::Group(GroupId::new("stapled")),
        kind: SecondaryStructure::Helix,
    })
    .unwrap();
    spec
}

#[test]
fn a_disulfide_pattern_resolves_to_the_right_slots() {
    let spec = disulfide_spec();
    let lengths = Lengths::new(&spec, [4, 2]).unwrap();
    let instantiated = spec.instantiate(&lengths).unwrap();
    assert_eq!(instantiated.slot_count(), 4 + 1 + 2 + 1 + 2);

    let cysteines: Vec<SlotIndex> = instantiated
        .group(&GroupId::new("stapled"))
        .unwrap()
        .to_vec();
    assert_eq!(cysteines, [SlotIndex::new(4), SlotIndex::new(7)]);
    for slot in &cysteines {
        assert_eq!(
            instantiated.slot(*slot).unwrap().identity,
            Identity::Exactly(comp("CYS"))
        );
    }

    let resolved = instantiated.constraints();
    assert_eq!(resolved.len(), 2);
    match &resolved[0] {
        Constraint::Distance {
            first,
            second,
            bounds,
        } => {
            assert_eq!(
                (first.slot, second.slot),
                (SlotIndex::new(4), SlotIndex::new(7))
            );
            assert_eq!(first.atom, AtomName::new("SG"));
            assert!(bounds.contains(2.0));
        }
        other => panic!("expected a distance constraint, got {other:?}"),
    }
    match &resolved[1] {
        Constraint::SecondaryStructure { scope, kind } => {
            assert_eq!(scope, &cysteines);
            assert_eq!(*kind, SecondaryStructure::Helix);
        }
        other => panic!("expected a secondary structure constraint, got {other:?}"),
    }
}

#[test]
fn constrained_geometry_is_restricted_rather_than_free() {
    let spec = disulfide_spec();
    let instantiated = spec.instantiate(&Lengths::minimum(&spec)).unwrap();
    let cysteine = instantiated.group(&GroupId::new("stapled")).unwrap()[0];
    let determinacy = instantiated.determinacy(cysteine).unwrap();
    assert_eq!(determinacy.identity, Level::Fixed);
    assert_eq!(determinacy.geometry, Level::Restricted);
    assert_eq!(determinacy.corner(), None);
    assert!(determinacy.is_partial());

    let unconstrained = instantiated.determinacy(SlotIndex::new(0)).unwrap();
    assert_eq!(unconstrained.corner(), Some(Corner::DeNovo));
    assert!(!instantiated.is_constrained(SlotIndex::new(0)).unwrap());
    assert!(instantiated.is_constrained(cysteine).unwrap());
}

#[test]
fn a_restricted_identity_is_neither_free_nor_fixed() {
    let mut spec = Spec::new();
    let acidic = Identity::one_of([comp("ASP"), comp("GLU")]).unwrap();
    spec.push_segment(Segment::new(
        SegmentId::new("patch"),
        [Element::Slot(Slot::new(acidic, Geometry::Free))],
    ))
    .unwrap();
    let determinacy = spec.determinacy(&address("patch", 0)).unwrap();
    assert_eq!(determinacy.identity, Level::Restricted);
    assert_eq!(determinacy.corner(), None);
    let instantiated = spec.instantiate(&Lengths::minimum(&spec)).unwrap();
    let slot = &instantiated.slots()[0];
    assert!(slot.identity.allows(&comp("ASP")));
    assert!(!slot.identity.allows(&comp("ALA")));
}

#[test]
fn every_slot_maps_back_to_its_source_element() {
    let spec = disulfide_spec();
    let instantiated = spec
        .instantiate(&Lengths::new(&spec, [5, 4]).unwrap())
        .unwrap();
    for index in 0..instantiated.slot_count() as u32 {
        let index = SlotIndex::new(index);
        let source = instantiated.source(index).unwrap();
        assert_eq!(instantiated.resolve(&source).unwrap(), [index]);
        let origin = instantiated.origin(index).unwrap();
        let element = spec
            .segment(&SegmentId::new("loop"))
            .unwrap()
            .element(origin.element)
            .unwrap();
        assert_eq!(instantiated.slot(index).unwrap(), element.slot());
    }
}

#[test]
fn a_whole_element_address_resolves_to_every_slot_it_produced() {
    let spec = disulfide_spec();
    let instantiated = spec
        .instantiate(&Lengths::new(&spec, [3, 4].to_vec()).unwrap())
        .unwrap();
    let run = instantiated.resolve(&address("loop", 0)).unwrap();
    assert_eq!(run, (0..3).map(SlotIndex::new).collect::<Vec<_>>());
    let tail = instantiated.resolve(&address("loop", 4)).unwrap();
    assert_eq!(tail.len(), 4);
    assert_eq!(instantiated.resolve(&at("loop", 4, 3)).unwrap(), [tail[3]]);
    assert!(matches!(
        instantiated.resolve(&at("loop", 0, 3)),
        Err(IrError::OffsetOutOfRange { offset: 3, .. })
    ));
}

#[test]
fn a_zero_length_run_produces_no_slots() {
    let mut spec = Spec::new();
    spec.push_segment(Segment::new(
        SegmentId::new("linker"),
        [
            Element::Slot(Slot::context(comp("GLY"), anchor(1))),
            Element::run(Slot::de_novo(), Length::range(0, 4).unwrap()),
        ],
    ))
    .unwrap();
    let instantiated = spec.instantiate(&Lengths::minimum(&spec)).unwrap();
    assert_eq!(instantiated.slot_count(), 1);
    assert!(
        instantiated
            .resolve(&address("linker", 1))
            .unwrap()
            .is_empty()
    );
    assert!(matches!(
        instantiated.resolve(&at("linker", 1, 0)),
        Err(IrError::OffsetOutOfRange { .. })
    ));
}

#[test]
fn structural_mistakes_are_rejected_when_they_are_declared() {
    let mut spec = binder_against_target();
    assert_eq!(
        spec.push_segment(Segment::new(SegmentId::new("binder"), [])),
        Err(IrError::DuplicateSegment(SegmentId::new("binder")))
    );
    assert_eq!(
        spec.determinacy(&address("missing", 0)),
        Err(IrError::UnknownSegment(SegmentId::new("missing")))
    );
    assert!(matches!(
        spec.determinacy(&address("binder", 7)),
        Err(IrError::UnknownElement { .. })
    ));
    assert!(matches!(
        spec.define_group(GroupId::new("bad"), [at("binder", 0, 99)]),
        Err(IrError::OffsetOutOfRange { offset: 99, .. })
    ));
    assert!(matches!(
        spec.push_constraint(Constraint::<Symbolic>::SecondaryStructure {
            scope: Scope::Group(GroupId::new("undefined")),
            kind: SecondaryStructure::Helix,
        }),
        Err(IrError::UnknownGroup(_))
    ));
}

#[test]
fn a_distance_constraint_needs_an_address_that_names_one_slot() {
    let mut spec = binder_against_target();
    let whole_run = Constraint::<Symbolic>::Distance {
        first: AtomAddress::new(address("binder", 0), AtomName::new("CA")),
        second: AtomAddress::new(address("target", 0), AtomName::new("CA")),
        bounds: Bounds::at_most(8.0).unwrap(),
    };
    assert!(matches!(
        spec.push_constraint(whole_run),
        Err(IrError::AddressNotSingular { slots: 12, .. })
    ));
    let one_slot = Constraint::<Symbolic>::Distance {
        first: AtomAddress::new(at("binder", 0, 3), AtomName::new("CA")),
        second: AtomAddress::new(address("target", 0), AtomName::new("CA")),
        bounds: Bounds::at_most(8.0).unwrap(),
    };
    spec.push_constraint(one_slot).unwrap();
    assert_eq!(spec.constraints().len(), 1);
}

#[test]
fn chosen_lengths_are_checked_against_their_ranges() {
    let spec = disulfide_spec();
    assert_eq!(
        Lengths::new(&spec, [4]),
        Err(IrError::LengthCountMismatch {
            expected: 2,
            actual: 1
        })
    );
    assert_eq!(
        Lengths::new(&spec, [9, 2]),
        Err(IrError::LengthOutOfRange {
            chosen: 9,
            min: 3,
            max: 5
        })
    );
    assert_eq!(Lengths::minimum(&spec).counts(), [3, 1]);
    assert_eq!(Lengths::maximum(&spec).counts(), [5, 4]);
}

#[test]
fn a_token_layout_maps_only_the_slots_it_covers() {
    let spec = binder_against_target();
    let instantiated = spec.instantiate(&Lengths::minimum(&spec)).unwrap();
    let binder = instantiated
        .segment_slots(&SegmentId::new("binder"))
        .unwrap();
    let layout = TokenLayout::new(instantiated.slot_count(), binder.iter()).unwrap();
    assert_eq!(layout.token_count(), binder.len());
    assert_eq!(layout.token_for(binder.start()).unwrap().get(), 0);
    assert_eq!(
        layout.token_for(SlotIndex::new(0)),
        Err(IrError::SlotNotTokenized(SlotIndex::new(0)))
    );
}

#[test]
fn an_offset_past_the_instantiated_length_fails_at_instantiation() {
    let mut spec = Spec::new();
    spec.push_segment(Segment::new(
        SegmentId::new("loop"),
        [Element::run(Slot::de_novo(), Length::range(1, 6).unwrap())],
    ))
    .unwrap();
    spec.push_constraint(Constraint::<Symbolic>::Distance {
        first: AtomAddress::new(at("loop", 0, 0), AtomName::new("CA")),
        second: AtomAddress::new(at("loop", 0, 4), AtomName::new("CA")),
        bounds: Bounds::at_most(10.0).unwrap(),
    })
    .unwrap();

    let too_short = spec.instantiate(&Lengths::new(&spec, [3]).unwrap());
    assert!(matches!(
        too_short,
        Err(IrError::OffsetOutOfRange { offset: 4, .. })
    ));
    let long_enough = spec
        .instantiate(&Lengths::new(&spec, [5]).unwrap())
        .unwrap();
    assert_eq!(long_enough.constraints().len(), 1);
    assert_eq!(
        long_enough.constraints()[0].slots(),
        [SlotIndex::new(0), SlotIndex::new(4)]
    );
}

#[test]
fn symbolic_determinacy_follows_the_offset_the_constraint_names() {
    let mut spec = Spec::new();
    spec.push_segment(Segment::new(
        SegmentId::new("binder"),
        [Element::run(Slot::de_novo(), Length::exact(10).unwrap())],
    ))
    .unwrap();
    spec.push_constraint(Constraint::<Symbolic>::Distance {
        first: AtomAddress::new(at("binder", 0, 3), AtomName::new("CA")),
        second: AtomAddress::new(at("binder", 0, 4), AtomName::new("CA")),
        bounds: Bounds::at_most(8.0).unwrap(),
    })
    .unwrap();

    assert_eq!(
        spec.determinacy(&at("binder", 0, 3)).unwrap().geometry,
        Level::Restricted
    );
    assert_eq!(
        spec.determinacy(&at("binder", 0, 7)).unwrap().geometry,
        Level::Free
    );
    assert!(matches!(
        spec.determinacy(&address("binder", 0)),
        Err(IrError::MixedDeterminacy { .. })
    ));

    let instantiated = spec.instantiate(&Lengths::minimum(&spec)).unwrap();
    for index in 0..instantiated.slot_count() as u32 {
        let index = SlotIndex::new(index);
        let source = instantiated.source(index).unwrap();
        assert_eq!(
            spec.determinacy(&source).unwrap(),
            instantiated.determinacy(index).unwrap()
        );
    }
}

#[test]
fn a_constraint_over_a_whole_element_leaves_it_uniform() {
    let mut spec = Spec::new();
    spec.push_segment(Segment::new(
        SegmentId::new("helix"),
        [Element::run(Slot::de_novo(), Length::exact(6).unwrap())],
    ))
    .unwrap();
    spec.push_constraint(Constraint::<Symbolic>::SecondaryStructure {
        scope: Scope::Address(address("helix", 0)),
        kind: SecondaryStructure::Helix,
    })
    .unwrap();
    assert_eq!(
        spec.determinacy(&address("helix", 0)).unwrap().geometry,
        Level::Restricted
    );
    assert_eq!(
        spec.determinacy(&at("helix", 0, 5)).unwrap().geometry,
        Level::Restricted
    );
}

#[test]
fn length_choices_are_bound_to_the_spec_they_came_from() {
    let make = |elements: Vec<Element>| {
        let mut spec = Spec::new();
        spec.push_segment(Segment::new(SegmentId::new("a"), elements))
            .unwrap();
        spec
    };
    let run = || Element::run(Slot::de_novo(), Length::range(2, 9).unwrap());
    let fixed = || Element::Slot(Slot::de_novo());
    let run_first = make(vec![run(), fixed()]);
    let run_second = make(vec![fixed(), run()]);
    let lengths = Lengths::new(&run_first, [5]).unwrap();
    assert!(matches!(
        run_second.instantiate(&lengths),
        Err(IrError::LengthsFromAnotherSpec { .. })
    ));
    assert_eq!(run_first.instantiate(&lengths).unwrap().slot_count(), 6);
}

#[test]
fn an_absurd_specification_is_refused_before_allocating() {
    let mut spec = Spec::new();
    spec.push_segment(Segment::new(
        SegmentId::new("huge"),
        [Element::run(
            Slot::de_novo(),
            Length::range(1, u32::MAX).unwrap(),
        )],
    ))
    .unwrap();
    assert_eq!(
        spec.instantiate(&Lengths::maximum(&spec)),
        Err(IrError::TooManySlots(mmcif_ir::MAX_SLOTS))
    );
    assert!(spec.instantiate(&Lengths::minimum(&spec)).is_ok());
}

#[test]
fn a_singleton_restriction_reads_as_a_fixed_identity() {
    let mut spec = Spec::new();
    let one = Identity::OneOf(mmcif_ir::AlphabetSet::new([comp("ALA")]).unwrap());
    spec.push_segment(Segment::new(
        SegmentId::new("one"),
        [Element::Slot(Slot::new(one, Geometry::Free))],
    ))
    .unwrap();
    assert_eq!(
        spec.determinacy(&address("one", 0)).unwrap().corner(),
        Some(Corner::Flexible)
    );
}

#[test]
fn a_group_naming_the_same_slots_twice_resolves_them_once() {
    let mut spec = Spec::new();
    spec.push_segment(Segment::new(
        SegmentId::new("loop"),
        [
            Element::run(Slot::de_novo(), Length::exact(3).unwrap()),
            Element::Slot(Slot::de_novo()),
        ],
    ))
    .unwrap();
    spec.define_group(
        GroupId::new("twice"),
        [
            address("loop", 0),
            address("loop", 1),
            address("loop", 0),
            at("loop", 0, 1),
        ],
    )
    .unwrap();
    spec.push_constraint(Constraint::<Symbolic>::SecondaryStructure {
        scope: Scope::Group(GroupId::new("twice")),
        kind: SecondaryStructure::Strand,
    })
    .unwrap();

    let instantiated = spec.instantiate(&Lengths::minimum(&spec)).unwrap();
    let expected: Vec<SlotIndex> = (0..4).map(SlotIndex::new).collect();
    assert_eq!(
        instantiated.group(&GroupId::new("twice")).unwrap(),
        expected
    );
    match &instantiated.constraints()[0] {
        Constraint::SecondaryStructure { scope, .. } => assert_eq!(scope, &expected),
        other => panic!("expected a secondary structure constraint, got {other:?}"),
    }
}
