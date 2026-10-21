use std::collections::BTreeMap;

use mmcif_ccd::{AtomStereo, BondStereo, CcdIndex, Component, build_index};
use mmcif_cif::{Document, DocumentBuilder, Value};
use proptest::prelude::*;

const BOND_ORDERS: [&str; 8] = ["SING", "DOUB", "TRIP", "QUAD", "AROM", "DELO", "PI", "POLY"];
const ATOM_TAGS: [&str; 7] = [
    "_chem_comp_atom.comp_id",
    "_chem_comp_atom.atom_id",
    "_chem_comp_atom.type_symbol",
    "_chem_comp_atom.charge",
    "_chem_comp_atom.pdbx_aromatic_flag",
    "_chem_comp_atom.pdbx_leaving_atom_flag",
    "_chem_comp_atom.pdbx_stereo_config",
];
const BOND_TAGS: [&str; 6] = [
    "_chem_comp_bond.comp_id",
    "_chem_comp_bond.atom_id_1",
    "_chem_comp_bond.atom_id_2",
    "_chem_comp_bond.value_order",
    "_chem_comp_bond.pdbx_aromatic_flag",
    "_chem_comp_bond.pdbx_stereo_config",
];

#[derive(Debug, Clone)]
struct AtomSpec {
    prefix: String,
    element: String,
    charge: Option<i8>,
    aromatic: Option<bool>,
    leaving: Option<bool>,
    stereo: Option<&'static str>,
}

#[derive(Debug, Clone)]
struct BondSpec {
    atoms: [u16; 2],
    order: &'static str,
    aromatic: Option<bool>,
    stereo: Option<&'static str>,
}

#[derive(Debug, Clone)]
struct ComponentSpec {
    name: Option<String>,
    kind: Option<&'static str>,
    formula: Option<String>,
    parent: Option<String>,
    one_letter_code: Option<String>,
    formal_charge: Option<i32>,
    release_status: Option<&'static str>,
    atoms: Vec<AtomSpec>,
    bonds: Vec<BondSpec>,
}

impl ComponentSpec {
    fn atom_name(&self, index: usize) -> String {
        format!("{}{index}", self.atoms[index].prefix)
    }
}

fn atom_spec() -> impl Strategy<Value = AtomSpec> {
    (
        "[A-Z][A-Z0-9']{0,2}",
        "[A-Z][A-Z]?",
        proptest::option::of(-4i8..=4),
        proptest::option::of(any::<bool>()),
        proptest::option::of(any::<bool>()),
        proptest::option::of(prop::sample::select(vec!["R", "S", "N"])),
    )
        .prop_map(
            |(prefix, element, charge, aromatic, leaving, stereo)| AtomSpec {
                prefix,
                element,
                charge,
                aromatic,
                leaving,
                stereo,
            },
        )
}

fn bond_spec(atom_count: u16) -> impl Strategy<Value = BondSpec> {
    (
        [0..atom_count, 0..atom_count],
        prop::sample::select(BOND_ORDERS.to_vec()),
        proptest::option::of(any::<bool>()),
        proptest::option::of(prop::sample::select(vec!["E", "Z", "N"])),
    )
        .prop_map(|(atoms, order, aromatic, stereo)| BondSpec {
            atoms,
            order,
            aromatic,
            stereo,
        })
}

fn component_spec() -> impl Strategy<Value = ComponentSpec> {
    let header = (
        proptest::option::of("[ -~]{0,20}"),
        proptest::option::of(prop::sample::select(vec![
            "L-PEPTIDE LINKING",
            "NON-POLYMER",
            "D-saccharide, beta linking",
        ])),
        proptest::option::of("[A-Z][a-z]?[0-9]{0,2}( [A-Z][a-z]?[0-9]{0,2}){0,3}"),
        proptest::option::of("[A-Z]{3}"),
        proptest::option::of("[A-Z]{1,3}"),
        proptest::option::of(-3i32..=3),
        proptest::option::of(prop::sample::select(vec!["REL", "OBS", "HOLD"])),
    );
    (header, prop::collection::vec(atom_spec(), 0..8)).prop_flat_map(|(header, atoms)| {
        let bonds = match u16::try_from(atoms.len()) {
            Ok(count) if count > 0 => prop::collection::vec(bond_spec(count), 0..8).boxed(),
            _ => Just(Vec::new()).boxed(),
        };
        let (name, kind, formula, parent, one_letter_code, formal_charge, release_status) = header;
        bonds.prop_map(move |bonds| ComponentSpec {
            name: name.clone(),
            kind,
            formula: formula.clone(),
            parent: parent.clone(),
            one_letter_code: one_letter_code.clone(),
            formal_charge,
            release_status,
            atoms: atoms.clone(),
            bonds,
        })
    })
}

fn optional(text: Option<&str>) -> Value<'_> {
    text.map_or(Value::UNKNOWN, Value::text)
}

fn flag(flag: Option<bool>) -> Value<'static> {
    optional(flag.map(|flag| if flag { "Y" } else { "N" }))
}

fn write_component(builder: &mut DocumentBuilder, id: &str, spec: &ComponentSpec) {
    let charge = spec.formal_charge.map(|charge| charge.to_string());
    builder.start_block(id).unwrap();
    let pairs = [
        ("_chem_comp.id", Value::text(id)),
        ("_chem_comp.name", optional(spec.name.as_deref())),
        ("_chem_comp.type", optional(spec.kind)),
        ("_chem_comp.formula", optional(spec.formula.as_deref())),
        (
            "_chem_comp.mon_nstd_parent_comp_id",
            optional(spec.parent.as_deref()),
        ),
        (
            "_chem_comp.one_letter_code",
            optional(spec.one_letter_code.as_deref()),
        ),
        ("_chem_comp.pdbx_formal_charge", optional(charge.as_deref())),
        (
            "_chem_comp.pdbx_release_status",
            optional(spec.release_status),
        ),
    ];
    for (tag, value) in pairs {
        builder.pair(tag, value).unwrap();
    }
    if !spec.atoms.is_empty() {
        builder.start_loop(ATOM_TAGS).unwrap();
        for (index, atom) in spec.atoms.iter().enumerate() {
            let name = spec.atom_name(index);
            let charge = atom.charge.map(|charge| charge.to_string());
            let row = [
                Value::text(id),
                Value::text(&name),
                Value::text(&atom.element),
                optional(charge.as_deref()),
                flag(atom.aromatic),
                flag(atom.leaving),
                optional(atom.stereo),
            ];
            row.into_iter()
                .for_each(|value| builder.push(value).unwrap());
        }
    }
    if !spec.bonds.is_empty() {
        builder.start_loop(BOND_TAGS).unwrap();
        for bond in &spec.bonds {
            let [first, second] = bond.atoms.map(|index| spec.atom_name(usize::from(index)));
            let row = [
                Value::text(id),
                Value::text(&first),
                Value::text(&second),
                Value::text(bond.order),
                flag(bond.aromatic),
                optional(bond.stereo),
            ];
            row.into_iter()
                .for_each(|value| builder.push(value).unwrap());
        }
    }
}

fn render(specs: &BTreeMap<String, ComponentSpec>) -> String {
    let mut builder = DocumentBuilder::new();
    for (id, spec) in specs.iter().rev() {
        write_component(&mut builder, id, spec);
    }
    builder.finish().unwrap().to_cif_string()
}

fn atom_stereo(code: Option<&str>) -> AtomStereo {
    match code {
        Some("R") => AtomStereo::R,
        Some("S") => AtomStereo::S,
        _ => AtomStereo::None,
    }
}

fn bond_stereo(code: Option<&str>) -> BondStereo {
    match code {
        Some("E") => BondStereo::E,
        Some("Z") => BondStereo::Z,
        _ => BondStereo::None,
    }
}

fn assert_matches_spec(component: &Component, id: &str, spec: &ComponentSpec) {
    assert_eq!(component.id(), id);
    assert_eq!(component.name(), spec.name.as_deref().unwrap_or_default());
    assert_eq!(component.kind(), spec.kind.unwrap_or_default());
    assert_eq!(component.formula(), spec.formula.as_deref());
    assert_eq!(component.parent(), spec.parent.as_deref());
    assert_eq!(component.one_letter_code(), spec.one_letter_code.as_deref());
    assert_eq!(component.formal_charge(), spec.formal_charge.unwrap_or(0));
    assert_eq!(
        component.release_status(),
        spec.release_status.unwrap_or_default()
    );
    assert_eq!(component.atoms().len(), spec.atoms.len());
    for (index, (atom, expected)) in component.atoms().iter().zip(&spec.atoms).enumerate() {
        assert_eq!(atom.name(), spec.atom_name(index));
        assert_eq!(atom.element(), expected.element);
        assert_eq!(atom.charge(), expected.charge.unwrap_or(0));
        assert_eq!(atom.is_aromatic(), expected.aromatic.unwrap_or(false));
        assert_eq!(atom.is_leaving(), expected.leaving.unwrap_or(false));
        assert_eq!(atom.stereo(), atom_stereo(expected.stereo));
    }
    assert_eq!(component.bonds().len(), spec.bonds.len());
    for (bond, expected) in component.bonds().iter().zip(&spec.bonds) {
        assert_eq!(bond.atom_indices(), expected.atoms);
        assert_eq!(bond.order().code(), expected.order);
        assert_eq!(bond.is_aromatic(), expected.aromatic.unwrap_or(false));
        assert_eq!(bond.stereo(), bond_stereo(expected.stereo));
    }
}

proptest! {
    #[test]
    fn components_round_trip_through_the_index(
        specs in prop::collection::btree_map("[A-Z0-9]{1,8}", component_spec(), 1..6)
    ) {
        let text = render(&specs);
        let mut out = Vec::new();
        let summary = build_index(&text, "proptest", &mut out).unwrap();
        prop_assert_eq!(summary.components(), specs.len());
        let index = CcdIndex::from_bytes(out).unwrap();
        prop_assert_eq!(index.release(), "proptest");
        prop_assert!(index.ids().eq(specs.keys().map(String::as_str)));
        let document = Document::parse(text).unwrap();
        for block in document.blocks() {
            let id = block.name();
            let direct = Component::from_block(block).unwrap();
            assert_matches_spec(&direct, id, &specs[id]);
            prop_assert_eq!(index.get(id).unwrap(), Some(direct));
        }
    }
}
