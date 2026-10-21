use crate::codec::{Reader, put_opt_str, put_str, put_varint};
use crate::component::{
    AtomStereo, BondOrder, BondStereo, Component, ComponentAtom, ComponentBond, MAX_ATOMS,
};
use crate::error::{CcdError, corrupt};

const MIN_ATOM_LEN: usize = 4;
const MIN_BOND_LEN: usize = 6;
const AROMATIC_BIT: u8 = 1;
const LEAVING_BIT: u8 = 1 << 1;
const ATOM_STEREO_SHIFT: u8 = 2;
const BOND_STEREO_SHIFT: u8 = 1;

pub(crate) fn encode_record(component: &Component, out: &mut Vec<u8>) {
    put_str(out, &component.id);
    put_str(out, &component.name);
    put_str(out, &component.kind);
    put_opt_str(out, component.formula.as_deref());
    put_opt_str(out, component.parent.as_deref());
    out.push(u8::from(component.one_letter_code.is_some()));
    if let Some(code) = component.one_letter_code {
        out.extend_from_slice(&u32::from(code).to_le_bytes());
    }
    out.extend_from_slice(&component.formal_charge.to_le_bytes());
    put_str(out, &component.release_status);
    put_varint(out, component.atoms.len() as u64);
    for atom in &component.atoms {
        put_str(out, &atom.name);
        put_str(out, &atom.element);
        out.extend_from_slice(&atom.charge.to_le_bytes());
        out.push(atom_flags(atom));
    }
    put_varint(out, component.bonds.len() as u64);
    for bond in &component.bonds {
        out.extend_from_slice(&bond.atoms[0].to_le_bytes());
        out.extend_from_slice(&bond.atoms[1].to_le_bytes());
        out.push(order_code(bond.order));
        out.push(bond_flags(bond));
    }
}

pub(crate) fn decode_record(payload: &[u8]) -> Result<Component, CcdError> {
    let mut reader = Reader::new(payload);
    let id = reader.str()?.to_owned();
    let name = reader.str()?.to_owned();
    let kind = reader.str()?.to_owned();
    let formula = reader.opt_str()?.map(str::to_owned);
    let parent = reader.opt_str()?.map(str::to_owned);
    let one_letter_code = if reader.bool()? {
        Some(char::from_u32(reader.u32()?).ok_or(corrupt("invalid one-letter code"))?)
    } else {
        None
    };
    let formal_charge = reader.i32()?;
    let release_status = reader.str()?.to_owned();
    let atoms = decode_atoms(&mut reader)?;
    let bonds = decode_bonds(&mut reader, atoms.len())?;
    reader.finish()?;
    Ok(Component {
        id,
        name,
        kind,
        one_letter_code,
        parent,
        formula,
        formal_charge,
        release_status,
        atoms,
        bonds,
    })
}

fn decode_atoms(reader: &mut Reader<'_>) -> Result<Vec<ComponentAtom>, CcdError> {
    let count = reader.count(MIN_ATOM_LEN)?;
    if count > MAX_ATOMS {
        return Err(corrupt("too many atoms"));
    }
    let mut atoms = Vec::with_capacity(count);
    for _ in 0..count {
        let name = reader.str()?.to_owned();
        let element = reader.str()?.to_owned();
        let charge = reader.i8()?;
        let flags = reader.u8()?;
        let stereo = match flags >> ATOM_STEREO_SHIFT {
            0 => AtomStereo::None,
            1 => AtomStereo::R,
            2 => AtomStereo::S,
            _ => return Err(corrupt("invalid atom flags")),
        };
        atoms.push(ComponentAtom {
            name,
            element,
            charge,
            leaving: flags & LEAVING_BIT != 0,
            aromatic: flags & AROMATIC_BIT != 0,
            stereo,
        });
    }
    Ok(atoms)
}

fn decode_bonds(
    reader: &mut Reader<'_>,
    atom_count: usize,
) -> Result<Vec<ComponentBond>, CcdError> {
    let count = reader.count(MIN_BOND_LEN)?;
    let mut bonds = Vec::with_capacity(count);
    for _ in 0..count {
        let atoms = [reader.u16()?, reader.u16()?];
        if atoms.iter().any(|&index| usize::from(index) >= atom_count) {
            return Err(corrupt("bond atom index out of range"));
        }
        let order = *BondOrder::ALL
            .get(usize::from(reader.u8()?))
            .ok_or(corrupt("invalid bond order"))?;
        let flags = reader.u8()?;
        let stereo = match flags >> BOND_STEREO_SHIFT {
            0 => BondStereo::None,
            1 => BondStereo::E,
            2 => BondStereo::Z,
            _ => return Err(corrupt("invalid bond flags")),
        };
        bonds.push(ComponentBond {
            atoms,
            order,
            aromatic: flags & AROMATIC_BIT != 0,
            stereo,
        });
    }
    Ok(bonds)
}

fn atom_flags(atom: &ComponentAtom) -> u8 {
    let stereo = match atom.stereo {
        AtomStereo::None => 0,
        AtomStereo::R => 1,
        AtomStereo::S => 2,
    };
    (u8::from(atom.aromatic) * AROMATIC_BIT)
        | (u8::from(atom.leaving) * LEAVING_BIT)
        | (stereo << ATOM_STEREO_SHIFT)
}

fn bond_flags(bond: &ComponentBond) -> u8 {
    let stereo = match bond.stereo {
        BondStereo::None => 0,
        BondStereo::E => 1,
        BondStereo::Z => 2,
    };
    (u8::from(bond.aromatic) * AROMATIC_BIT) | (stereo << BOND_STEREO_SHIFT)
}

fn order_code(order: BondOrder) -> u8 {
    match order {
        BondOrder::Single => 0,
        BondOrder::Double => 1,
        BondOrder::Triple => 2,
        BondOrder::Quadruple => 3,
        BondOrder::Aromatic => 4,
        BondOrder::Delocalized => 5,
        BondOrder::Pi => 6,
        BondOrder::Polymeric => 7,
    }
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;

    pub(crate) fn sample_component(id: &str) -> Component {
        let atom = |name: &str, element: &str, charge, stereo| ComponentAtom {
            name: name.to_owned(),
            element: element.to_owned(),
            charge,
            leaving: name.starts_with('H'),
            aromatic: element == "N",
            stereo,
        };
        let bond = |a, b, order, stereo| ComponentBond {
            atoms: [a, b],
            order,
            aromatic: order == BondOrder::Aromatic,
            stereo,
        };
        Component {
            id: id.to_owned(),
            name: format!("component {id}"),
            kind: "L-PEPTIDE LINKING".to_owned(),
            one_letter_code: Some('\u{e9}'),
            parent: Some("ALA".to_owned()),
            formula: None,
            formal_charge: -2,
            release_status: "REL".to_owned(),
            atoms: vec![
                atom("CA", "C", 0, AtomStereo::S),
                atom("N1", "N", -1, AtomStereo::R),
                atom("HXT", "H", 1, AtomStereo::None),
                atom("FE", "FE", 3, AtomStereo::None),
            ],
            bonds: vec![
                bond(0, 1, BondOrder::Aromatic, BondStereo::E),
                bond(0, 2, BondOrder::Single, BondStereo::None),
                bond(1, 3, BondOrder::Polymeric, BondStereo::Z),
                bond(3, 3, BondOrder::Pi, BondStereo::None),
            ],
        }
    }

    fn encoded(component: &Component) -> Vec<u8> {
        let mut out = Vec::new();
        encode_record(component, &mut out);
        out
    }

    #[test]
    fn records_round_trip() {
        let component = sample_component("HEM");
        assert_eq!(decode_record(&encoded(&component)).unwrap(), component);
    }

    #[test]
    fn every_bond_order_round_trips() {
        let mut component = sample_component("X");
        for order in BondOrder::ALL {
            component.bonds[0].order = order;
            assert_eq!(decode_record(&encoded(&component)).unwrap(), component);
        }
    }

    #[test]
    fn empty_components_round_trip() {
        let mut component = sample_component("ZN");
        component.atoms.clear();
        component.bonds.clear();
        component.one_letter_code = None;
        assert_eq!(decode_record(&encoded(&component)).unwrap(), component);
    }

    #[test]
    fn truncated_records_are_rejected() {
        let bytes = encoded(&sample_component("ATP"));
        for len in 0..bytes.len() {
            assert!(decode_record(&bytes[..len]).is_err(), "length {len}");
        }
    }

    #[test]
    fn trailing_bytes_are_rejected() {
        let mut bytes = encoded(&sample_component("ATP"));
        bytes.push(0);
        assert!(decode_record(&bytes).is_err());
    }

    #[test]
    fn invalid_bond_fields_are_rejected() {
        let bytes = encoded(&sample_component("ATP"));
        let last_bond = bytes.len() - MIN_BOND_LEN;
        let corruptions = [(last_bond, 9), (last_bond + 4, 8), (last_bond + 5, 6)];
        for (offset, value) in corruptions {
            let mut corrupted = bytes.clone();
            corrupted[offset] = value;
            assert!(
                matches!(
                    decode_record(&corrupted),
                    Err(CcdError::CorruptIndex { .. })
                ),
                "offset {offset}"
            );
        }
    }

    #[test]
    fn flipped_bytes_never_panic() {
        let bytes = encoded(&sample_component("NAG"));
        for offset in 0..bytes.len() {
            for mask in [0x01, 0x80, 0xFF] {
                let mut corrupted = bytes.clone();
                corrupted[offset] ^= mask;
                let outcome = decode_record(&corrupted);
                assert!(
                    matches!(outcome, Ok(_) | Err(CcdError::CorruptIndex { .. })),
                    "offset {offset}"
                );
            }
        }
    }
}
