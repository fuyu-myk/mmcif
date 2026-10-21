use std::collections::HashMap;

use mmcif_cif::{Block, Category};

use crate::component::{
    AtomStereo, BondOrder, BondStereo, Component, ComponentAtom, ComponentBond, MAX_ATOMS,
};
use crate::error::CcdError;

struct Table<'a> {
    block: &'a str,
    name: &'static str,
    category: Option<Category<'a>>,
}

impl<'a> Table<'a> {
    fn new(block: Block<'a>, name: &'static str) -> Self {
        Table {
            block: block.name(),
            name,
            category: block.category(name),
        }
    }

    fn row_count(&self) -> usize {
        self.category.as_ref().map_or(0, Category::row_count)
    }

    fn text(&self, item: &str, row: usize) -> Option<&'a str> {
        self.category.as_ref()?.column(item)?.get(row)?.as_str()
    }

    fn tag(&self, item: &str) -> String {
        format!("_{}.{}", self.name, item)
    }

    fn required(&self, item: &str, row: usize) -> Result<&'a str, CcdError> {
        self.text(item, row).ok_or_else(|| CcdError::MissingItem {
            block: self.block.to_owned(),
            item: self.tag(item),
        })
    }

    fn parsed<T>(
        &self,
        item: &str,
        row: usize,
        parse: impl FnOnce(&str) -> Option<T>,
    ) -> Result<Option<T>, CcdError> {
        let Some(text) = self.text(item, row) else {
            return Ok(None);
        };
        parse(text).map(Some).ok_or_else(|| CcdError::InvalidValue {
            block: self.block.to_owned(),
            item: self.tag(item),
            value: text.to_owned(),
        })
    }

    fn flag(&self, item: &str, row: usize) -> Result<bool, CcdError> {
        Ok(self.parsed(item, row, parse_flag)?.unwrap_or(false))
    }
}

fn parse_flag(text: &str) -> Option<bool> {
    match text {
        "Y" | "y" => Some(true),
        "N" | "n" => Some(false),
        _ => None,
    }
}

impl Component {
    pub fn from_block(block: Block<'_>) -> Result<Component, CcdError> {
        let chem_comp = Table::new(block, "chem_comp");
        let id = chem_comp.required("id", 0)?.to_owned();
        let atoms = read_atoms(&Table::new(block, "chem_comp_atom"))?;
        let bonds = read_bonds(&Table::new(block, "chem_comp_bond"), &atoms)?;
        let owned = |item| chem_comp.text(item, 0).map(str::to_owned);
        Ok(Component {
            id,
            name: owned("name").unwrap_or_default(),
            kind: owned("type").unwrap_or_default(),
            one_letter_code: owned("one_letter_code"),
            parent: owned("mon_nstd_parent_comp_id"),
            formula: owned("formula"),
            formal_charge: chem_comp
                .parsed("pdbx_formal_charge", 0, |text| text.parse().ok())?
                .unwrap_or(0),
            release_status: owned("pdbx_release_status").unwrap_or_default(),
            atoms,
            bonds,
        })
    }
}

fn read_atoms(table: &Table<'_>) -> Result<Vec<ComponentAtom>, CcdError> {
    let count = table.row_count();
    if count > MAX_ATOMS {
        return Err(CcdError::TooManyAtoms {
            block: table.block.to_owned(),
            count,
        });
    }
    (0..count).map(|row| read_atom(table, row)).collect()
}

fn read_atom(table: &Table<'_>, row: usize) -> Result<ComponentAtom, CcdError> {
    Ok(ComponentAtom {
        name: table.required("atom_id", row)?.to_owned(),
        element: table.required("type_symbol", row)?.to_owned(),
        charge: table
            .parsed("charge", row, |text| text.parse().ok())?
            .unwrap_or(0),
        leaving: table.flag("pdbx_leaving_atom_flag", row)?,
        aromatic: table.flag("pdbx_aromatic_flag", row)?,
        stereo: table
            .parsed("pdbx_stereo_config", row, AtomStereo::from_code)?
            .unwrap_or_default(),
    })
}

fn read_bonds(table: &Table<'_>, atoms: &[ComponentAtom]) -> Result<Vec<ComponentBond>, CcdError> {
    let mut indices = HashMap::with_capacity(atoms.len());
    for (atom, index) in atoms.iter().zip(0..=u16::MAX) {
        indices.entry(atom.name.as_str()).or_insert(index);
    }
    (0..table.row_count())
        .map(|row| read_bond(table, row, &indices))
        .collect()
}

fn read_bond(
    table: &Table<'_>,
    row: usize,
    indices: &HashMap<&str, u16>,
) -> Result<ComponentBond, CcdError> {
    let atom = |item| {
        let name = table.required(item, row)?;
        indices
            .get(name)
            .copied()
            .ok_or_else(|| CcdError::UnknownBondAtom {
                block: table.block.to_owned(),
                atom: name.to_owned(),
            })
    };
    let atoms = [atom("atom_id_1")?, atom("atom_id_2")?];
    let code = table.required("value_order", row)?;
    let order = BondOrder::from_code(code).ok_or_else(|| CcdError::UnknownBondOrder {
        block: table.block.to_owned(),
        order: code.to_owned(),
    })?;
    Ok(ComponentBond {
        atoms,
        order,
        aromatic: table.flag("pdbx_aromatic_flag", row)?,
        stereo: table
            .parsed("pdbx_stereo_config", row, BondStereo::from_code)?
            .unwrap_or_default(),
    })
}

#[cfg(test)]
mod tests {
    use mmcif_cif::Document;

    use super::*;

    const SAMPLE: &str = "\
data_TST
_chem_comp.id TST
_chem_comp.name 'test compound'
_chem_comp.type 'L-PEPTIDE LINKING'
_chem_comp.formula 'C2 H4'
_chem_comp.mon_nstd_parent_comp_id ALA
_chem_comp.pdbx_formal_charge -1
_chem_comp.one_letter_code X
_chem_comp.pdbx_release_status REL
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.charge
_chem_comp_atom.pdbx_aromatic_flag
_chem_comp_atom.pdbx_leaving_atom_flag
_chem_comp_atom.pdbx_stereo_config
TST C1 C 0 Y N R
TST \"O1'\" O -1 N Y S
TST FE FE ? ? ? ?
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.value_order
_chem_comp_bond.pdbx_aromatic_flag
_chem_comp_bond.pdbx_stereo_config
TST C1 \"O1'\" DOUB Y E
TST C1 FE SING N N
";

    fn component(text: &str) -> Result<Component, CcdError> {
        let document = Document::parse(text).unwrap();
        Component::from_block(document.block(0).unwrap())
    }

    #[test]
    fn reads_component_header() {
        let tst = component(SAMPLE).unwrap();
        assert_eq!(tst.id(), "TST");
        assert_eq!(tst.name(), "test compound");
        assert_eq!(tst.kind(), "L-PEPTIDE LINKING");
        assert_eq!(tst.formula(), Some("C2 H4"));
        assert_eq!(tst.parent(), Some("ALA"));
        assert_eq!(tst.formal_charge(), -1);
        assert_eq!(tst.one_letter_code(), Some("X"));
        assert_eq!(tst.release_status(), "REL");
    }

    #[test]
    fn reads_atoms_and_bonds() {
        let tst = component(SAMPLE).unwrap();
        let names: Vec<_> = tst.atoms().iter().map(ComponentAtom::name).collect();
        assert_eq!(names, ["C1", "O1'", "FE"]);
        let oxygen = &tst.atoms()[1];
        assert_eq!(oxygen.element(), "O");
        assert_eq!(oxygen.charge(), -1);
        assert!(oxygen.is_leaving() && !oxygen.is_aromatic());
        assert_eq!(oxygen.stereo(), AtomStereo::S);
        let iron = &tst.atoms()[2];
        assert_eq!((iron.charge(), iron.is_leaving()), (0, false));
        assert_eq!(iron.stereo(), AtomStereo::None);
        let bonds = tst.bonds();
        assert_eq!(bonds[0].atom_indices(), [0, 1]);
        assert_eq!(bonds[0].order(), BondOrder::Double);
        assert!(bonds[0].is_aromatic());
        assert_eq!(bonds[0].stereo(), BondStereo::E);
        assert_eq!(bonds[1].atom_indices(), [0, 2]);
        assert_eq!(bonds[1].order(), BondOrder::Single);
    }

    #[test]
    fn missing_optional_items_use_defaults() {
        let bare = component("data_X\n_chem_comp.id X\n_chem_comp.formula ?\n").unwrap();
        assert_eq!(bare.id(), "X");
        assert_eq!(
            (bare.name(), bare.kind(), bare.release_status()),
            ("", "", "")
        );
        assert_eq!((bare.formula(), bare.parent()), (None, None));
        assert_eq!((bare.one_letter_code(), bare.formal_charge()), (None, 0));
        assert!(bare.atoms().is_empty() && bare.bonds().is_empty());
    }

    #[test]
    fn requires_component_id() {
        let error = component("data_X\n_chem_comp.name x\n").unwrap_err();
        assert!(matches!(
            error,
            CcdError::MissingItem { ref block, ref item } if block == "X" && item == "_chem_comp.id"
        ));
        let error = component("data_X\n_other.id X\n").unwrap_err();
        assert!(matches!(error, CcdError::MissingItem { .. }));
    }

    #[test]
    fn requires_atom_names_and_elements() {
        let error = component(&SAMPLE.replace("TST FE FE", "TST FE ?")).unwrap_err();
        assert!(matches!(
            error,
            CcdError::MissingItem { ref item, .. } if item == "_chem_comp_atom.type_symbol"
        ));
    }

    #[test]
    fn rejects_unknown_bond_order() {
        let error = component(&SAMPLE.replace("DOUB", "DOUBLE")).unwrap_err();
        assert!(matches!(
            error,
            CcdError::UnknownBondOrder { ref order, .. } if order == "DOUBLE"
        ));
    }

    #[test]
    fn rejects_bonds_to_unknown_atoms() {
        let error = component(&SAMPLE.replace("TST C1 FE", "TST C1 ZN")).unwrap_err();
        assert!(matches!(
            error,
            CcdError::UnknownBondAtom { ref block, ref atom } if block == "TST" && atom == "ZN"
        ));
    }

    #[test]
    fn rejects_malformed_values() {
        let cases = [
            (
                "C 0 Y N R",
                "C 0 maybe N R",
                "_chem_comp_atom.pdbx_aromatic_flag",
            ),
            (
                "C 0 Y N R",
                "C 0 Y N Q",
                "_chem_comp_atom.pdbx_stereo_config",
            ),
            ("C 0 Y N R", "C 300 Y N R", "_chem_comp_atom.charge"),
            ("DOUB Y E", "DOUB Y R", "_chem_comp_bond.pdbx_stereo_config"),
            ("charge -1", "charge one", "_chem_comp.pdbx_formal_charge"),
        ];
        for (from, to, expected) in cases {
            let error = component(&SAMPLE.replace(from, to)).unwrap_err();
            assert!(
                matches!(error, CcdError::InvalidValue { ref item, .. } if item == expected),
                "{to}: {error}"
            );
        }
    }
}
