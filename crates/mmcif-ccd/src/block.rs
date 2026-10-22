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

    fn missing(&self, item: &str) -> CcdError {
        CcdError::MissingItem {
            block: self.block.to_owned(),
            item: self.tag(item),
        }
    }

    fn required(&self, item: &str, row: usize) -> Result<&'a str, CcdError> {
        self.text(item, row).ok_or_else(|| self.missing(item))
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
    /// Reads a component from one data block of the Chemical Component
    /// Dictionary.
    ///
    /// The block must define one component: its header is the one row of
    /// `_chem_comp`, its atoms are the rows of `_chem_comp_atom` and its bonds
    /// the rows of `_chem_comp_bond`. An atom or bond row without a `comp_id`
    /// belongs to the component. An item that is missing or holds `?` or `.`
    /// is absent. An absent name, type or release status reads as an empty
    /// string, an absent formula, parent or one-letter code as `None`, an
    /// absent charge as 0, an absent flag as `false`, and an absent stereo
    /// configuration as [`AtomStereo::None`] or [`BondStereo::None`].
    ///
    /// To read a block that defines several components, such as a structure
    /// file with its component definitions, use [`Component::all_from_block`].
    ///
    /// # Errors
    ///
    /// Returns [`CcdError::ManyComponents`] if `_chem_comp` has more than one
    /// row; [`CcdError::MissingItem`] if `_chem_comp.id` is absent, an atom
    /// lacks `atom_id` or `type_symbol`, or a bond lacks `atom_id_1`,
    /// `atom_id_2` or `value_order`; [`CcdError::UnknownComponent`] if an atom
    /// or bond names another component; [`CcdError::TooManyAtoms`] if the
    /// component has more than 65,536 atoms; [`CcdError::InvalidValue`] if the
    /// formal charge is not an integer that fits in an `i32`, an atom charge is
    /// not an integer that fits in an `i8`, an aromatic or leaving atom flag
    /// is not `Y` or `N` in either case, or a stereo configuration is not a
    /// code that [`AtomStereo::from_code`] or [`BondStereo::from_code`]
    /// accepts; [`CcdError::UnknownBondAtom`] if a bond names an atom that is
    /// not one of the component's atoms; and [`CcdError::UnknownBondOrder`] if
    /// a bond order is not a code that [`BondOrder::from_code`] accepts.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_ccd::{CcdError, Component};
    /// use mmcif_cif::Document;
    ///
    /// let document = Document::parse(include_str!("../tests/fixtures/components-sample.cif"))?;
    /// let alanine = Component::from_block(document.block(0).unwrap())?;
    /// assert_eq!(alanine.id(), "ALA");
    /// assert_eq!(alanine.one_letter_code(), Some("A"));
    ///
    /// let unnamed = Document::parse("data_X\n_chem_comp.name x\n")?;
    /// let error = Component::from_block(unnamed.block(0).unwrap()).unwrap_err();
    /// assert!(matches!(error, CcdError::MissingItem { item, .. } if item == "_chem_comp.id"));
    ///
    /// let structure = Document::parse("data_1ABC\nloop_\n_chem_comp.id\nHOH\nZN\n")?;
    /// let error = Component::from_block(structure.block(0).unwrap()).unwrap_err();
    /// assert!(matches!(error, CcdError::ManyComponents { count: 2, .. }));
    /// # Ok::<(), CcdError>(())
    /// ```
    pub fn from_block(block: Block<'_>) -> Result<Component, CcdError> {
        let chem_comp = Table::new(block, "chem_comp");
        let count = chem_comp.row_count();
        if count > 1 {
            return Err(CcdError::ManyComponents {
                block: block.name().to_owned(),
                count,
            });
        }
        let id = chem_comp.required("id", 0)?;
        let components = read_components(block, &chem_comp, &[id])?;
        Ok(components
            .into_iter()
            .next()
            .expect("one component is read for each id"))
    }

    /// Reads every component that one data block defines, in `_chem_comp`
    /// order.
    ///
    /// Each atom and bond row belongs to the component its `comp_id` names; a
    /// row without a `comp_id` belongs to the only component when the block
    /// defines one. A block without `_chem_comp` defines no components. Items
    /// are read as [`Component::from_block`] reads them.
    ///
    /// # Errors
    ///
    /// Returns [`CcdError::DuplicateId`] if two rows of `_chem_comp` have the
    /// same id; [`CcdError::MissingItem`] if an atom or bond row has no
    /// `comp_id` while the block defines more than one component; and
    /// otherwise the errors of [`Component::from_block`] other than
    /// [`CcdError::ManyComponents`].
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_ccd::{CcdError, Component};
    /// use mmcif_cif::Document;
    ///
    /// let text = "\
    /// data_1ABC
    /// loop_
    /// _chem_comp.id
    /// HOH
    /// ZN
    /// loop_
    /// _chem_comp_atom.comp_id
    /// _chem_comp_atom.atom_id
    /// _chem_comp_atom.type_symbol
    /// HOH O O
    /// ZN ZN ZN
    /// ";
    /// let document = Document::parse(text)?;
    /// let components = Component::all_from_block(document.block(0).unwrap())?;
    ///
    /// let ids: Vec<&str> = components.iter().map(Component::id).collect();
    /// assert_eq!(ids, ["HOH", "ZN"]);
    /// assert_eq!(components[1].atoms()[0].name(), "ZN");
    /// # Ok::<(), CcdError>(())
    /// ```
    pub fn all_from_block(block: Block<'_>) -> Result<Vec<Component>, CcdError> {
        let chem_comp = Table::new(block, "chem_comp");
        let ids = (0..chem_comp.row_count())
            .map(|row| chem_comp.required("id", row))
            .collect::<Result<Vec<_>, _>>()?;
        read_components(block, &chem_comp, &ids)
    }
}

fn read_components(
    block: Block<'_>,
    chem_comp: &Table<'_>,
    ids: &[&str],
) -> Result<Vec<Component>, CcdError> {
    let positions = component_positions(ids)?;
    let atoms = Table::new(block, "chem_comp_atom");
    let bonds = Table::new(block, "chem_comp_bond");
    let atom_rows = rows_by_component(&atoms, &positions)?;
    let bond_rows = rows_by_component(&bonds, &positions)?;
    ids.iter()
        .zip(atom_rows.iter().zip(&bond_rows))
        .enumerate()
        .map(|(row, (id, (atom_rows, bond_rows)))| {
            let component_atoms = read_atoms(&atoms, atom_rows)?;
            let component_bonds = read_bonds(&bonds, bond_rows, &component_atoms)?;
            read_header(chem_comp, row, id, component_atoms, component_bonds)
        })
        .collect()
}

fn component_positions<'a>(ids: &[&'a str]) -> Result<HashMap<&'a str, usize>, CcdError> {
    let mut positions = HashMap::with_capacity(ids.len());
    for (position, id) in ids.iter().enumerate() {
        if positions.insert(*id, position).is_some() {
            return Err(CcdError::DuplicateId {
                id: (*id).to_owned(),
            });
        }
    }
    Ok(positions)
}

fn rows_by_component(
    table: &Table<'_>,
    positions: &HashMap<&str, usize>,
) -> Result<Vec<Vec<usize>>, CcdError> {
    let mut rows = vec![Vec::new(); positions.len()];
    for row in 0..table.row_count() {
        let position = match table.text("comp_id", row) {
            Some(component) => {
                *positions
                    .get(component)
                    .ok_or_else(|| CcdError::UnknownComponent {
                        block: table.block.to_owned(),
                        item: table.tag("comp_id"),
                        component: component.to_owned(),
                    })?
            }
            None if positions.len() == 1 => 0,
            None => return Err(table.missing("comp_id")),
        };
        rows[position].push(row);
    }
    Ok(rows)
}

fn read_header(
    chem_comp: &Table<'_>,
    row: usize,
    id: &str,
    atoms: Vec<ComponentAtom>,
    bonds: Vec<ComponentBond>,
) -> Result<Component, CcdError> {
    let owned = |item| chem_comp.text(item, row).map(str::to_owned);
    Ok(Component {
        id: id.to_owned(),
        name: owned("name").unwrap_or_default(),
        kind: owned("type").unwrap_or_default(),
        one_letter_code: owned("one_letter_code"),
        parent: owned("mon_nstd_parent_comp_id"),
        formula: owned("formula"),
        formal_charge: chem_comp
            .parsed("pdbx_formal_charge", row, |text| text.parse().ok())?
            .unwrap_or(0),
        release_status: owned("pdbx_release_status").unwrap_or_default(),
        atoms,
        bonds,
    })
}

fn read_atoms(table: &Table<'_>, rows: &[usize]) -> Result<Vec<ComponentAtom>, CcdError> {
    if rows.len() > MAX_ATOMS {
        return Err(CcdError::TooManyAtoms {
            block: table.block.to_owned(),
            count: rows.len(),
        });
    }
    rows.iter().map(|&row| read_atom(table, row)).collect()
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

fn read_bonds(
    table: &Table<'_>,
    rows: &[usize],
    atoms: &[ComponentAtom],
) -> Result<Vec<ComponentBond>, CcdError> {
    let mut indices = HashMap::with_capacity(atoms.len());
    for (atom, index) in atoms.iter().zip(0..=u16::MAX) {
        indices.entry(atom.name.as_str()).or_insert(index);
    }
    rows.iter()
        .map(|&row| read_bond(table, row, &indices))
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

    const TWO_COMPONENTS: &str = "\
data_1ABC
loop_
_chem_comp.id
_chem_comp.name
TST 'test compound'
ZN 'zinc ion'
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
TST C1 C
ZN ZN ZN
TST O1 O
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.value_order
TST C1 O1 DOUB
";

    fn all_components(text: &str) -> Result<Vec<Component>, CcdError> {
        let document = Document::parse(text).unwrap();
        Component::all_from_block(document.block(0).unwrap())
    }

    #[test]
    fn rejects_a_block_with_many_components() {
        let error = component(TWO_COMPONENTS).unwrap_err();
        assert!(matches!(
            error,
            CcdError::ManyComponents { ref block, count: 2 } if block == "1ABC"
        ));
    }

    #[test]
    fn rejects_atoms_and_bonds_of_another_component() {
        let error = component(&SAMPLE.replace("TST FE FE", "OTH FE FE")).unwrap_err();
        assert!(matches!(
            error,
            CcdError::UnknownComponent { ref item, ref component, .. }
                if item == "_chem_comp_atom.comp_id" && component == "OTH"
        ));
        let error = component(&SAMPLE.replace("TST C1 FE", "OTH C1 FE")).unwrap_err();
        assert!(matches!(
            error,
            CcdError::UnknownComponent { ref item, ref component, .. }
                if item == "_chem_comp_bond.comp_id" && component == "OTH"
        ));
    }

    #[test]
    fn atoms_without_a_component_id_belong_to_the_only_component() {
        let tst = component(&SAMPLE.replace("TST FE FE", "? FE FE")).unwrap();
        assert_eq!(tst.atoms().len(), 3);
    }

    #[test]
    fn reads_every_component_with_its_own_atoms_and_bonds() {
        let components = all_components(TWO_COMPONENTS).unwrap();
        let ids: Vec<_> = components.iter().map(Component::id).collect();
        assert_eq!(ids, ["TST", "ZN"]);
        let names: Vec<_> = components[0]
            .atoms()
            .iter()
            .map(ComponentAtom::name)
            .collect();
        assert_eq!(names, ["C1", "O1"]);
        assert_eq!(components[0].bonds()[0].atom_indices(), [0, 1]);
        assert_eq!(components[1].name(), "zinc ion");
        assert_eq!(components[1].atoms()[0].element(), "ZN");
        assert!(components[1].bonds().is_empty());
    }

    #[test]
    fn reading_every_component_agrees_with_reading_one() {
        let one = component(SAMPLE).unwrap();
        assert_eq!(all_components(SAMPLE).unwrap(), [one]);
        assert!(
            all_components("data_X\n_cell.length_a 1\n")
                .unwrap()
                .is_empty()
        );
    }

    #[test]
    fn reading_every_component_rejects_ambiguous_rows() {
        let error = all_components(&TWO_COMPONENTS.replace("ZN ZN ZN", "OTH ZN ZN")).unwrap_err();
        assert!(matches!(
            error,
            CcdError::UnknownComponent { ref component, .. } if component == "OTH"
        ));
        let error = all_components(&TWO_COMPONENTS.replace("ZN ZN ZN", "? ZN ZN")).unwrap_err();
        assert!(matches!(
            error,
            CcdError::MissingItem { ref item, .. } if item == "_chem_comp_atom.comp_id"
        ));
        let error =
            all_components(&TWO_COMPONENTS.replace("ZN 'zinc ion'", "TST 'again'")).unwrap_err();
        assert!(matches!(error, CcdError::DuplicateId { ref id } if id == "TST"));
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
