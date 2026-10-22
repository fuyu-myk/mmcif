pub(crate) const MAX_ATOMS: usize = u16::MAX as usize + 1;

/// A chemical component from the dictionary, with its atoms and the bonds
/// between them.
///
/// # Examples
///
/// ```
/// use mmcif_ccd::{CcdIndex, build_index};
///
/// # let dictionary = include_str!("../tests/fixtures/components-sample.cif");
/// # let mut bytes = Vec::new();
/// # build_index(dictionary, "2026-10-20", &mut bytes)?;
/// let index = CcdIndex::from_bytes(bytes)?;
/// let alanine = index.get("ALA")?.unwrap();
///
/// assert_eq!(alanine.name(), "ALANINE");
/// assert_eq!(alanine.kind(), "L-PEPTIDE LINKING");
/// assert_eq!(alanine.formula(), Some("C3 H7 N O2"));
/// assert_eq!((alanine.atoms().len(), alanine.bonds().len()), (13, 12));
/// # Ok::<(), mmcif_ccd::CcdError>(())
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Component {
    pub(crate) id: String,
    pub(crate) name: String,
    pub(crate) kind: String,
    pub(crate) one_letter_code: Option<String>,
    pub(crate) parent: Option<String>,
    pub(crate) formula: Option<String>,
    pub(crate) formal_charge: i32,
    pub(crate) release_status: String,
    pub(crate) atoms: Vec<ComponentAtom>,
    pub(crate) bonds: Vec<ComponentBond>,
}

impl Component {
    /// Returns the component id, such as `ALA`.
    pub fn id(&self) -> &str {
        &self.id
    }

    /// Returns the component name, such as `ALANINE`, or an empty string if
    /// the dictionary gives none.
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Returns the component type, such as `L-PEPTIDE LINKING`, or an empty
    /// string if the dictionary gives none.
    pub fn kind(&self) -> &str {
        &self.kind
    }

    /// Returns the one-letter code, such as `A` for alanine, if the
    /// dictionary gives one.
    pub fn one_letter_code(&self) -> Option<&str> {
        self.one_letter_code.as_deref()
    }

    /// Returns the id of the parent component, such as `MET` for
    /// selenomethionine, if the dictionary gives one.
    pub fn parent(&self) -> Option<&str> {
        self.parent.as_deref()
    }

    /// Returns the chemical formula, such as `C3 H7 N O2`, if the dictionary
    /// gives one.
    pub fn formula(&self) -> Option<&str> {
        self.formula.as_deref()
    }

    /// Returns the formal charge of the component, or 0 if the dictionary
    /// gives none.
    pub fn formal_charge(&self) -> i32 {
        self.formal_charge
    }

    /// Returns the release status, such as `REL`, or an empty string if the
    /// dictionary gives none.
    pub fn release_status(&self) -> &str {
        &self.release_status
    }

    /// Returns the atoms in the order the dictionary lists them.
    pub fn atoms(&self) -> &[ComponentAtom] {
        &self.atoms
    }

    /// Returns the bonds in the order the dictionary lists them.
    pub fn bonds(&self) -> &[ComponentBond] {
        &self.bonds
    }
}

/// An atom of a [`Component`].
///
/// # Examples
///
/// ```
/// use mmcif_ccd::{AtomStereo, CcdIndex, build_index};
///
/// # let dictionary = include_str!("../tests/fixtures/components-sample.cif");
/// # let mut bytes = Vec::new();
/// # build_index(dictionary, "2026-10-20", &mut bytes)?;
/// let index = CcdIndex::from_bytes(bytes)?;
/// let alanine = index.get("ALA")?.unwrap();
///
/// let alpha = &alanine.atoms()[1];
/// assert_eq!((alpha.name(), alpha.element()), ("CA", "C"));
/// assert_eq!(alpha.stereo(), AtomStereo::S);
/// assert!(!alpha.is_leaving());
/// # Ok::<(), mmcif_ccd::CcdError>(())
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ComponentAtom {
    pub(crate) name: String,
    pub(crate) element: String,
    pub(crate) charge: i8,
    pub(crate) leaving: bool,
    pub(crate) aromatic: bool,
    pub(crate) stereo: AtomStereo,
}

impl ComponentAtom {
    /// Returns the atom name, such as `CA`.
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Returns the element symbol as the dictionary writes it, such as `C` or
    /// `FE`.
    pub fn element(&self) -> &str {
        &self.element
    }

    /// Returns the formal charge of the atom, or 0 if the dictionary gives
    /// none.
    pub fn charge(&self) -> i8 {
        self.charge
    }

    /// Returns `true` if the dictionary flags the atom as a leaving atom.
    pub fn is_leaving(&self) -> bool {
        self.leaving
    }

    /// Returns `true` if the dictionary flags the atom as aromatic.
    pub fn is_aromatic(&self) -> bool {
        self.aromatic
    }

    /// Returns the stereo configuration of the atom.
    pub fn stereo(&self) -> AtomStereo {
        self.stereo
    }
}

/// A bond between two atoms of a [`Component`].
///
/// # Examples
///
/// ```
/// use mmcif_ccd::{BondOrder, CcdIndex, build_index};
///
/// # let dictionary = include_str!("../tests/fixtures/components-sample.cif");
/// # let mut bytes = Vec::new();
/// # build_index(dictionary, "2026-10-20", &mut bytes)?;
/// let index = CcdIndex::from_bytes(bytes)?;
/// let alanine = index.get("ALA")?.unwrap();
///
/// let carbonyl = alanine.bonds()[6];
/// let names = carbonyl.atom_indices().map(|i| alanine.atoms()[usize::from(i)].name());
/// assert_eq!(names, ["C", "O"]);
/// assert_eq!(carbonyl.order(), BondOrder::Double);
/// # Ok::<(), mmcif_ccd::CcdError>(())
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ComponentBond {
    pub(crate) atoms: [u16; 2],
    pub(crate) order: BondOrder,
    pub(crate) aromatic: bool,
    pub(crate) stereo: BondStereo,
}

impl ComponentBond {
    /// Returns the positions of the two bonded atoms in
    /// [`Component::atoms`], in the order the dictionary lists them.
    pub fn atom_indices(&self) -> [u16; 2] {
        self.atoms
    }

    /// Returns the bond order.
    pub fn order(&self) -> BondOrder {
        self.order
    }

    /// Returns `true` if the dictionary flags the bond as aromatic.
    pub fn is_aromatic(&self) -> bool {
        self.aromatic
    }

    /// Returns the stereo configuration of the bond.
    pub fn stereo(&self) -> BondStereo {
        self.stereo
    }
}

/// The order of a bond, as given by `_chem_comp_bond.value_order`.
///
/// # Examples
///
/// ```
/// use mmcif_ccd::BondOrder;
///
/// assert_eq!(BondOrder::from_code("DOUB"), Some(BondOrder::Double));
/// assert_eq!(BondOrder::Double.code(), "DOUB");
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum BondOrder {
    /// A single bond, code `SING`.
    Single,
    /// A double bond, code `DOUB`.
    Double,
    /// A triple bond, code `TRIP`.
    Triple,
    /// A quadruple bond, code `QUAD`.
    Quadruple,
    /// An aromatic bond, code `AROM`.
    Aromatic,
    /// A delocalised bond, code `DELO`.
    Delocalized,
    /// A pi bond, code `PI`.
    Pi,
    /// A polymeric bond, code `POLY`.
    Polymeric,
}

impl BondOrder {
    /// Every bond order, in declaration order.
    pub const ALL: [BondOrder; 8] = [
        BondOrder::Single,
        BondOrder::Double,
        BondOrder::Triple,
        BondOrder::Quadruple,
        BondOrder::Aromatic,
        BondOrder::Delocalized,
        BondOrder::Pi,
        BondOrder::Polymeric,
    ];

    /// Returns the bond order for a dictionary code such as `SING`, ignoring
    /// ASCII case, or `None` if the code is not one of the codes listed on
    /// the variants.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_ccd::BondOrder;
    ///
    /// assert_eq!(BondOrder::from_code("SING"), Some(BondOrder::Single));
    /// assert_eq!(BondOrder::from_code("arom"), Some(BondOrder::Aromatic));
    /// assert_eq!(BondOrder::from_code("SINGLE"), None);
    /// ```
    pub fn from_code(code: &str) -> Option<BondOrder> {
        BondOrder::ALL
            .into_iter()
            .find(|order| order.code().eq_ignore_ascii_case(code))
    }

    /// Returns the dictionary code for the bond order, in upper case.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_ccd::BondOrder;
    ///
    /// assert_eq!(BondOrder::Triple.code(), "TRIP");
    /// assert_eq!(BondOrder::Pi.code(), "PI");
    /// ```
    pub fn code(self) -> &'static str {
        match self {
            BondOrder::Single => "SING",
            BondOrder::Double => "DOUB",
            BondOrder::Triple => "TRIP",
            BondOrder::Quadruple => "QUAD",
            BondOrder::Aromatic => "AROM",
            BondOrder::Delocalized => "DELO",
            BondOrder::Pi => "PI",
            BondOrder::Polymeric => "POLY",
        }
    }
}

/// The stereo configuration of an atom, as given by
/// `_chem_comp_atom.pdbx_stereo_config`.
///
/// The default is [`AtomStereo::None`].
///
/// # Examples
///
/// ```
/// use mmcif_ccd::AtomStereo;
///
/// assert_eq!(AtomStereo::from_code("S"), Some(AtomStereo::S));
/// assert_eq!(AtomStereo::default(), AtomStereo::None);
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub enum AtomStereo {
    /// No stereo configuration, code `N`.
    #[default]
    None,
    /// The R configuration, code `R`.
    R,
    /// The S configuration, code `S`.
    S,
}

impl AtomStereo {
    /// Returns the configuration for a code of `N`, `R` or `S` in either
    /// case, or `None` for any other code.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_ccd::AtomStereo;
    ///
    /// assert_eq!(AtomStereo::from_code("R"), Some(AtomStereo::R));
    /// assert_eq!(AtomStereo::from_code("n"), Some(AtomStereo::None));
    /// assert_eq!(AtomStereo::from_code("E"), None);
    /// ```
    pub fn from_code(code: &str) -> Option<AtomStereo> {
        match code {
            "N" | "n" => Some(AtomStereo::None),
            "R" | "r" => Some(AtomStereo::R),
            "S" | "s" => Some(AtomStereo::S),
            _ => None,
        }
    }
}

/// The stereo configuration of a bond, as given by
/// `_chem_comp_bond.pdbx_stereo_config`.
///
/// The default is [`BondStereo::None`].
///
/// # Examples
///
/// ```
/// use mmcif_ccd::BondStereo;
///
/// assert_eq!(BondStereo::from_code("Z"), Some(BondStereo::Z));
/// assert_eq!(BondStereo::default(), BondStereo::None);
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub enum BondStereo {
    /// No stereo configuration, code `N`.
    #[default]
    None,
    /// The E configuration, code `E`.
    E,
    /// The Z configuration, code `Z`.
    Z,
}

impl BondStereo {
    /// Returns the configuration for a code of `N`, `E` or `Z` in either
    /// case, or `None` for any other code.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_ccd::BondStereo;
    ///
    /// assert_eq!(BondStereo::from_code("E"), Some(BondStereo::E));
    /// assert_eq!(BondStereo::from_code("z"), Some(BondStereo::Z));
    /// assert_eq!(BondStereo::from_code("R"), None);
    /// ```
    pub fn from_code(code: &str) -> Option<BondStereo> {
        match code {
            "N" | "n" => Some(BondStereo::None),
            "E" | "e" => Some(BondStereo::E),
            "Z" | "z" => Some(BondStereo::Z),
            _ => None,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_every_bond_order_code() {
        let codes = ["SING", "DOUB", "TRIP", "QUAD", "AROM", "DELO", "PI", "POLY"];
        for (code, order) in codes.into_iter().zip(BondOrder::ALL) {
            assert_eq!(BondOrder::from_code(code), Some(order));
            assert_eq!(order.code(), code);
        }
    }

    #[test]
    fn bond_order_codes_ignore_case() {
        assert_eq!(BondOrder::from_code("sing"), Some(BondOrder::Single));
        assert_eq!(BondOrder::from_code("Arom"), Some(BondOrder::Aromatic));
        assert_eq!(BondOrder::from_code("pi"), Some(BondOrder::Pi));
    }

    #[test]
    fn rejects_unknown_bond_orders() {
        for code in ["", "SINGLE", "1", "?", "SIN", "DOUBX"] {
            assert_eq!(BondOrder::from_code(code), None, "{code}");
        }
    }

    #[test]
    fn parses_stereo_codes() {
        assert_eq!(AtomStereo::from_code("R"), Some(AtomStereo::R));
        assert_eq!(AtomStereo::from_code("s"), Some(AtomStereo::S));
        assert_eq!(AtomStereo::from_code("N"), Some(AtomStereo::None));
        assert_eq!(AtomStereo::from_code("E"), None);
        assert_eq!(BondStereo::from_code("E"), Some(BondStereo::E));
        assert_eq!(BondStereo::from_code("z"), Some(BondStereo::Z));
        assert_eq!(BondStereo::from_code("N"), Some(BondStereo::None));
        assert_eq!(BondStereo::from_code("R"), None);
    }
}
