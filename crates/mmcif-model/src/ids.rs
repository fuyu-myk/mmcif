use std::borrow::Borrow;
use std::fmt;

macro_rules! string_id {
    ($(#[$meta:meta])* $name:ident) => {
        $(#[$meta])*
        #[derive(Debug, Clone, PartialEq, Eq, Hash, PartialOrd, Ord)]
        pub struct $name(Box<str>);

        impl $name {
            /// Creates an id from `value`.
            pub fn new(value: impl Into<Box<str>>) -> Self {
                $name(value.into())
            }

            /// Returns the id as a string slice.
            pub fn as_str(&self) -> &str {
                &self.0
            }
        }

        impl fmt::Display for $name {
            fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                f.write_str(&self.0)
            }
        }

        impl From<&str> for $name {
            fn from(value: &str) -> Self {
                $name::new(value)
            }
        }

        impl Borrow<str> for $name {
            fn borrow(&self) -> &str {
                &self.0
            }
        }

        impl PartialEq<str> for $name {
            fn eq(&self, other: &str) -> bool {
                &*self.0 == other
            }
        }

        impl PartialEq<&str> for $name {
            fn eq(&self, other: &&str) -> bool {
                &*self.0 == *other
            }
        }

        impl PartialEq<$name> for str {
            fn eq(&self, other: &$name) -> bool {
                self == &*other.0
            }
        }

        impl PartialEq<$name> for &str {
            fn eq(&self, other: &$name) -> bool {
                *self == &*other.0
            }
        }
    };
}

string_id!(
    /// A label chain id, from `_atom_site.label_asym_id`.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_model::LabelAsymId;
    ///
    /// let asym = LabelAsymId::new("B");
    /// assert_eq!(asym.as_str(), "B");
    /// assert_eq!(asym, LabelAsymId::from("B"));
    /// ```
    LabelAsymId
);
string_id!(
    /// An author chain id, from `_atom_site.auth_asym_id` or the chain
    /// column of a PDB record.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_model::AuthAsymId;
    ///
    /// let chain = AuthAsymId::new("H");
    /// assert_eq!(chain.to_string(), "H");
    /// ```
    AuthAsymId
);
string_id!(
    /// An entity id, from `_entity.id`.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_model::EntityId;
    ///
    /// let entity = EntityId::new("1");
    /// assert_eq!(entity.as_str(), "1");
    /// ```
    EntityId
);
string_id!(
    /// A chemical component id, such as `ALA` or `HOH`, from
    /// `_atom_site.label_comp_id`.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_model::CompId;
    ///
    /// let comp = CompId::new("HOH");
    /// assert_eq!(comp.as_str(), "HOH");
    /// ```
    CompId
);
string_id!(
    /// An atom name within a residue, such as `CA`, from
    /// `_atom_site.label_atom_id`.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_model::AtomName;
    ///
    /// let name = AtomName::new("CA");
    /// assert_eq!(name.as_str(), "CA");
    /// ```
    AtomName
);
string_id!(
    /// A transformation operator id, from `_pdbx_struct_oper_list.id`.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_model::OperatorId;
    ///
    /// let operator = OperatorId::new("X0");
    /// assert_eq!(operator.as_str(), "X0");
    /// ```
    OperatorId
);
string_id!(
    /// An assembly id, from `_pdbx_struct_assembly.id`.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_model::AssemblyId;
    ///
    /// let assembly = AssemblyId::new("1");
    /// assert_eq!(assembly.as_str(), "1");
    /// ```
    AssemblyId
);

/// A label residue number, from `_atom_site.label_seq_id`.
///
/// Label numbers start at 1 and number the positions of an entity's polymer
/// sequence.
///
/// # Examples
///
/// ```
/// use mmcif_model::LabelSeqId;
///
/// let seq = LabelSeqId::new(105).unwrap();
/// assert_eq!(seq.get(), 105);
/// assert_eq!(LabelSeqId::new(0), None);
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct LabelSeqId(u32);

impl LabelSeqId {
    /// Creates a label residue number, or returns `None` if `value` is 0.
    pub fn new(value: u32) -> Option<Self> {
        (value > 0).then_some(LabelSeqId(value))
    }

    /// Returns the number.
    pub fn get(self) -> u32 {
        self.0
    }
}

impl fmt::Display for LabelSeqId {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.0.fmt(f)
    }
}

/// The author-assigned residue number, from `_atom_site.auth_seq_id`.
///
/// Author numbers may be zero or negative.
///
/// # Examples
///
/// ```
/// use mmcif_model::AuthSeqId;
///
/// let seq = AuthSeqId::new(-3);
/// assert_eq!(seq.get(), -3);
/// assert_eq!(seq.to_string(), "-3");
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct AuthSeqId(i32);

impl AuthSeqId {
    /// Creates an author residue number.
    pub fn new(value: i32) -> Self {
        AuthSeqId(value)
    }

    /// Returns the number.
    pub fn get(self) -> i32 {
        self.0
    }
}

impl fmt::Display for AuthSeqId {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.0.fmt(f)
    }
}

/// An author insertion code, from `_atom_site.pdbx_PDB_ins_code`.
///
/// [`InsCode::NONE`] is the default, sorts before every code and displays as
/// an empty string.
///
/// # Examples
///
/// ```
/// use mmcif_model::InsCode;
///
/// let code = InsCode::new('A');
/// assert_eq!(code.get(), Some('A'));
/// assert_eq!(code.to_string(), "A");
/// assert_eq!(InsCode::NONE.to_string(), "");
/// assert!(InsCode::NONE < code);
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Default)]
pub struct InsCode(Option<char>);

impl InsCode {
    /// No insertion code.
    pub const NONE: InsCode = InsCode(None);

    /// Creates an insertion code.
    pub fn new(code: char) -> Self {
        InsCode(Some(code))
    }

    /// Returns the code, or `None` for [`InsCode::NONE`].
    pub fn get(self) -> Option<char> {
        self.0
    }
}

impl fmt::Display for InsCode {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self.0 {
            Some(code) => write!(f, "{code}"),
            None => Ok(()),
        }
    }
}

/// An alternate location indicator, from `_atom_site.label_alt_id`.
///
/// [`AltLoc::NONE`] marks an atom with no alternate location. It is the
/// default and displays as an empty string.
///
/// # Examples
///
/// ```
/// use mmcif_model::AltLoc;
///
/// let altloc = AltLoc::new('B');
/// assert_eq!(altloc.get(), Some('B'));
/// assert!(!altloc.is_none());
/// assert!(AltLoc::NONE.is_none());
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Default)]
pub struct AltLoc(Option<char>);

impl AltLoc {
    /// No alternate location.
    pub const NONE: AltLoc = AltLoc(None);

    /// Creates an alternate location indicator.
    pub fn new(code: char) -> Self {
        AltLoc(Some(code))
    }

    /// Returns the indicator, or `None` for [`AltLoc::NONE`].
    pub fn get(self) -> Option<char> {
        self.0
    }

    /// Returns `true` for [`AltLoc::NONE`].
    pub fn is_none(self) -> bool {
        self.0.is_none()
    }
}

impl fmt::Display for AltLoc {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self.0 {
            Some(code) => write!(f, "{code}"),
            None => Ok(()),
        }
    }
}

/// A model number, from `_atom_site.pdbx_PDB_model_num` or a PDB `MODEL`
/// record.
///
/// # Examples
///
/// ```
/// use mmcif_model::ModelNumber;
///
/// assert_eq!(ModelNumber::new(2).get(), 2);
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct ModelNumber(i32);

impl ModelNumber {
    /// Creates a model number.
    pub fn new(value: i32) -> Self {
        ModelNumber(value)
    }

    /// Returns the number.
    pub fn get(self) -> i32 {
        self.0
    }
}

impl fmt::Display for ModelNumber {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.0.fmt(f)
    }
}

/// A chemical element symbol of one or two letters.
///
/// The symbol is held with its first letter in upper case and its second in
/// lower case. Reading does not check the letters against the periodic table;
/// [`Element::is_known`] does. Each element, deuterium `D` and tritium `T` has
/// a constant, such as [`Element::FE`], and [`Element::UNKNOWN`] is `X`.
///
/// # Examples
///
/// ```
/// use mmcif_model::Element;
///
/// let iron = Element::parse("FE").unwrap();
/// assert_eq!(iron, Element::FE);
/// assert_eq!(iron.symbol(), "Fe");
/// assert_eq!(iron.to_string(), "Fe");
/// assert!(!iron.is_hydrogen());
///
/// match Element::parse("c").unwrap() {
///     Element::C => {}
///     other => panic!("expected carbon, got {other}"),
/// }
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct Element([u8; 2]);

impl Element {
    /// Parses an element symbol, ignoring case.
    ///
    /// Returns `None` if `symbol` is not one or two ASCII letters.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_model::Element;
    ///
    /// assert_eq!(Element::parse("c").unwrap().symbol(), "C");
    /// assert_eq!(Element::parse("ZN").unwrap().symbol(), "Zn");
    /// assert_eq!(Element::parse("C1"), None);
    /// assert_eq!(Element::parse(""), None);
    /// ```
    pub fn parse(symbol: &str) -> Option<Element> {
        match symbol.as_bytes() {
            [first] if first.is_ascii_alphabetic() => {
                Some(Element([first.to_ascii_uppercase(), 0]))
            }
            [first, second] if first.is_ascii_alphabetic() && second.is_ascii_alphabetic() => {
                Some(Element([
                    first.to_ascii_uppercase(),
                    second.to_ascii_lowercase(),
                ]))
            }
            _ => None,
        }
    }

    /// Returns the symbol, such as `C` or `Fe`.
    pub fn symbol(&self) -> &str {
        let len = if self.0[1] == 0 { 1 } else { 2 };
        std::str::from_utf8(&self.0[..len]).expect("element symbols are ASCII letters")
    }

    /// Returns `true` for hydrogen and its isotopes, `H`, `D` and `T`.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_model::Element;
    ///
    /// assert!(Element::parse("D").unwrap().is_hydrogen());
    /// assert!(!Element::parse("Hg").unwrap().is_hydrogen());
    /// ```
    pub fn is_hydrogen(&self) -> bool {
        matches!(self.symbol(), "H" | "D" | "T")
    }

    /// Returns `true` if the symbol names one of the 118 elements, or
    /// deuterium `D` or tritium `T`.
    ///
    /// # Examples
    ///
    /// ```
    /// use mmcif_model::Element;
    ///
    /// assert!(Element::parse("Og").unwrap().is_known());
    /// assert!(Element::D.is_known());
    /// assert!(!Element::UNKNOWN.is_known());
    /// assert!(!Element::parse("Qq").unwrap().is_known());
    /// ```
    pub fn is_known(&self) -> bool {
        KNOWN_ELEMENTS.contains(self)
    }

    /// The symbol `X`, which files give for an atom of unknown element, as in
    /// a `UNX` residue.
    pub const UNKNOWN: Element = Element::from_symbol("X");

    const fn from_symbol(symbol: &str) -> Element {
        match symbol.as_bytes() {
            [first] => Element([*first, 0]),
            [first, second] => Element([*first, *second]),
            _ => panic!("element symbols have one or two letters"),
        }
    }
}

macro_rules! elements {
    ($($constant:ident $symbol:literal),* $(,)?) => {
        impl Element {
            $(
                #[doc = concat!("The element with symbol `", $symbol, "`.")]
                pub const $constant: Element = Element::from_symbol($symbol);
            )*
        }

        const KNOWN_ELEMENTS: [Element; [$($symbol),*].len()] = [$(Element::$constant),*];
    };
}

elements!(
    H "H", HE "He", LI "Li", BE "Be", B "B", C "C", N "N", O "O", F "F", NE "Ne", NA "Na",
    MG "Mg", AL "Al", SI "Si", P "P", S "S", CL "Cl", AR "Ar", K "K", CA "Ca", SC "Sc", TI "Ti",
    V "V", CR "Cr", MN "Mn", FE "Fe", CO "Co", NI "Ni", CU "Cu", ZN "Zn", GA "Ga", GE "Ge",
    AS "As", SE "Se", BR "Br", KR "Kr", RB "Rb", SR "Sr", Y "Y", ZR "Zr", NB "Nb", MO "Mo",
    TC "Tc", RU "Ru", RH "Rh", PD "Pd", AG "Ag", CD "Cd", IN "In", SN "Sn", SB "Sb", TE "Te",
    I "I", XE "Xe", CS "Cs", BA "Ba", LA "La", CE "Ce", PR "Pr", ND "Nd", PM "Pm", SM "Sm",
    EU "Eu", GD "Gd", TB "Tb", DY "Dy", HO "Ho", ER "Er", TM "Tm", YB "Yb", LU "Lu", HF "Hf",
    TA "Ta", W "W", RE "Re", OS "Os", IR "Ir", PT "Pt", AU "Au", HG "Hg", TL "Tl", PB "Pb",
    BI "Bi", PO "Po", AT "At", RN "Rn", FR "Fr", RA "Ra", AC "Ac", TH "Th", PA "Pa", U "U",
    NP "Np", PU "Pu", AM "Am", CM "Cm", BK "Bk", CF "Cf", ES "Es", FM "Fm", MD "Md", NO "No",
    LR "Lr", RF "Rf", DB "Db", SG "Sg", BH "Bh", HS "Hs", MT "Mt", DS "Ds", RG "Rg", CN "Cn",
    NH "Nh", FL "Fl", MC "Mc", LV "Lv", TS "Ts", OG "Og", D "D", T "T"
);

impl fmt::Display for Element {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.symbol())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn label_seq_ids_start_at_one() {
        assert!(LabelSeqId::new(0).is_none());
        assert_eq!(LabelSeqId::new(7).unwrap().get(), 7);
    }

    #[test]
    fn elements_normalise_case() {
        assert_eq!(Element::parse("FE").unwrap().symbol(), "Fe");
        assert_eq!(Element::parse("c").unwrap().symbol(), "C");
        assert!(Element::parse("D").unwrap().is_hydrogen());
        assert!(Element::parse("").is_none());
        assert!(Element::parse("C1").is_none());
        assert!(Element::parse("ABC").is_none());
    }

    #[test]
    fn element_constants_match_parsed_symbols() {
        assert_eq!(Element::parse("c"), Some(Element::C));
        assert_eq!(Element::parse("FE"), Some(Element::FE));
        assert_eq!(Element::FE.symbol(), "Fe");
        assert_eq!(Element::OG.symbol(), "Og");
        assert_eq!(Element::parse("x"), Some(Element::UNKNOWN));
        assert!(matches!(Element::parse("Zn").unwrap(), Element::ZN));
    }

    #[test]
    fn known_elements_are_the_periodic_table_and_hydrogen_isotopes() {
        assert_eq!(KNOWN_ELEMENTS.len(), 120);
        assert!(KNOWN_ELEMENTS.iter().all(Element::is_known));
        assert!(Element::D.is_known() && Element::T.is_known());
        assert!(!Element::UNKNOWN.is_known());
        assert!(!Element::parse("Qq").unwrap().is_known());
        let mut symbols: Vec<&str> = KNOWN_ELEMENTS.iter().map(Element::symbol).collect();
        symbols.sort_unstable();
        symbols.dedup();
        assert_eq!(symbols.len(), 120);
    }

    #[test]
    fn ids_compare_with_strings() {
        let chain = AuthAsymId::new("A");
        assert!(chain == "A");
        assert!(chain == *"A");
        assert!(chain != "B");
        assert!("A" == chain);
        assert!(*"A" == chain);
        assert!(CompId::new("HOH") == "HOH");
    }

    #[test]
    fn missing_codes_display_empty() {
        assert_eq!(InsCode::NONE.to_string(), "");
        assert_eq!(InsCode::new('A').to_string(), "A");
        assert_eq!(AltLoc::new('B').get(), Some('B'));
    }
}
