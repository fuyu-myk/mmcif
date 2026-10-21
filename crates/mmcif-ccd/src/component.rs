pub(crate) const MAX_ATOMS: usize = u16::MAX as usize + 1;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Component {
    pub(crate) id: String,
    pub(crate) name: String,
    pub(crate) kind: String,
    pub(crate) one_letter_code: Option<char>,
    pub(crate) parent: Option<String>,
    pub(crate) formula: Option<String>,
    pub(crate) formal_charge: i32,
    pub(crate) release_status: String,
    pub(crate) atoms: Vec<ComponentAtom>,
    pub(crate) bonds: Vec<ComponentBond>,
}

impl Component {
    pub fn id(&self) -> &str {
        &self.id
    }

    pub fn name(&self) -> &str {
        &self.name
    }

    pub fn kind(&self) -> &str {
        &self.kind
    }

    pub fn one_letter_code(&self) -> Option<char> {
        self.one_letter_code
    }

    pub fn parent(&self) -> Option<&str> {
        self.parent.as_deref()
    }

    pub fn formula(&self) -> Option<&str> {
        self.formula.as_deref()
    }

    pub fn formal_charge(&self) -> i32 {
        self.formal_charge
    }

    pub fn release_status(&self) -> &str {
        &self.release_status
    }

    pub fn atoms(&self) -> &[ComponentAtom] {
        &self.atoms
    }

    pub fn bonds(&self) -> &[ComponentBond] {
        &self.bonds
    }
}

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
    pub fn name(&self) -> &str {
        &self.name
    }

    pub fn element(&self) -> &str {
        &self.element
    }

    pub fn charge(&self) -> i8 {
        self.charge
    }

    pub fn is_leaving(&self) -> bool {
        self.leaving
    }

    pub fn is_aromatic(&self) -> bool {
        self.aromatic
    }

    pub fn stereo(&self) -> AtomStereo {
        self.stereo
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ComponentBond {
    pub(crate) atoms: [u16; 2],
    pub(crate) order: BondOrder,
    pub(crate) aromatic: bool,
    pub(crate) stereo: BondStereo,
}

impl ComponentBond {
    pub fn atom_indices(&self) -> [u16; 2] {
        self.atoms
    }

    pub fn order(&self) -> BondOrder {
        self.order
    }

    pub fn is_aromatic(&self) -> bool {
        self.aromatic
    }

    pub fn stereo(&self) -> BondStereo {
        self.stereo
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum BondOrder {
    Single,
    Double,
    Triple,
    Quadruple,
    Aromatic,
    Delocalized,
    Pi,
    Polymeric,
}

impl BondOrder {
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

    pub fn from_code(code: &str) -> Option<BondOrder> {
        BondOrder::ALL
            .into_iter()
            .find(|order| order.code().eq_ignore_ascii_case(code))
    }

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

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub enum AtomStereo {
    #[default]
    None,
    R,
    S,
}

impl AtomStereo {
    pub fn from_code(code: &str) -> Option<AtomStereo> {
        match code {
            "N" | "n" => Some(AtomStereo::None),
            "R" | "r" => Some(AtomStereo::R),
            "S" | "s" => Some(AtomStereo::S),
            _ => None,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub enum BondStereo {
    #[default]
    None,
    E,
    Z,
}

impl BondStereo {
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
