use std::collections::HashMap;
use std::path::Path;

use crate::entity::{Entity, EntityKind};
use crate::error::ReadError;
use crate::ids::{
    AltLoc, AtomName, AuthAsymId, AuthSeqId, CompId, Element, EntityId, InsCode, LabelAsymId,
    ModelNumber,
};
use crate::numbering::Numbering;
use crate::structure::{Atom, Chain, Model, Residue, Structure};

const WATER_NAMES: [&str; 4] = ["HOH", "WAT", "DOD", "H2O"];

struct Line<'a> {
    number: usize,
    text: &'a str,
}

impl<'a> Line<'a> {
    fn error(&self, problem: &'static str) -> ReadError {
        ReadError::Pdb {
            line: self.number,
            problem,
        }
    }

    fn columns(&self, start: usize, end: usize) -> Result<&'a str, ReadError> {
        let len = self.text.len();
        self.text
            .get(start.min(len)..end.min(len))
            .ok_or_else(|| self.error("non-ASCII text in fixed columns"))
    }

    fn field(&self, start: usize, end: usize) -> Result<&'a str, ReadError> {
        Ok(self.columns(start, end)?.trim())
    }

    fn code(&self, column: usize) -> Result<Option<char>, ReadError> {
        Ok(self
            .columns(column, column + 1)?
            .chars()
            .next()
            .filter(|c| *c != ' '))
    }

    fn number<T: std::str::FromStr>(
        &self,
        start: usize,
        end: usize,
        problem: &'static str,
    ) -> Result<T, ReadError> {
        self.field(start, end)?
            .parse()
            .map_err(|_| self.error(problem))
    }

    fn optional<T>(
        &self,
        start: usize,
        end: usize,
        problem: &'static str,
        decode: impl FnOnce(&str) -> Option<T>,
    ) -> Result<Option<T>, ReadError> {
        let raw = self.columns(start, end)?;
        if raw.trim().is_empty() {
            return Ok(None);
        }
        decode(raw).map(Some).ok_or_else(|| self.error(problem))
    }

    fn number_or<T: std::str::FromStr>(
        &self,
        start: usize,
        end: usize,
        default: T,
        problem: &'static str,
    ) -> Result<T, ReadError> {
        match self.field(start, end)? {
            "" => Ok(default),
            text => text.parse().map_err(|_| self.error(problem)),
        }
    }
}

pub(crate) fn decode_hybrid36(field: &str, width: u32) -> Option<i64> {
    if let Ok(value) = field.trim().parse::<i64>() {
        return Some(value);
    }
    if field.len() != width as usize {
        return None;
    }
    let first = field.chars().next()?;
    let digits_valid = |upper: bool| {
        field.chars().all(|c| {
            c.is_ascii_digit()
                || if upper {
                    c.is_ascii_uppercase()
                } else {
                    c.is_ascii_lowercase()
                }
        })
    };
    let raw = i64::from_str_radix(field, 36).ok()?;
    let block = 10 * 36_i64.pow(width - 1);
    let decimal_limit = 10_i64.pow(width);
    if first.is_ascii_uppercase() && digits_valid(true) {
        Some(raw - block + decimal_limit)
    } else if first.is_ascii_lowercase() && digits_valid(false) {
        Some(raw - block + decimal_limit + 26 * 36_i64.pow(width - 1))
    } else {
        None
    }
}

fn parse_charge(field: &str) -> Option<i8> {
    let bytes = field.as_bytes();
    match bytes {
        [digit @ b'0'..=b'9', sign @ (b'+' | b'-')] => {
            let magnitude = (digit - b'0') as i8;
            Some(if *sign == b'-' { -magnitude } else { magnitude })
        }
        _ => None,
    }
}

fn infer_element(name: &str, hetero: bool) -> Option<Element> {
    let first = name.chars().next()?;
    if first == ' ' || first.is_ascii_digit() {
        return Element::parse(name.get(1..2)?);
    }
    if hetero {
        let letters: String = name
            .chars()
            .take(2)
            .take_while(char::is_ascii_alphabetic)
            .collect();
        return Element::parse(&letters);
    }
    Element::parse(&first.to_string())
}

struct PdbResidue {
    chain: AuthAsymId,
    segment: u32,
    residue: Residue,
}

struct PdbModel {
    number: i32,
    residues: Vec<PdbResidue>,
    terminations: HashMap<AuthAsymId, u32>,
}

impl PdbModel {
    fn new(number: i32) -> Self {
        PdbModel {
            number,
            residues: Vec::new(),
            terminations: HashMap::new(),
        }
    }

    fn add_atom(&mut self, line: &Line<'_>) -> Result<(), ReadError> {
        let hetero = line.text.starts_with("HETATM");
        let name_columns = line.columns(12, 16)?;
        let chain = AuthAsymId::new(line.code(21)?.map_or(String::new(), String::from));
        let auth_seq = decode_hybrid36(line.columns(22, 26)?, 4)
            .and_then(|value| i32::try_from(value).ok())
            .ok_or_else(|| line.error("invalid residue number"))?;
        let ins_code = line.code(26)?.map_or(InsCode::NONE, InsCode::new);
        let comp = line.field(17, 21)?;
        let element_field = line.field(76, 78)?;
        let element = if element_field.is_empty() {
            infer_element(name_columns, hetero)
        } else {
            Element::parse(element_field)
        }
        .ok_or_else(|| line.error("missing or invalid element"))?;
        let atom = Atom {
            serial: line.optional(6, 11, "invalid atom serial", |field| {
                decode_hybrid36(field, 5)
            })?,
            name: AtomName::new(name_columns.trim()),
            element,
            altloc: line.code(16)?.map_or(AltLoc::NONE, AltLoc::new),
            position: [
                line.number(30, 38, "invalid x coordinate")?,
                line.number(38, 46, "invalid y coordinate")?,
                line.number(46, 54, "invalid z coordinate")?,
            ],
            occupancy: line.number_or(54, 60, 1.0, "invalid occupancy")?,
            b_factor: line.number_or(60, 66, 0.0, "invalid temperature factor")?,
            charge: line.optional(78, 80, "invalid formal charge", parse_charge)?,
        };
        let segment = self.segment(&chain);
        let continues = self.residues.last().is_some_and(|last| {
            last.chain == chain
                && last.segment == segment
                && last.residue.auth_seq.get() == auth_seq
                && last.residue.ins_code == ins_code
                && last.residue.comp.as_str() == comp
        });
        if !continues {
            let mut residue =
                Residue::new(CompId::new(comp), None, AuthSeqId::new(auth_seq), ins_code);
            residue.hetero = hetero;
            self.residues.push(PdbResidue {
                segment,
                chain,
                residue,
            });
        }
        self.residues
            .last_mut()
            .expect("a residue was just ensured")
            .residue
            .atoms
            .push(atom);
        Ok(())
    }

    fn terminate(&mut self, line: &Line<'_>) -> Result<(), ReadError> {
        let chain = match line.code(21)? {
            Some(code) => AuthAsymId::new(code.to_string()),
            None => match self.residues.last() {
                Some(last) => last.chain.clone(),
                None => return Ok(()),
            },
        };
        *self.terminations.entry(chain).or_default() += 1;
        Ok(())
    }

    fn segment(&self, chain: &AuthAsymId) -> u32 {
        self.terminations.get(chain).copied().unwrap_or(0)
    }
}

#[derive(Clone, Copy, PartialEq, Eq, Hash)]
enum Role {
    Polymer,
    NonPolymer,
    Water,
}

fn roles(model: &PdbModel) -> Vec<Role> {
    let mut last_atom_record: HashMap<(&AuthAsymId, u32), usize> = HashMap::new();
    for (index, entry) in model.residues.iter().enumerate() {
        if !entry.residue.hetero {
            last_atom_record.insert((&entry.chain, entry.segment), index);
        }
    }
    model
        .residues
        .iter()
        .enumerate()
        .map(|(index, entry)| {
            let closed = entry.segment < model.segment(&entry.chain);
            let in_polymer = last_atom_record
                .get(&(&entry.chain, entry.segment))
                .is_some_and(|&last| closed || index <= last);
            if WATER_NAMES.contains(&entry.residue.comp.as_str()) {
                Role::Water
            } else if in_polymer {
                Role::Polymer
            } else {
                Role::NonPolymer
            }
        })
        .collect()
}

fn label_asym_name(index: usize) -> String {
    let mut name = String::new();
    let mut rest = index;
    loop {
        name.push((b'A' + (rest % 26) as u8) as char);
        if rest < 26 {
            return name;
        }
        rest = rest / 26 - 1;
    }
}

fn build_model(model: PdbModel) -> (Model, Vec<Role>) {
    let roles = roles(&model);
    let mut chains: Vec<(Role, Chain)> = Vec::new();
    let mut grouped_chains: HashMap<(Role, AuthAsymId, u32), usize> = HashMap::new();
    for (entry, role) in model.residues.into_iter().zip(&roles) {
        let group = match role {
            Role::Polymer => Some((*role, entry.chain.clone(), entry.segment)),
            Role::Water => Some((*role, entry.chain.clone(), 0)),
            Role::NonPolymer => None,
        };
        let index = match group {
            Some(key) => *grouped_chains.entry(key).or_insert_with(|| {
                chains.push((
                    *role,
                    Chain::new(LabelAsymId::new(""), entry.chain.clone(), None),
                ));
                chains.len() - 1
            }),
            None => {
                chains.push((
                    *role,
                    Chain::new(LabelAsymId::new(""), entry.chain.clone(), None),
                ));
                chains.len() - 1
            }
        };
        chains[index].1.residues.push(entry.residue);
    }
    let rank = |role: &Role| match role {
        Role::Polymer => 0,
        Role::NonPolymer => 1,
        Role::Water => 2,
    };
    chains.sort_by_key(|(role, _)| rank(role));
    let chain_roles = chains.iter().map(|(role, _)| *role).collect();
    let chains = chains
        .into_iter()
        .enumerate()
        .map(|(index, (_, chain))| Chain {
            label_asym: LabelAsymId::new(label_asym_name(index)),
            ..chain
        })
        .collect();
    (
        Model::new(ModelNumber::new(model.number), chains),
        chain_roles,
    )
}

fn assign_entities(model: &mut Model, roles: &[Role]) -> Vec<Entity> {
    let mut entities: Vec<Entity> = Vec::new();
    let mut by_content: HashMap<(u8, Vec<CompId>), usize> = HashMap::new();
    for (chain, role) in model.chains.iter_mut().zip(roles) {
        let (code, kind, content) = match role {
            Role::Polymer => (
                0,
                EntityKind::Polymer,
                chain.residues.iter().map(|r| r.comp.clone()).collect(),
            ),
            Role::NonPolymer => (
                1,
                EntityKind::NonPolymer,
                vec![chain.residues[0].comp.clone()],
            ),
            Role::Water => (2, EntityKind::Water, Vec::new()),
        };
        let index = *by_content.entry((code, content)).or_insert_with(|| {
            entities.push(Entity {
                id: EntityId::new((entities.len() + 1).to_string()),
                kind,
                description: None,
                polymer: None,
                branch: None,
                asym_ids: Vec::new(),
            });
            entities.len() - 1
        });
        entities[index].asym_ids.push(chain.label_asym.clone());
        chain.entity = Some(entities[index].id.clone());
    }
    entities
}

fn build_numbering(model: &Model, roles: &[Role]) -> Numbering {
    let mut numbering = Numbering::default();
    for (chain, role) in model.chains.iter().zip(roles) {
        for residue in &chain.residues {
            let auth = chain.auth_key(residue);
            match role {
                Role::Polymer => {
                    numbering.insert_unnumbered_polymer(chain.label_asym.clone(), auth)
                }
                _ => numbering.insert_non_polymer(chain.label_asym.clone(), auth),
            }
        }
    }
    numbering
}

pub fn parse_pdb(text: &str) -> Result<Structure, ReadError> {
    let mut entry_id = None;
    let mut finished: Vec<PdbModel> = Vec::new();
    let mut current: Option<PdbModel> = None;
    for (index, text) in text.lines().enumerate() {
        let line = Line {
            number: index + 1,
            text,
        };
        let record = line.columns(0, 6)?.trim_end();
        match record {
            "ATOM" | "HETATM" => current
                .get_or_insert_with(|| PdbModel::new(finished.len() as i32 + 1))
                .add_atom(&line)?,
            "MODEL" => {
                let number = line
                    .field(6, 80)?
                    .parse()
                    .map_err(|_| line.error("invalid model number"))?;
                finished.extend(current.replace(PdbModel::new(number)));
            }
            "ENDMDL" => finished.extend(current.take()),
            "TER" => {
                if let Some(model) = current.as_mut() {
                    model.terminate(&line)?;
                }
            }
            "HEADER" => entry_id = Some(line.field(62, 66)?.to_owned()).filter(|id| !id.is_empty()),
            "END" => break,
            _ => {}
        }
    }
    finished.extend(current.take());
    finished.retain(|model| !model.residues.is_empty());

    let mut models = Vec::with_capacity(finished.len());
    let mut entities = Vec::new();
    let mut numbering = Numbering::default();
    let mut label_entities: HashMap<LabelAsymId, EntityId> = HashMap::new();
    for (index, pdb_model) in finished.into_iter().enumerate() {
        let (mut model, roles) = build_model(pdb_model);
        if index == 0 {
            entities = assign_entities(&mut model, &roles);
            numbering = build_numbering(&model, &roles);
            label_entities = model
                .chains
                .iter()
                .filter_map(|chain| Some((chain.label_asym.clone(), chain.entity.clone()?)))
                .collect();
        } else {
            for chain in &mut model.chains {
                chain.entity = label_entities.get(&chain.label_asym).cloned();
            }
        }
        models.push(model);
    }
    Ok(Structure {
        entry_id,
        entities,
        models,
        operators: Vec::new(),
        assemblies: Vec::new(),
        numbering,
    })
}

pub fn read_pdb(path: impl AsRef<Path>) -> Result<Structure, ReadError> {
    let bytes = std::fs::read(path).map_err(mmcif_cif::Error::from)?;
    let text = String::from_utf8(bytes).map_err(|error| {
        let valid = &error.as_bytes()[..error.utf8_error().valid_up_to()];
        ReadError::Pdb {
            line: valid.iter().filter(|&&b| b == b'\n').count() + 1,
            problem: "text is not valid UTF-8",
        }
    })?;
    parse_pdb(&text)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn decodes_hybrid36_numbers() {
        assert_eq!(decode_hybrid36(" 9999", 5), Some(9999));
        assert_eq!(decode_hybrid36("  -12", 5), Some(-12));
        assert_eq!(decode_hybrid36("A0000", 5), Some(100_000));
        assert_eq!(decode_hybrid36("A000", 4), Some(10_000));
        assert_eq!(decode_hybrid36("ZZZZ", 4), Some(1_223_055));
        assert_eq!(decode_hybrid36("a000", 4), Some(1_223_056));
        assert_eq!(decode_hybrid36("A0a0", 4), None);
        assert_eq!(decode_hybrid36("    ", 4), None);
    }

    #[test]
    fn parses_formal_charges() {
        assert_eq!(parse_charge("2+"), Some(2));
        assert_eq!(parse_charge("1-"), Some(-1));
        assert_eq!(parse_charge(""), None);
    }

    #[test]
    fn infers_elements_from_atom_names() {
        assert_eq!(infer_element(" CA ", false).unwrap().symbol(), "C");
        assert_eq!(infer_element("HG21", false).unwrap().symbol(), "H");
        assert_eq!(infer_element("FE  ", true).unwrap().symbol(), "Fe");
        assert_eq!(infer_element("1HB ", false).unwrap().symbol(), "H");
    }

    #[test]
    fn names_label_chains_like_the_archive() {
        assert_eq!(label_asym_name(0), "A");
        assert_eq!(label_asym_name(25), "Z");
        assert_eq!(label_asym_name(26), "AA");
        assert_eq!(label_asym_name(27), "BA");
        assert_eq!(label_asym_name(26 + 26 * 26), "AAA");
    }
}
