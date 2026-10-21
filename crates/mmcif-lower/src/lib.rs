mod boltzgen;
mod capability;
mod ranges;
mod rfd3;

use mmcif_spec::{Design, Diagnostic, SourceId};

pub use capability::Capabilities;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Backend {
    BoltzGen,
    Rfd3,
}

impl Backend {
    pub fn named(name: &str) -> Option<Backend> {
        match name {
            "boltzgen" => Some(Backend::BoltzGen),
            "rfd3" => Some(Backend::Rfd3),
            _ => None,
        }
    }

    pub fn name(self) -> &'static str {
        self.capabilities().name
    }

    pub fn capabilities(self) -> &'static Capabilities {
        match self {
            Backend::BoltzGen => &boltzgen::CAPABILITIES,
            Backend::Rfd3 => &rfd3::CAPABILITIES,
        }
    }

    pub fn names() -> [&'static str; 2] {
        ["boltzgen", "rfd3"]
    }
}

pub fn lower(
    backend: Backend,
    source: SourceId,
    design: &Design,
) -> Result<String, Vec<Diagnostic>> {
    match backend {
        Backend::BoltzGen => boltzgen::emit(source, design),
        Backend::Rfd3 => rfd3::emit(source, design),
    }
}
