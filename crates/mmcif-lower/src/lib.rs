//! Backend input written from a design. The output format tracks each backend and may change until 1.0.

mod boltzgen;
mod capability;
mod ranges;
mod rfd3;
mod slots;
mod text;

use mmcif_spec::{Design, Diagnostic, SourceId};

pub use capability::Capabilities;

const ALL: [Backend; 2] = [Backend::BoltzGen, Backend::Rfd3];

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum Backend {
    BoltzGen,
    Rfd3,
}

#[derive(Debug, thiserror::Error)]
#[error("{} cannot express this specification", backend.name())]
pub struct Refusal {
    pub backend: Backend,
    problems: Vec<Diagnostic>,
}

impl Refusal {
    pub fn problems(&self) -> &[Diagnostic] {
        &self.problems
    }
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

    pub fn names() -> Vec<&'static str> {
        ALL.iter().map(|backend| backend.name()).collect()
    }

    pub fn all() -> &'static [Backend] {
        &ALL
    }

    pub fn schema(self) -> &'static str {
        self.capabilities().schema
    }
}

pub fn lower(backend: Backend, source: SourceId, design: &Design) -> Result<String, Refusal> {
    let emitted = match backend {
        Backend::BoltzGen => boltzgen::emit(source, design),
        Backend::Rfd3 => rfd3::emit(source, design),
    };
    emitted.map_err(|problems| Refusal { backend, problems })
}
