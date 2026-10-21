use std::sync::Arc;

use crate::model::{AtomIndex, ChainIndex, Model, ResidueIndex};
use crate::view::{AtomRef, ChainRef, ResidueRef};

#[derive(Debug, Clone)]
pub struct ChainHandle {
    model: Arc<Model>,
    index: ChainIndex,
}

#[derive(Debug, Clone)]
pub struct ResidueHandle {
    model: Arc<Model>,
    index: ResidueIndex,
}

#[derive(Debug, Clone)]
pub struct AtomHandle {
    model: Arc<Model>,
    index: AtomIndex,
}

impl ChainHandle {
    pub fn new(model: Arc<Model>, index: ChainIndex) -> Option<Self> {
        model.chain_at(index)?;
        Some(ChainHandle { model, index })
    }

    pub fn model(&self) -> &Arc<Model> {
        &self.model
    }

    pub fn view(&self) -> ChainRef<'_> {
        ChainRef::new(&self.model, self.index.0)
    }
}

impl ResidueHandle {
    pub fn new(model: Arc<Model>, index: ResidueIndex) -> Option<Self> {
        model.residue_at(index)?;
        Some(ResidueHandle { model, index })
    }

    pub fn model(&self) -> &Arc<Model> {
        &self.model
    }

    pub fn view(&self) -> ResidueRef<'_> {
        ResidueRef::new(&self.model, self.index.0)
    }
}

impl AtomHandle {
    pub fn new(model: Arc<Model>, index: AtomIndex) -> Option<Self> {
        model.atom_at(index)?;
        Some(AtomHandle { model, index })
    }

    pub fn model(&self) -> &Arc<Model> {
        &self.model
    }

    pub fn view(&self) -> AtomRef<'_> {
        AtomRef::new(&self.model, self.index.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parse_mmcif;

    #[test]
    fn handles_outlive_the_structure_and_reject_bad_indices() {
        let text = "data_x\nloop_\n_atom_site.label_asym_id\n_atom_site.label_comp_id\n_atom_site.label_atom_id\n\
                    _atom_site.type_symbol\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n\
                    _atom_site.Cartn_z\nA GLY N N 1 1 2 3\nA GLY CA C 1 4 5 6\n";
        let model = Arc::clone(&parse_mmcif(text).unwrap().models()[0]);
        let atom = AtomHandle::new(Arc::clone(&model), AtomIndex(1)).unwrap();
        let residue = ResidueHandle::new(Arc::clone(&model), ResidueIndex(0)).unwrap();
        let chain = ChainHandle::new(model, ChainIndex(0)).unwrap();
        assert_eq!(atom.view().name().as_str(), "CA");
        assert_eq!(atom.view().position(), [4.0, 5.0, 6.0]);
        assert_eq!(residue.view().atom_count(), 2);
        assert_eq!(chain.view().label_asym().as_str(), "A");
        assert!(AtomHandle::new(Arc::clone(chain.model()), AtomIndex(2)).is_none());
    }
}
