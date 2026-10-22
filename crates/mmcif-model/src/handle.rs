use std::sync::Arc;

use crate::model::{AtomIndex, ChainIndex, Model, ResidueIndex};
use crate::view::{AtomRef, ChainRef, ResidueRef};

/// An owned reference to a chain, holding a shared pointer to its [`Model`].
///
/// # Examples
///
/// ```
/// use std::sync::Arc;
///
/// use mmcif_model::ChainHandle;
///
/// # let structure = mmcif_model::parse_mmcif(include_str!("../tests/fixtures/pathologies.cif"))?;
/// let model = Arc::clone(structure.first_model().unwrap());
/// let index = model.chain("C").unwrap().index();
///
/// let handle = ChainHandle::new(model, index).unwrap();
/// assert_eq!(handle.view().label_asym().as_str(), "C");
/// # Ok::<(), mmcif_model::ReadError>(())
/// ```
#[derive(Debug, Clone)]
pub struct ChainHandle {
    model: Arc<Model>,
    index: ChainIndex,
}

/// An owned reference to a residue, holding a shared pointer to its
/// [`Model`].
///
/// # Examples
///
/// ```
/// use std::sync::Arc;
///
/// use mmcif_model::ResidueHandle;
///
/// # let structure = mmcif_model::parse_mmcif(include_str!("../tests/fixtures/pathologies.cif"))?;
/// let model = Arc::clone(structure.first_model().unwrap());
/// let index = model.chain("C").unwrap().residues().next().unwrap().index();
///
/// let handle = ResidueHandle::new(model, index).unwrap();
/// assert_eq!(handle.view().comp().as_str(), "ZN");
/// # Ok::<(), mmcif_model::ReadError>(())
/// ```
#[derive(Debug, Clone)]
pub struct ResidueHandle {
    model: Arc<Model>,
    index: ResidueIndex,
}

/// An owned reference to an atom, holding a shared pointer to its [`Model`].
///
/// # Examples
///
/// ```
/// use std::sync::Arc;
///
/// use mmcif_model::AtomHandle;
///
/// # let structure = mmcif_model::parse_mmcif(include_str!("../tests/fixtures/pathologies.cif"))?;
/// let model = Arc::clone(structure.first_model().unwrap());
/// let index = model.chain("C").unwrap().atoms().next().unwrap().index();
///
/// let handle = AtomHandle::new(model, index).unwrap();
/// assert_eq!(handle.view().position(), [0.0, 0.0, 5.0]);
/// # Ok::<(), mmcif_model::ReadError>(())
/// ```
#[derive(Debug, Clone)]
pub struct AtomHandle {
    model: Arc<Model>,
    index: AtomIndex,
}

impl ChainHandle {
    /// Creates a handle to the chain at `index` in `model`, or returns `None`
    /// if `index` is out of range.
    pub fn new(model: Arc<Model>, index: ChainIndex) -> Option<Self> {
        model.chain_at(index)?;
        Some(ChainHandle { model, index })
    }

    /// Returns the model the chain belongs to.
    pub fn model(&self) -> &Arc<Model> {
        &self.model
    }

    /// Returns a view of the chain.
    pub fn view(&self) -> ChainRef<'_> {
        ChainRef::new(&self.model, self.index.0)
    }
}

impl ResidueHandle {
    /// Creates a handle to the residue at `index` in `model`, or returns
    /// `None` if `index` is out of range.
    pub fn new(model: Arc<Model>, index: ResidueIndex) -> Option<Self> {
        model.residue_at(index)?;
        Some(ResidueHandle { model, index })
    }

    /// Returns the model the residue belongs to.
    pub fn model(&self) -> &Arc<Model> {
        &self.model
    }

    /// Returns a view of the residue.
    pub fn view(&self) -> ResidueRef<'_> {
        ResidueRef::new(&self.model, self.index.0)
    }
}

impl AtomHandle {
    /// Creates a handle to the atom at `index` in `model`, or returns `None`
    /// if `index` is out of range.
    pub fn new(model: Arc<Model>, index: AtomIndex) -> Option<Self> {
        model.atom_at(index)?;
        Some(AtomHandle { model, index })
    }

    /// Returns the model the atom belongs to.
    pub fn model(&self) -> &Arc<Model> {
        &self.model
    }

    /// Returns a view of the atom.
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
