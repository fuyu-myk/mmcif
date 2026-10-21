use std::sync::Arc;

use mmcif_model::{
    AltLoc, AltlocPolicy, AtomHandle, AuthSeqId, ChainHandle, InsCode, LabelSeqId, Model,
    ResidueHandle,
};
use numpy::{PyArray1, PyArray2, PyArrayMethods};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

fn code(value: Option<char>) -> Option<String> {
    value.map(String::from)
}

fn coordinates<'py>(
    py: Python<'py>,
    positions: &[[f64; 3]],
) -> PyResult<Bound<'py, PyArray2<f64>>> {
    PyArray1::from_slice(py, positions.as_flattened()).reshape([positions.len(), 3])
}

fn ins_code(value: Option<char>) -> InsCode {
    value.map_or(InsCode::NONE, InsCode::new)
}

fn same_model(a: &Arc<Model>, b: &Arc<Model>) -> bool {
    Arc::ptr_eq(a, b)
}

fn identity(model: &Arc<Model>, index: usize) -> isize {
    (Arc::as_ptr(model) as usize as isize)
        .wrapping_mul(31)
        .wrapping_add(index as isize)
}

fn parse_policy(policy: &str) -> PyResult<AltlocPolicy> {
    match policy {
        "first" => Ok(AltlocPolicy::First),
        "highest_occupancy" => Ok(AltlocPolicy::HighestOccupancy),
        other => {
            let mut chars = other.chars();
            match (chars.next(), chars.next()) {
                (Some(altloc), None) => Ok(AltlocPolicy::Prefer(AltLoc::new(altloc))),
                _ => Err(PyValueError::new_err(format!(
                    "policy must be 'first', 'highest_occupancy' or a single altloc character, not {other:?}"
                ))),
            }
        }
    }
}

#[pyclass(name = "Model", module = "mmcif_io", frozen)]
pub struct PyModel {
    pub(crate) inner: Arc<Model>,
}

impl PyModel {
    fn chain_object(&self, index: mmcif_model::ChainIndex) -> PyChain {
        PyChain {
            inner: ChainHandle::new(Arc::clone(&self.inner), index)
                .expect("index comes from this model"),
        }
    }
}

#[pymethods]
impl PyModel {
    #[getter]
    fn number(&self) -> i32 {
        self.inner.number().get()
    }

    #[getter]
    fn atom_count(&self) -> usize {
        self.inner.atom_count()
    }

    #[getter]
    fn residue_count(&self) -> usize {
        self.inner.residue_count()
    }

    #[getter]
    fn chain_count(&self) -> usize {
        self.inner.chain_count()
    }

    #[getter]
    fn chains(&self) -> Vec<PyChain> {
        self.inner
            .chains()
            .map(|chain| self.chain_object(chain.index()))
            .collect()
    }

    #[pyo3(name = "chain")]
    fn find_chain(&self, label_asym: &str) -> Option<PyChain> {
        self.inner
            .chain(label_asym)
            .map(|chain| self.chain_object(chain.index()))
    }

    #[getter]
    fn residues(&self) -> Vec<PyResidue> {
        self.inner
            .residues()
            .map(|residue| PyResidue::new(&self.inner, residue.index()))
            .collect()
    }

    fn atoms(&self) -> Vec<PyAtom> {
        self.inner
            .atoms()
            .map(|atom| PyAtom::new(&self.inner, atom.index()))
            .collect()
    }

    #[getter]
    fn coordinates<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f64>>> {
        coordinates(py, self.inner.positions())
    }

    #[getter]
    fn elements(&self) -> Vec<String> {
        self.inner
            .elements()
            .iter()
            .map(|element| element.symbol().to_owned())
            .collect()
    }

    #[getter]
    fn atom_names(&self) -> Vec<String> {
        self.inner
            .atoms()
            .map(|atom| atom.name().to_string())
            .collect()
    }

    #[getter]
    fn altlocs(&self) -> Vec<Option<String>> {
        self.inner
            .altlocs()
            .iter()
            .map(|altloc| code(altloc.get()))
            .collect()
    }

    #[getter]
    fn occupancies<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        PyArray1::from_slice(py, self.inner.occupancies())
    }

    #[getter]
    fn b_factors<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        PyArray1::from_slice(py, self.inner.b_factors())
    }

    #[getter]
    fn residue_indices<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<u32>> {
        let mut indices = vec![0u32; self.inner.atom_count()];
        for residue in self.inner.residues() {
            let index = residue.index().get() as u32;
            indices[residue.atom_range()].fill(index);
        }
        PyArray1::from_vec(py, indices)
    }

    #[getter]
    fn chain_indices<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<u32>> {
        let mut indices = vec![0u32; self.inner.atom_count()];
        for chain in self.inner.chains() {
            let index = chain.index().get() as u32;
            indices[chain.atom_range()].fill(index);
        }
        PyArray1::from_vec(py, indices)
    }

    #[pyo3(signature = (policy = "highest_occupancy"))]
    fn single_conformer(&self, py: Python<'_>, policy: &str) -> PyResult<PyModel> {
        let policy = parse_policy(policy)?;
        let model = Arc::clone(&self.inner);
        let selected = py.detach(move || model.single_conformer(policy));
        Ok(PyModel {
            inner: Arc::new(selected),
        })
    }

    fn __len__(&self) -> usize {
        self.inner.atom_count()
    }

    fn __repr__(&self) -> String {
        format!(
            "<Model {} with {} chains, {} atoms>",
            self.inner.number(),
            self.inner.chain_count(),
            self.inner.atom_count()
        )
    }
}

#[pyclass(name = "Chain", module = "mmcif_io", frozen)]
pub struct PyChain {
    inner: ChainHandle,
}

#[pymethods]
impl PyChain {
    #[getter]
    fn index(&self) -> usize {
        self.inner.view().index().get()
    }

    #[getter]
    fn model(&self) -> PyModel {
        PyModel {
            inner: Arc::clone(self.inner.model()),
        }
    }

    #[getter]
    fn label_asym(&self) -> String {
        self.inner.view().label_asym().to_string()
    }

    #[getter]
    fn auth_asym(&self) -> String {
        self.inner.view().auth_asym().to_string()
    }

    #[getter]
    fn entity_id(&self) -> Option<String> {
        self.inner.view().entity().map(|entity| entity.to_string())
    }

    #[getter]
    fn operators(&self) -> Vec<String> {
        self.inner
            .view()
            .operators()
            .iter()
            .map(|operator| operator.to_string())
            .collect()
    }

    #[getter]
    fn residues(&self) -> Vec<PyResidue> {
        let model = self.inner.model();
        self.inner
            .view()
            .residues()
            .map(|residue| PyResidue::new(model, residue.index()))
            .collect()
    }

    #[getter]
    fn atom_range(&self) -> (usize, usize) {
        let range = self.inner.view().atom_range();
        (range.start, range.end)
    }

    #[getter]
    fn coordinates<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f64>>> {
        coordinates(py, self.inner.view().positions())
    }

    fn residues_at_label(&self, seq: u32) -> Vec<PyResidue> {
        let Some(seq) = LabelSeqId::new(seq) else {
            return Vec::new();
        };
        let model = self.inner.model();
        self.inner
            .view()
            .residues_at_label(seq)
            .map(|residue| PyResidue::new(model, residue.index()))
            .collect()
    }

    #[pyo3(signature = (seq, ins_code = None))]
    fn residues_at_auth(&self, seq: i32, ins_code: Option<char>) -> Vec<PyResidue> {
        let model = self.inner.model();
        self.inner
            .view()
            .residues_at_auth(AuthSeqId::new(seq), self::ins_code(ins_code))
            .map(|residue| PyResidue::new(model, residue.index()))
            .collect()
    }

    fn __len__(&self) -> usize {
        self.inner.view().residue_count()
    }

    fn __eq__(&self, other: &PyChain) -> bool {
        same_model(self.inner.model(), other.inner.model()) && self.index() == other.index()
    }

    fn __hash__(&self) -> isize {
        identity(self.inner.model(), self.index())
    }

    fn __repr__(&self) -> String {
        let view = self.inner.view();
        format!(
            "<Chain {} (auth {}) with {} residues>",
            view.label_asym(),
            view.auth_asym(),
            view.residue_count()
        )
    }
}

#[pyclass(name = "Residue", module = "mmcif_io", frozen)]
pub struct PyResidue {
    inner: ResidueHandle,
}

impl PyResidue {
    fn new(model: &Arc<Model>, index: mmcif_model::ResidueIndex) -> PyResidue {
        PyResidue {
            inner: ResidueHandle::new(Arc::clone(model), index)
                .expect("index comes from this model"),
        }
    }
}

#[pymethods]
impl PyResidue {
    #[getter]
    fn index(&self) -> usize {
        self.inner.view().index().get()
    }

    #[getter]
    fn chain(&self) -> PyChain {
        PyChain {
            inner: ChainHandle::new(
                Arc::clone(self.inner.model()),
                self.inner.view().chain().index(),
            )
            .expect("index comes from this model"),
        }
    }

    #[getter]
    fn comp(&self) -> String {
        self.inner.view().comp().to_string()
    }

    #[getter]
    fn label_seq(&self) -> Option<u32> {
        self.inner.view().label_seq().map(LabelSeqId::get)
    }

    #[getter]
    fn auth_seq(&self) -> i32 {
        self.inner.view().auth_seq().get()
    }

    #[getter]
    fn ins_code(&self) -> Option<String> {
        code(self.inner.view().ins_code().get())
    }

    #[getter]
    fn is_hetero(&self) -> bool {
        self.inner.view().is_hetero()
    }

    #[getter]
    fn atoms(&self) -> Vec<PyAtom> {
        let model = self.inner.model();
        self.inner
            .view()
            .atoms()
            .map(|atom| PyAtom::new(model, atom.index()))
            .collect()
    }

    fn atom(&self, name: &str) -> Option<PyAtom> {
        let atom = self.inner.view().atom(name)?;
        Some(PyAtom::new(self.inner.model(), atom.index()))
    }

    #[getter]
    fn altlocs(&self) -> Vec<String> {
        self.inner
            .view()
            .altlocs()
            .into_iter()
            .filter_map(|altloc| code(altloc.get()))
            .collect()
    }

    #[getter]
    fn atom_range(&self) -> (usize, usize) {
        let range = self.inner.view().atom_range();
        (range.start, range.end)
    }

    #[getter]
    fn coordinates<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f64>>> {
        coordinates(py, self.inner.view().positions())
    }

    #[getter]
    fn key(&self) -> (String, Option<u32>, i32, Option<String>, String) {
        let key = self.inner.view().key();
        (
            key.asym.to_string(),
            key.label_seq.map(LabelSeqId::get),
            key.auth_seq.get(),
            code(key.ins_code.get()),
            key.comp.to_string(),
        )
    }

    #[getter]
    fn auth_key(&self) -> (String, i32, Option<String>) {
        let key = self.inner.view().auth_key();
        (
            key.chain.to_string(),
            key.seq.get(),
            code(key.ins_code.get()),
        )
    }

    fn __len__(&self) -> usize {
        self.inner.view().atom_count()
    }

    fn __eq__(&self, other: &PyResidue) -> bool {
        same_model(self.inner.model(), other.inner.model()) && self.index() == other.index()
    }

    fn __hash__(&self) -> isize {
        identity(self.inner.model(), self.index())
    }

    fn __repr__(&self) -> String {
        format!("<Residue {}>", self.inner.view().key())
    }
}

#[pyclass(name = "Atom", module = "mmcif_io", frozen)]
pub struct PyAtom {
    inner: AtomHandle,
}

impl PyAtom {
    fn new(model: &Arc<Model>, index: mmcif_model::AtomIndex) -> PyAtom {
        PyAtom {
            inner: AtomHandle::new(Arc::clone(model), index).expect("index comes from this model"),
        }
    }
}

#[pymethods]
impl PyAtom {
    #[getter]
    fn index(&self) -> usize {
        self.inner.view().index().get()
    }

    #[getter]
    fn residue(&self) -> PyResidue {
        PyResidue::new(self.inner.model(), self.inner.view().residue().index())
    }

    #[getter]
    fn serial(&self) -> Option<i64> {
        self.inner.view().serial()
    }

    #[getter]
    fn name(&self) -> String {
        self.inner.view().name().to_string()
    }

    #[getter]
    fn element(&self) -> String {
        self.inner.view().element().symbol().to_owned()
    }

    #[getter]
    fn altloc(&self) -> Option<String> {
        code(self.inner.view().altloc().get())
    }

    #[getter]
    fn position(&self) -> (f64, f64, f64) {
        let [x, y, z] = self.inner.view().position();
        (x, y, z)
    }

    #[getter]
    fn occupancy(&self) -> f64 {
        self.inner.view().occupancy()
    }

    #[getter]
    fn b_factor(&self) -> f64 {
        self.inner.view().b_factor()
    }

    #[getter]
    fn charge(&self) -> Option<i8> {
        self.inner.view().charge()
    }

    fn __eq__(&self, other: &PyAtom) -> bool {
        same_model(self.inner.model(), other.inner.model()) && self.index() == other.index()
    }

    fn __hash__(&self) -> isize {
        identity(self.inner.model(), self.index())
    }

    fn __repr__(&self) -> String {
        let view = self.inner.view();
        format!(
            "<Atom {} {} at {:?}>",
            view.name(),
            view.element(),
            view.position()
        )
    }
}
