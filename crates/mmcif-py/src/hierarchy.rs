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

/// One model of a structure: its chains, residues and atoms.
///
/// Per-atom data is available as lists and NumPy arrays in atom order.
/// ``len(model)`` is the number of chains, and two ``Model`` objects are equal
/// when they refer to the same model.
///
/// Examples
/// --------
/// >>> model = mmcif_io.read("1crn.bcif").models[0]
/// >>> model.number, model.chain_count, model.residue_count, model.atom_count
/// (1, 1, 46, 327)
/// >>> model.coordinates.shape
/// (327, 3)
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
    /// The model number, from ``_atom_site.pdbx_PDB_model_num``.
    #[getter]
    fn number(&self) -> i32 {
        self.inner.number().get()
    }

    /// The number of atoms in the model.
    #[getter]
    fn atom_count(&self) -> usize {
        self.inner.atom_count()
    }

    /// The number of residues in the model.
    #[getter]
    fn residue_count(&self) -> usize {
        self.inner.residue_count()
    }

    /// The number of chains in the model.
    #[getter]
    fn chain_count(&self) -> usize {
        self.inner.chain_count()
    }

    /// The chains of the model, in file order.
    #[getter]
    fn chains(&self) -> Vec<PyChain> {
        self.inner
            .chains()
            .map(|chain| self.chain_object(chain.index()))
            .collect()
    }

    /// Return the chain with label asym id ``label_asym`` that carries no assembly
    /// operators, or ``None``.
    ///
    /// Examples
    /// --------
    /// >>> model = mmcif_io.read("1crn.bcif").models[0]
    /// >>> model.chain("A")
    /// <Chain A (auth A) with 46 residues>
    /// >>> model.chain("Z") is None
    /// True
    #[pyo3(name = "chain")]
    fn find_chain(&self, label_asym: &str) -> Option<PyChain> {
        self.inner
            .chain(label_asym)
            .map(|chain| self.chain_object(chain.index()))
    }

    /// Return the residues of the model, in file order.
    ///
    /// Examples
    /// --------
    /// >>> model = mmcif_io.read("1crn.bcif").models[0]
    /// >>> residues = model.residues()
    /// >>> len(residues), residues[0].comp
    /// (46, 'THR')
    fn residues(&self) -> Vec<PyResidue> {
        self.inner
            .residues()
            .map(|residue| PyResidue::new(&self.inner, residue.index()))
            .collect()
    }

    /// Return the atoms of the model, in file order.
    ///
    /// Examples
    /// --------
    /// >>> model = mmcif_io.read("1crn.bcif").models[0]
    /// >>> model.atoms()[0]
    /// <Atom N N at [17.047, 14.099, 3.625]>
    fn atoms(&self) -> Vec<PyAtom> {
        self.inner
            .atoms()
            .map(|atom| PyAtom::new(&self.inner, atom.index()))
            .collect()
    }

    /// The atom positions, a float64 NumPy array of shape ``(atom_count, 3)``.
    #[getter]
    fn coordinates<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f64>>> {
        coordinates(py, self.inner.positions())
    }

    /// The element symbol of each atom.
    #[getter]
    fn elements(&self) -> Vec<String> {
        self.inner
            .elements()
            .iter()
            .map(|element| element.symbol().to_owned())
            .collect()
    }

    /// The name of each atom.
    #[getter]
    fn atom_names(&self) -> Vec<String> {
        self.inner
            .atoms()
            .map(|atom| atom.name().to_string())
            .collect()
    }

    /// The alternate location of each atom, or ``None`` for an atom without one.
    #[getter]
    fn altlocs(&self) -> Vec<Option<String>> {
        self.inner
            .altlocs()
            .iter()
            .map(|altloc| code(altloc.get()))
            .collect()
    }

    /// The occupancy of each atom, a float64 NumPy array.
    ///
    /// An atom whose file gives no occupancy has 1.0.
    #[getter]
    fn occupancies<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        PyArray1::from_slice(py, self.inner.occupancies())
    }

    /// The isotropic B-factor of each atom, a float64 NumPy array.
    ///
    /// An atom whose file gives no B-factor has NaN.
    #[getter]
    fn b_factors<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        PyArray1::from_slice(py, self.inner.b_factors())
    }

    /// The index in ``residues()`` of the residue each atom belongs to, a uint32
    /// NumPy array.
    #[getter]
    fn residue_indices<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<u32>> {
        PyArray1::from_slice(py, self.inner.atom_residue_indices())
    }

    /// The index in ``chains`` of the chain each atom belongs to, a uint32 NumPy
    /// array.
    #[getter]
    fn chain_indices<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<u32>> {
        let mut indices = vec![0u32; self.inner.atom_count()];
        for chain in self.inner.chains() {
            let index = chain.index().get() as u32;
            indices[chain.atom_range()].fill(index);
        }
        PyArray1::from_vec(py, indices)
    }

    /// Return a copy of the model that keeps one alternate location at each residue
    /// position.
    ///
    /// ``policy`` chooses the alternate location at each position:
    ///
    /// * ``"highest_occupancy"``: the one with the highest mean atom occupancy, taking
    ///   the first on a tie.
    /// * ``"first"``: the one that comes first among the atoms of the position.
    /// * a single character: that alternate location where the position has it, and
    ///   the first one otherwise.
    ///
    /// Atoms with no alternate location are always kept, and residues left with no
    /// atoms are removed.
    ///
    /// Raises
    /// ------
    /// ValueError
    ///     If ``policy`` is none of these.
    ///
    /// Examples
    /// --------
    /// >>> model = mmcif_io.read("example.cif").models[0]
    /// >>> model.atom_count
    /// 16
    /// >>> model.single_conformer().atom_count
    /// 13
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
        self.inner.chain_count()
    }

    fn __eq__(&self, other: &PyModel) -> bool {
        same_model(&self.inner, &other.inner)
    }

    fn __hash__(&self) -> isize {
        identity(&self.inner, 0)
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

/// A chain of a model.
///
/// ``len(chain)`` is the number of residues, and two ``Chain`` objects are equal
/// when they refer to the same chain of the same model.
///
/// Examples
/// --------
/// >>> chain = mmcif_io.read("1crn.bcif").models[0].chains[0]
/// >>> chain.label_asym, chain.auth_asym, chain.entity_id
/// ('A', 'A', '1')
/// >>> len(chain)
/// 46
#[pyclass(name = "Chain", module = "mmcif_io", frozen)]
pub struct PyChain {
    inner: ChainHandle,
}

#[pymethods]
impl PyChain {
    /// The position of the chain in ``Model.chains``.
    #[getter]
    fn index(&self) -> usize {
        self.inner.view().index().get()
    }

    /// The model the chain belongs to.
    #[getter]
    fn model(&self) -> PyModel {
        PyModel {
            inner: Arc::clone(self.inner.model()),
        }
    }

    /// The label asym id, from ``_atom_site.label_asym_id``.
    #[getter]
    fn label_asym(&self) -> String {
        self.inner.view().label_asym().to_string()
    }

    /// The author chain id, from ``_atom_site.auth_asym_id``.
    #[getter]
    fn auth_asym(&self) -> String {
        self.inner.view().auth_asym().to_string()
    }

    /// The id of the entity the chain belongs to, or ``None``.
    #[getter]
    fn entity_id(&self) -> Option<String> {
        self.inner.view().entity().map(|entity| entity.to_string())
    }

    /// The ids of the assembly operators that placed the chain, in the order they
    /// were applied. Empty for a chain read from a file.
    #[getter]
    fn operators(&self) -> Vec<String> {
        self.inner
            .view()
            .operators()
            .iter()
            .map(|operator| operator.to_string())
            .collect()
    }

    /// The residues of the chain, in file order.
    #[getter]
    fn residues(&self) -> Vec<PyResidue> {
        let model = self.inner.model();
        self.inner
            .view()
            .residues()
            .map(|residue| PyResidue::new(model, residue.index()))
            .collect()
    }

    /// The ``(start, end)`` range of the chain's atoms in the model's atom order,
    /// with ``end`` exclusive.
    #[getter]
    fn atom_range(&self) -> (usize, usize) {
        let range = self.inner.view().atom_range();
        (range.start, range.end)
    }

    /// The positions of the chain's atoms, a float64 NumPy array of shape
    /// ``(n, 3)``.
    #[getter]
    fn coordinates<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f64>>> {
        coordinates(py, self.inner.view().positions())
    }

    /// Return the residues at label sequence number ``seq``, in file order.
    ///
    /// The list holds more than one residue where the chain is heterogeneous at that
    /// position, and is empty if no residue has the number.
    ///
    /// Examples
    /// --------
    /// >>> chain = mmcif_io.read("example.cif").models[0].chain("AA")
    /// >>> chain.residues_at_label(2)
    /// [<Residue AA/2/100/SER>, <Residue AA/2/100/THR>]
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

    /// Return the residues at author residue number ``seq`` with insertion code
    /// ``ins_code``, in file order.
    ///
    /// Examples
    /// --------
    /// >>> chain = mmcif_io.read("example.cif").models[0].chain("AA")
    /// >>> chain.residues_at_auth(100, "A")
    /// [<Residue AA/3/100A/ALA>]
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

/// A residue of a chain.
///
/// ``len(residue)`` is the number of atoms, and two ``Residue`` objects are equal
/// when they refer to the same residue of the same model.
///
/// Examples
/// --------
/// >>> residue = mmcif_io.read("1crn.bcif").models[0].chains[0].residues[2]
/// >>> residue.comp, residue.label_seq, residue.auth_seq
/// ('CYS', 3, 3)
/// >>> residue.atom("SG").position
/// (15.885, 12.426, 10.016)
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
    /// The position of the residue in ``Model.residues()``.
    #[getter]
    fn index(&self) -> usize {
        self.inner.view().index().get()
    }

    /// The chain the residue belongs to.
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

    /// The chemical component id, from ``_atom_site.label_comp_id``.
    #[getter]
    fn comp(&self) -> String {
        self.inner.view().comp().to_string()
    }

    /// The label sequence number, from ``_atom_site.label_seq_id``, or ``None`` if
    /// the residue has none.
    #[getter]
    fn label_seq(&self) -> Option<u32> {
        self.inner.view().label_seq().map(LabelSeqId::get)
    }

    /// The author residue number, from ``_atom_site.auth_seq_id``.
    #[getter]
    fn auth_seq(&self) -> i32 {
        self.inner.view().auth_seq().get()
    }

    /// The insertion code, from ``_atom_site.pdbx_PDB_ins_code``, or ``None``.
    #[getter]
    fn ins_code(&self) -> Option<String> {
        code(self.inner.view().ins_code().get())
    }

    /// ``True`` if the residue's first atom was read from a ``HETATM`` record.
    #[getter]
    fn is_hetero(&self) -> bool {
        self.inner.view().is_hetero()
    }

    /// The atoms of the residue, in file order.
    #[getter]
    fn atoms(&self) -> Vec<PyAtom> {
        let model = self.inner.model();
        self.inner
            .view()
            .atoms()
            .map(|atom| PyAtom::new(model, atom.index()))
            .collect()
    }

    /// Return the first atom of the residue named ``name``, or ``None``.
    ///
    /// Examples
    /// --------
    /// >>> residue = mmcif_io.read("1crn.bcif").models[0].chains[0].residues[2]
    /// >>> residue.atom("CA").name
    /// 'CA'
    /// >>> residue.atom("OXT") is None
    /// True
    fn atom(&self, name: &str) -> Option<PyAtom> {
        let atom = self.inner.view().atom(name)?;
        Some(PyAtom::new(self.inner.model(), atom.index()))
    }

    /// The distinct alternate locations of the residue's atoms, in order of first
    /// appearance.
    #[getter]
    fn altlocs(&self) -> Vec<String> {
        self.inner
            .view()
            .altlocs()
            .into_iter()
            .filter_map(|altloc| code(altloc.get()))
            .collect()
    }

    /// The ``(start, end)`` range of the residue's atoms in the model's atom
    /// order, with ``end`` exclusive.
    #[getter]
    fn atom_range(&self) -> (usize, usize) {
        let range = self.inner.view().atom_range();
        (range.start, range.end)
    }

    /// The positions of the residue's atoms, a float64 NumPy array of shape
    /// ``(n, 3)``.
    #[getter]
    fn coordinates<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f64>>> {
        coordinates(py, self.inner.view().positions())
    }

    /// The residue key ``(label_asym, label_seq, auth_seq, ins_code, comp)``.
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

    /// The author address ``(auth_asym, auth_seq, ins_code)``.
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

/// An atom of a residue.
///
/// Two ``Atom`` objects are equal when they refer to the same atom of the same
/// model.
///
/// Examples
/// --------
/// >>> atom = mmcif_io.read("1crn.bcif").models[0].atoms()[1]
/// >>> atom.name, atom.element, atom.serial
/// ('CA', 'C', 2)
/// >>> atom.position
/// (16.967, 12.784, 4.338)
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
    /// The position of the atom in the model's atom order.
    #[getter]
    fn index(&self) -> usize {
        self.inner.view().index().get()
    }

    /// The residue the atom belongs to.
    #[getter]
    fn residue(&self) -> PyResidue {
        PyResidue::new(self.inner.model(), self.inner.view().residue().index())
    }

    /// The atom serial number, from ``_atom_site.id``, or ``None``.
    #[getter]
    fn serial(&self) -> Option<i64> {
        self.inner.view().serial()
    }

    /// The atom name, from ``_atom_site.label_atom_id``.
    #[getter]
    fn name(&self) -> String {
        self.inner.view().name().to_string()
    }

    /// The element symbol, from ``_atom_site.type_symbol``.
    #[getter]
    fn element(&self) -> String {
        self.inner.view().element().symbol().to_owned()
    }

    /// The alternate location, from ``_atom_site.label_alt_id``, or ``None``.
    #[getter]
    fn altloc(&self) -> Option<String> {
        code(self.inner.view().altloc().get())
    }

    /// The ``(x, y, z)`` position.
    #[getter]
    fn position(&self) -> (f64, f64, f64) {
        let [x, y, z] = self.inner.view().position();
        (x, y, z)
    }

    /// The occupancy, or 1.0 if the file gives none.
    #[getter]
    fn occupancy(&self) -> f64 {
        self.inner.view().occupancy()
    }

    /// The isotropic B-factor, from ``_atom_site.B_iso_or_equiv``, or ``None`` if
    /// the file gives none.
    #[getter]
    fn b_factor(&self) -> Option<f64> {
        self.inner.view().b_factor()
    }

    /// The formal charge, from ``_atom_site.pdbx_formal_charge``, or ``None``.
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
