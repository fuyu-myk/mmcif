"""Read mmCIF, BinaryCIF and PDB files into a typed structure model.

``read`` returns a ``Structure`` holding models, chains, residues and
atoms, with coordinates as NumPy arrays. ``Document`` gives the CIF
tables of a file, ``NeighbourGrid`` answers distance queries over
coordinates, and ``CcdIndex`` reads a compiled Chemical Component
Dictionary.

The examples assume ``import mmcif_io``.

Examples
--------
>>> structure = mmcif_io.read("1crn.bcif")
>>> model = structure.models[0]
>>> model.atom_count
327
>>> model.coordinates.shape
(327, 3)
"""

from collections.abc import Sequence
from os import PathLike
from typing import final

import numpy as np
from numpy.typing import NDArray

_Path = str | PathLike[str]
_Point = Sequence[float] | NDArray[np.float64]
_ResidueKey = tuple[str, int | None, int, str | None, str, tuple[str, ...]]
_AuthAddress = tuple[str, int, str | None]
_AtomRow = tuple[str, str, int, bool, bool, str]
_BondRow = tuple[int, int, str, bool, str]

__all__ = [
    "ReadError",
    "NumberingError",
    "AssemblyError",
    "Document",
    "Block",
    "Category",
    "Structure",
    "Entity",
    "Assembly",
    "Model",
    "Chain",
    "Residue",
    "Atom",
    "NeighbourGrid",
    "CcdIndex",
    "Component",
    "read",
    "read_mmcif",
    "read_bcif",
    "read_pdb",
    "parse_mmcif",
    "parse_bcif",
    "parse_pdb",
]

class ReadError(ValueError):
    """Raised when a file, text or bytes cannot be read as mmCIF, BinaryCIF, PDB
    or a component index.

    A subclass of ``ValueError``.
    """
class NumberingError(ValueError):
    """Raised when a residue address has no counterpart in the other numbering
    scheme.

    A subclass of ``ValueError``.
    """
class AssemblyError(ValueError):
    """Raised when a biological assembly cannot be built.

    A subclass of ``ValueError``.
    """

def read(path: _Path, *, first_model_only: bool = False, lenient: bool = False) -> Structure:
    """Read the structure in the file at ``path``.

    The format follows the file extension: ``.cif`` and ``.mmcif`` are read as
    mmCIF, ``.bcif`` as BinaryCIF, and ``.pdb``, ``.ent`` and assembly files such
    as ``.pdb1`` as PDB. A gzip-compressed file, such as ``1crn.cif.gz``, is
    decompressed. ``path`` is a string or a path-like object.

    With ``first_model_only``, only the first model is read. With ``lenient``,
    a malformed optional category is skipped instead of raising, and its error is
    listed in ``Structure.skipped``. The optional categories are the entities,
    the numbering schemes, the operators and the assemblies, and in a PDB file the
    ``REMARK 350`` assemblies.

    Raises
    ------
    OSError
        If the file cannot be read.
    ReadError
        If the extension is not one of these, or the contents are not valid for
        the format.

    Examples
    --------
    >>> mmcif_io.read("1crn.bcif")
    <Structure 1CRN with 1 models>
    >>> mmcif_io.read("example.pdb").entry_id
    '1ABC'
    >>> len(mmcif_io.read("example.pdb", first_model_only=True).models)
    1
    """
def read_mmcif(
    path: _Path, *, first_model_only: bool = False, lenient: bool = False
) -> Structure:
    """Read the structure in the mmCIF file at ``path``.

    ``first_model_only`` and ``lenient`` work as in ``read``.

    Raises
    ------
    OSError
        If the file cannot be read.
    ReadError
        If the contents are not valid for the format, or lack the categories a
        structure needs.

    Examples
    --------
    >>> structure = mmcif_io.read_mmcif("1crn.cif")
    >>> structure.entry_id
    '1CRN'
    """
def read_bcif(
    path: _Path, *, first_model_only: bool = False, lenient: bool = False
) -> Structure:
    """Read the structure in the BinaryCIF file at ``path``.

    ``first_model_only`` and ``lenient`` work as in ``read``.

    Raises
    ------
    OSError
        If the file cannot be read.
    ReadError
        If the contents are not valid for the format, or lack the categories a
        structure needs.

    Examples
    --------
    >>> structure = mmcif_io.read_bcif("1crn.bcif")
    >>> structure.models[0].atom_count
    327
    """
def read_pdb(
    path: _Path, *, first_model_only: bool = False, lenient: bool = False
) -> Structure:
    """Read the structure in the PDB file at ``path``.

    ``first_model_only`` and ``lenient`` work as in ``read``.

    Raises
    ------
    OSError
        If the file cannot be read.
    ReadError
        If the contents are not valid for the format, or lack the categories a
        structure needs.

    Examples
    --------
    >>> structure = mmcif_io.read_pdb("example.pdb")
    >>> len(structure.models)
    2
    """
def parse_mmcif(
    text: str, *, first_model_only: bool = False, lenient: bool = False
) -> Structure:
    """Parse a structure from mmCIF text.

    ``first_model_only`` and ``lenient`` work as in ``read``.

    Raises
    ------
    ReadError
        If ``text`` is not valid mmCIF, or lacks the categories a structure
        needs.

    Examples
    --------
    >>> with open("1crn.cif") as file:
    ...     structure = mmcif_io.parse_mmcif(file.read())
    >>> structure.entry_id
    '1CRN'
    """
def parse_bcif(
    data: bytes, *, first_model_only: bool = False, lenient: bool = False
) -> Structure:
    """Parse a structure from BinaryCIF bytes.

    ``first_model_only`` and ``lenient`` work as in ``read``.

    Raises
    ------
    ReadError
        If ``data`` is not valid BinaryCIF, or lacks the categories a structure
        needs.

    Examples
    --------
    >>> with open("1crn.bcif", "rb") as file:
    ...     structure = mmcif_io.parse_bcif(file.read())
    >>> structure.entry_id
    '1CRN'
    """
def parse_pdb(
    text: str, *, first_model_only: bool = False, lenient: bool = False
) -> Structure:
    """Parse a structure from PDB text.

    ``first_model_only`` and ``lenient`` work as in ``read``.

    Raises
    ------
    ReadError
        If ``text`` is not valid PDB.

    Examples
    --------
    >>> with open("example.pdb") as file:
    ...     structure = mmcif_io.parse_pdb(file.read())
    >>> structure.entry_id
    '1ABC'
    """
@final
class Structure:
    """A macromolecular structure: its models, entities, biological assemblies and
    residue numbering.

    Read one with ``read`` or one of the other readers.

    Examples
    --------
    >>> structure = mmcif_io.read("1crn.bcif")
    >>> structure.entry_id
    '1CRN'
    >>> structure.models
    [<Model 1 with 1 chains, 327 atoms>]
    >>> structure.entities
    [<Entity 1 polymer>]
    """
    @staticmethod
    def from_document(document: Document) -> Structure:
        """Build a structure from a parsed CIF document.

        The first data block of ``document`` is read.

        Raises
        ------
        ReadError
            If the document has no data block, or lacks the categories a structure
            needs.

        Examples
        --------
        >>> document = mmcif_io.Document.open("1crn.cif")
        >>> mmcif_io.Structure.from_document(document).entry_id
        '1CRN'
        """
    @property
    def entry_id(self) -> str | None:
        """The ``_entry.id`` of the structure, or ``None`` if the file gives none."""
    @property
    def skipped(self) -> list[str]:
        """The errors of the malformed optional categories that a lenient read
        skipped, as messages, in the order they were read. Empty unless the
        structure was read with ``lenient=True``.

        Examples
        --------
        >>> mmcif_io.read("1crn.bcif", lenient=True).skipped
        []
        """
    @property
    def models(self) -> list[Model]:
        """The models of the structure, in file order."""
    def model(self, number: int) -> Model | None:
        """Return the model with model number ``number``, or ``None``.

        Examples
        --------
        >>> structure = mmcif_io.read("example.pdb")
        >>> structure.model(2).number
        2
        >>> structure.model(3) is None
        True
        """
    @property
    def entities(self) -> list[Entity]:
        """The entities of the structure, in file order."""
    def entity(self, id: str) -> Entity | None:
        """Return the entity with id ``id``, or ``None``.

        Examples
        --------
        >>> structure = mmcif_io.read("1crn.bcif")
        >>> structure.entity("1").description
        'CRAMBIN'
        >>> structure.entity("9") is None
        True
        """
    @property
    def assemblies(self) -> list[Assembly]:
        """The biological assemblies of the structure, in file order."""
    def assembly(self, id: str) -> Assembly | None:
        """Return the assembly with id ``id``, or ``None``.

        Examples
        --------
        >>> structure = mmcif_io.read("1crn.bcif")
        >>> structure.assembly("1").details
        'author_defined_assembly'
        >>> structure.assembly("9") is None
        True
        """
    def assemble(self, assembly_id: str, model_index: int = 0) -> Model:
        """Return a new model holding the assembly ``assembly_id``, built from the model
        at index ``model_index`` in ``models``.

        Each chain of the new model records the operators that placed it in
        ``Chain.operators``.

        Raises
        ------
        IndexError
            If there is no model at ``model_index``.
        AssemblyError
            If the structure has no assembly ``assembly_id``, or it cannot be built.

        Examples
        --------
        >>> structure = mmcif_io.read("1crn.bcif")
        >>> assembled = structure.assemble("1")
        >>> assembled.chains[0].operators
        ['1']
        """
    def assembly_is_deposited(self, assembly_id: str, tolerance: float = 1e-3) -> bool:
        """Return ``True`` if assembly ``assembly_id`` shows the deposited chains in
        their deposited arrangement.

        That is, it places no chain more than once, and places every chain by the
        same transform, to within ``tolerance`` in each rotation and translation
        element. Such an assembly may hold fewer chains than the deposited model, or
        move all of them together, but it shows nothing the deposited model does
        not.

        Raises
        ------
        AssemblyError
            If the structure has no assembly ``assembly_id``, or it names an operator
            the structure does not define.

        Examples
        --------
        >>> structure = mmcif_io.read("example.cif")
        >>> structure.assembly_is_deposited("1")
        True
        >>> structure.assembly_is_deposited("2")
        False
        """
    def to_auth(self, label_asym: str, label_seq: int) -> _AuthAddress:
        """Return the author address ``(chain, seq, ins_code)`` of the polymer residue
        at label address ``(label_asym, label_seq)``.

        ``ins_code`` is ``None`` when the residue has no insertion code.

        Raises
        ------
        NumberingError
            If ``label_seq`` is 0, or no polymer residue has that label address.

        Examples
        --------
        >>> structure = mmcif_io.read("example.cif")
        >>> structure.to_auth("AA", 3)
        ('H1', 100, 'A')
        """
    def to_label(self, chain: str, seq: int, ins_code: str | None = None) -> tuple[str, int]:
        """Return the label address ``(label_asym, label_seq)`` of the polymer residue
        at author address ``(chain, seq, ins_code)``.

        Raises
        ------
        NumberingError
            If no polymer residue has that author address.

        Examples
        --------
        >>> structure = mmcif_io.read("example.cif")
        >>> structure.to_label("H1", 100, "A")
        ('AA', 3)
        """

@final
class Entity:
    """An entity of a structure, from ``_entity``.

    Examples
    --------
    >>> entity = mmcif_io.read("1crn.bcif").entities[0]
    >>> entity.kind, entity.polymer_type, entity.description
    ('polymer', 'polypeptide(L)', 'CRAMBIN')
    >>> entity.sequence[:3]
    [['THR'], ['THR'], ['CYS']]
    """
    @property
    def id(self) -> str:
        """The ``_entity.id``."""
    @property
    def kind(self) -> str:
        """The entity type: ``"polymer"``, ``"non-polymer"``, ``"branched"``,
        ``"water"`` or ``"macrolide"``.
        """
    @property
    def description(self) -> str | None:
        """The ``_entity.pdbx_description``, or ``None``."""
    @property
    def polymer_type(self) -> str | None:
        """The ``_entity_poly.type``, such as ``"polypeptide(L)"``, or ``None`` if
        the entity is not a polymer.
        """
    @property
    def canonical_sequence(self) -> str | None:
        """The canonical one-letter sequence, from
        ``_entity_poly.pdbx_seq_one_letter_code_can``, or ``None``.
        """
    @property
    def sequence(self) -> list[list[str]]:
        """The component ids at each position of the polymer sequence. A position
        holds more than one id where the sequence is heterogeneous. Empty if
        the entity is not a polymer.
        """
    @property
    def asym_ids(self) -> list[str]:
        """The label asym ids of the chains that belong to the entity."""

@final
class Assembly:
    """A biological assembly, from ``_pdbx_struct_assembly`` and
    ``_pdbx_struct_assembly_gen``.

    Build its coordinates with ``Structure.assemble``.

    Examples
    --------
    >>> assembly = mmcif_io.read("1crn.bcif").assemblies[0]
    >>> assembly.id, assembly.details
    ('1', 'author_defined_assembly')
    >>> assembly.oligomeric_details, assembly.oligomeric_count
    ('monomeric', 1)
    >>> assembly.generators
    [('1', ['A'])]
    """
    @property
    def id(self) -> str:
        """The ``_pdbx_struct_assembly.id``."""
    @property
    def details(self) -> str | None:
        """The ``_pdbx_struct_assembly.details``, or ``None``."""
    @property
    def oligomeric_details(self) -> str | None:
        """The oligomeric state from ``_pdbx_struct_assembly.oligomeric_details``,
        such as ``"dimeric"``, or ``None``.
        """
    @property
    def oligomeric_count(self) -> int | None:
        """The number of polymer chains in the assembly from
        ``_pdbx_struct_assembly.oligomeric_count``, or ``None``.
        """
    @property
    def operator_count(self) -> int:
        """The number of operator combinations across all generators."""
    @property
    def generators(self) -> list[tuple[str, list[str]]]:
        """The ``(operator expression, label asym ids)`` pairs from
        ``_pdbx_struct_assembly_gen``, in file order.
        """

@final
class Model:
    """One model of a structure: its chains, residues and atoms.

    Per-atom data is available as lists and NumPy arrays in atom order.
    ``len(model)`` is the number of chains, and two ``Model`` objects are equal
    when they refer to the same model.

    Examples
    --------
    >>> model = mmcif_io.read("1crn.bcif").models[0]
    >>> model.number, model.chain_count, model.residue_count, model.atom_count
    (1, 1, 46, 327)
    >>> model.coordinates.shape
    (327, 3)
    """
    @property
    def number(self) -> int:
        """The model number, from ``_atom_site.pdbx_PDB_model_num``."""
    @property
    def atom_count(self) -> int:
        """The number of atoms in the model."""
    @property
    def residue_count(self) -> int:
        """The number of residues in the model."""
    @property
    def chain_count(self) -> int:
        """The number of chains in the model."""
    @property
    def chains(self) -> list[Chain]:
        """The chains of the model, in file order."""
    def chain(self, label_asym: str) -> Chain | None:
        """Return the chain with label asym id ``label_asym`` that carries no assembly
        operators, or ``None``.

        Examples
        --------
        >>> model = mmcif_io.read("1crn.bcif").models[0]
        >>> model.chain("A")
        <Chain A (auth A) with 46 residues>
        >>> model.chain("Z") is None
        True
        """
    def residues(self) -> list[Residue]:
        """Return the residues of the model, in file order.

        Examples
        --------
        >>> model = mmcif_io.read("1crn.bcif").models[0]
        >>> residues = model.residues()
        >>> len(residues), residues[0].comp
        (46, 'THR')
        """
    def atoms(self) -> list[Atom]:
        """Return the atoms of the model, in file order.

        Examples
        --------
        >>> model = mmcif_io.read("1crn.bcif").models[0]
        >>> model.atoms()[0]
        <Atom N N at [17.047, 14.099, 3.625]>
        """
    @property
    def coordinates(self) -> NDArray[np.float64]:
        """The atom positions, a float64 NumPy array of shape ``(atom_count, 3)``."""
    @property
    def elements(self) -> list[str]:
        """The element symbol of each atom."""
    @property
    def atom_names(self) -> list[str]:
        """The name of each atom."""
    @property
    def altlocs(self) -> list[str | None]:
        """The alternate location of each atom, or ``None`` for an atom without one."""
    @property
    def occupancies(self) -> NDArray[np.float64]:
        """The occupancy of each atom, a float64 NumPy array.

        An atom whose file gives no occupancy has 1.0.
        """
    @property
    def b_factors(self) -> NDArray[np.float64]:
        """The isotropic B-factor of each atom, a float64 NumPy array.

        An atom whose file gives no B-factor has NaN.
        """
    @property
    def residue_indices(self) -> NDArray[np.uint32]:
        """The index in ``residues()`` of the residue each atom belongs to, a uint32
        NumPy array.
        """
    @property
    def chain_indices(self) -> NDArray[np.uint32]:
        """The index in ``chains`` of the chain each atom belongs to, a uint32 NumPy
        array.
        """
    def single_conformer(self, policy: str = "highest_occupancy") -> Model:
        """Return a copy of the model that keeps one alternate location at each residue
        position.

        ``policy`` chooses the alternate location at each position:

        * ``"highest_occupancy"``: the one with the highest mean atom occupancy, taking
          the first on a tie.
        * ``"first"``: the one that comes first among the atoms of the position.
        * a single character: that alternate location where the position has it, and
          the first one otherwise.

        Atoms with no alternate location are always kept, and residues left with no
        atoms are removed.

        Raises
        ------
        ValueError
            If ``policy`` is none of these.

        Examples
        --------
        >>> model = mmcif_io.read("example.cif").models[0]
        >>> model.atom_count
        16
        >>> model.single_conformer().atom_count
        13
        """
    def __len__(self) -> int: ...
    def __eq__(self, other: object, /) -> bool: ...
    def __hash__(self) -> int: ...

@final
class Chain:
    """A chain of a model.

    ``len(chain)`` is the number of residues, and two ``Chain`` objects are equal
    when they refer to the same chain of the same model.

    Examples
    --------
    >>> chain = mmcif_io.read("1crn.bcif").models[0].chains[0]
    >>> chain.label_asym, chain.auth_asym, chain.entity_id
    ('A', 'A', '1')
    >>> len(chain)
    46
    """
    @property
    def index(self) -> int:
        """The position of the chain in ``Model.chains``."""
    @property
    def model(self) -> Model:
        """The model the chain belongs to."""
    @property
    def label_asym(self) -> str:
        """The label asym id, from ``_atom_site.label_asym_id``."""
    @property
    def auth_asym(self) -> str:
        """The author chain id, from ``_atom_site.auth_asym_id``."""
    @property
    def entity_id(self) -> str | None:
        """The id of the entity the chain belongs to, or ``None``."""
    @property
    def operators(self) -> list[str]:
        """The ids of the assembly operators that placed the chain, in the order they
        were applied. Empty for a chain read from a file.
        """
    @property
    def instance_id(self) -> str:
        """The id of the chain's copy, named as in the wwPDB assembly files.

        A chain read from a file, or placed by an identity transform, has its
        author chain id. Any other copy has the author chain id followed by a hyphen
        and its operator ids, joined by hyphens, as in ``A-2`` or ``A-X0-1``.

        Examples
        --------
        >>> structure = mmcif_io.read("example.cif")
        >>> [chain.instance_id for chain in structure.assemble("2").chains]
        ['H1', 'H1', 'H1-2', 'H1-2']
        """
    @property
    def residues(self) -> list[Residue]:
        """The residues of the chain, in file order."""
    @property
    def atom_range(self) -> tuple[int, int]:
        """The ``(start, end)`` range of the chain's atoms in the model's atom order,
        with ``end`` exclusive.
        """
    @property
    def coordinates(self) -> NDArray[np.float64]:
        """The positions of the chain's atoms, a float64 NumPy array of shape
        ``(n, 3)``.
        """
    def residues_at_label(self, seq: int) -> list[Residue]:
        """Return the residues at label sequence number ``seq``, in file order.

        The list holds more than one residue where the chain is heterogeneous at that
        position, and is empty if no residue has the number.

        Examples
        --------
        >>> chain = mmcif_io.read("example.cif").models[0].chain("AA")
        >>> chain.residues_at_label(2)
        [<Residue AA/2/100/SER>, <Residue AA/2/100/THR>]
        """
    def residues_at_auth(self, seq: int, ins_code: str | None = None) -> list[Residue]:
        """Return the residues at author residue number ``seq`` with insertion code
        ``ins_code``, in file order.

        Examples
        --------
        >>> chain = mmcif_io.read("example.cif").models[0].chain("AA")
        >>> chain.residues_at_auth(100, "A")
        [<Residue AA/3/100A/ALA>]
        """
    def __len__(self) -> int: ...
    def __eq__(self, other: object, /) -> bool: ...
    def __hash__(self) -> int: ...

@final
class Residue:
    """A residue of a chain.

    ``len(residue)`` is the number of atoms, and two ``Residue`` objects are equal
    when they refer to the same residue of the same model.

    Examples
    --------
    >>> residue = mmcif_io.read("1crn.bcif").models[0].chains[0].residues[2]
    >>> residue.comp, residue.label_seq, residue.auth_seq
    ('CYS', 3, 3)
    >>> residue.atom("SG").position
    (15.885, 12.426, 10.016)
    """
    @property
    def index(self) -> int:
        """The position of the residue in ``Model.residues()``."""
    @property
    def chain(self) -> Chain:
        """The chain the residue belongs to."""
    @property
    def comp(self) -> str:
        """The chemical component id, from ``_atom_site.label_comp_id``."""
    @property
    def label_seq(self) -> int | None:
        """The label sequence number, from ``_atom_site.label_seq_id``, or ``None`` if
        the residue has none.
        """
    @property
    def auth_seq(self) -> int:
        """The author residue number, from ``_atom_site.auth_seq_id``."""
    @property
    def ins_code(self) -> str | None:
        """The insertion code, from ``_atom_site.pdbx_PDB_ins_code``, or ``None``."""
    @property
    def is_hetero(self) -> bool:
        """``True`` if the residue's first atom was read from a ``HETATM`` record."""
    @property
    def atoms(self) -> list[Atom]:
        """The atoms of the residue, in file order."""
    def atom(self, name: str) -> Atom | None:
        """Return the first atom of the residue named ``name``, or ``None``.

        Examples
        --------
        >>> residue = mmcif_io.read("1crn.bcif").models[0].chains[0].residues[2]
        >>> residue.atom("CA").name
        'CA'
        >>> residue.atom("OXT") is None
        True
        """
    @property
    def altlocs(self) -> list[str]:
        """The distinct alternate locations of the residue's atoms, in order of first
        appearance.
        """
    @property
    def atom_range(self) -> tuple[int, int]:
        """The ``(start, end)`` range of the residue's atoms in the model's atom
        order, with ``end`` exclusive.
        """
    @property
    def coordinates(self) -> NDArray[np.float64]:
        """The positions of the residue's atoms, a float64 NumPy array of shape
        ``(n, 3)``.
        """
    @property
    def key(self) -> _ResidueKey:
        """The residue key ``(label_asym, label_seq, auth_seq, ins_code, comp,
        operators)``, where ``operators`` is the tuple of assembly operators that
        placed the residue's chain, empty in a deposited model.
        """
    @property
    def auth_key(self) -> _AuthAddress:
        """The author address ``(auth_asym, auth_seq, ins_code)``."""
    def __len__(self) -> int: ...
    def __eq__(self, other: object, /) -> bool: ...
    def __hash__(self) -> int: ...

@final
class Atom:
    """An atom of a residue.

    Two ``Atom`` objects are equal when they refer to the same atom of the same
    model.

    Examples
    --------
    >>> atom = mmcif_io.read("1crn.bcif").models[0].atoms()[1]
    >>> atom.name, atom.element, atom.serial
    ('CA', 'C', 2)
    >>> atom.position
    (16.967, 12.784, 4.338)
    """
    @property
    def index(self) -> int:
        """The position of the atom in the model's atom order."""
    @property
    def residue(self) -> Residue:
        """The residue the atom belongs to."""
    @property
    def serial(self) -> int | None:
        """The atom serial number, from ``_atom_site.id``, or ``None``."""
    @property
    def name(self) -> str:
        """The atom name, from ``_atom_site.label_atom_id``."""
    @property
    def element(self) -> str:
        """The element symbol, from ``_atom_site.type_symbol``."""
    @property
    def altloc(self) -> str | None:
        """The alternate location, from ``_atom_site.label_alt_id``, or ``None``."""
    @property
    def position(self) -> tuple[float, float, float]:
        """The ``(x, y, z)`` position."""
    @property
    def occupancy(self) -> float:
        """The occupancy, or 1.0 if the file gives none."""
    @property
    def b_factor(self) -> float | None:
        """The isotropic B-factor, from ``_atom_site.B_iso_or_equiv``, or ``None`` if
        the file gives none.
        """
    @property
    def charge(self) -> int | None:
        """The formal charge, from ``_atom_site.pdbx_formal_charge``, or ``None``."""
    def __eq__(self, other: object, /) -> bool: ...
    def __hash__(self) -> int: ...

@final
class NeighbourGrid:
    """A spatial index over points for distance queries.

    ``NeighbourGrid(points, cell_size=4.0)`` builds a grid over ``points``, an
    ``(n, 3)`` float64 NumPy array, with cubic cells of edge length ``cell_size``.
    The grid keeps its own copy of the points, and queries return indices into
    them. ``len(grid)`` is the number of points.

    Raises
    ------
    ValueError
        If ``points`` is not an ``(n, 3)`` array, a point has a NaN or infinite
        coordinate, or ``cell_size`` is not positive and finite.

    Examples
    --------
    >>> import numpy as np
    >>> grid = mmcif_io.NeighbourGrid(np.array([[0.0, 0, 0], [1, 0, 0], [5, 0, 0]]))
    >>> grid
    <NeighbourGrid with 3 points, cell size 4>
    >>> grid.within([0, 0, 0], 1.0).tolist()
    [0, 1]
    """
    def __new__(cls, points: NDArray[np.float64], cell_size: float = 4.0) -> NeighbourGrid: ...
    @staticmethod
    def for_model(model: Model, cell_size: float = 4.0) -> NeighbourGrid:
        """Build a grid over the atom positions of ``model``, with cubic cells of edge
        length ``cell_size``.

        Queries return atom indices, as in ``Model.atoms()``.

        Raises
        ------
        ValueError
            If ``cell_size`` is not positive and finite.

        Examples
        --------
        >>> model = mmcif_io.read("1crn.bcif").models[0]
        >>> mmcif_io.NeighbourGrid.for_model(model)
        <NeighbourGrid with 327 points, cell size 4>
        """
    @property
    def cell_size(self) -> float:
        """The edge length of a grid cell."""
    def within(self, point: _Point, radius: float) -> NDArray[np.uint32]:
        """Return the indices of the points within ``radius`` of ``point``, in ascending
        order, as a uint32 NumPy array.

        Raises
        ------
        ValueError
            If ``radius`` is negative or NaN, or ``point`` has a NaN or infinite
            coordinate.

        Examples
        --------
        >>> model = mmcif_io.read("1crn.bcif").models[0]
        >>> grid = mmcif_io.NeighbourGrid.for_model(model)
        >>> atoms = model.atoms()
        >>> [atoms[index].name for index in grid.within(atoms[1].position, 1.6)]
        ['N', 'CA', 'C', 'CB']
        """
    def nearest(self, point: _Point, max_radius: float = ...) -> tuple[int, float] | None:
        """Return ``(index, distance)`` for the point closest to ``point``, or ``None`` if
        no point lies within ``max_radius``.

        When several points are equally close, the one with the lowest index is
        returned.

        Raises
        ------
        ValueError
            If ``max_radius`` is negative or NaN, or ``point`` has a NaN or infinite
            coordinate.

        Examples
        --------
        >>> model = mmcif_io.read("1crn.bcif").models[0]
        >>> grid = mmcif_io.NeighbourGrid.for_model(model)
        >>> grid.nearest(model.atoms()[1].position)
        (1, 0.0)
        >>> grid.nearest([0, 0, 0], max_radius=1.0) is None
        True
        """
    def pairs_within(self, radius: float) -> NDArray[np.uint32]:
        """Return every pair of points within ``radius`` of each other, as a uint32 NumPy
        array of shape ``(n, 2)``.

        Each row is ``(lower, higher)`` by index, and the rows are in ascending order.

        Raises
        ------
        ValueError
            If ``radius`` is negative or NaN.

        Examples
        --------
        >>> model = mmcif_io.read("1crn.bcif").models[0]
        >>> grid = mmcif_io.NeighbourGrid.for_model(model)
        >>> grid.pairs_within(1.6)[:3].tolist()
        [[0, 1], [1, 2], [1, 4]]
        """
    def contacts(self, others: NDArray[np.float64], radius: float) -> NDArray[np.uint32]:
        """Return an ``(other, point)`` index pair for each point in ``others`` that lies
        within ``radius`` of a point in the grid, as a uint32 NumPy array of shape
        ``(n, 2)``.

        ``others`` is an ``(m, 3)`` float64 NumPy array. Rows are ordered by ``other``
        and then by ``point``.

        Raises
        ------
        ValueError
            If ``others`` is not an ``(m, 3)`` array, a point in ``others`` has a NaN
            or infinite coordinate, or ``radius`` is negative or NaN.

        Examples
        --------
        >>> import numpy as np
        >>> grid = mmcif_io.NeighbourGrid(np.array([[0.0, 0, 0], [1, 0, 0], [5, 0, 0]]))
        >>> grid.contacts(np.array([[4.5, 0, 0], [0.5, 0, 0]]), 1.0).tolist()
        [[0, 2], [1, 0], [1, 1]]
        """
    def __len__(self) -> int: ...

@final
class Document:
    """A CIF document: data blocks of categories, each a table of named items.

    ``len(document)`` is the number of data blocks.

    Examples
    --------
    >>> document = mmcif_io.Document.open("1crn.cif")
    >>> document
    <Document with 1 blocks>
    >>> block = document.blocks[0]
    >>> block["atom_site"].row_count
    327
    """
    @staticmethod
    def parse(text: str) -> Document:
        r"""Parse a CIF document from text.

        Raises
        ------
        ReadError
            If ``text`` is not valid CIF.

        Examples
        --------
        >>> document = mmcif_io.Document.parse("data_x\n_cell.length_a 10.0\n")
        >>> document.blocks[0]["cell"]["length_a"]
        ['10.0']
        """
    @staticmethod
    def open(path: _Path) -> Document:
        """Read the CIF document in the file at ``path``.

        Raises
        ------
        OSError
            If the file cannot be read.
        ReadError
            If the contents are not valid CIF.

        Examples
        --------
        >>> mmcif_io.Document.open("1crn.cif").blocks[0].name
        '1CRN'
        """
    @staticmethod
    def from_bcif(data: bytes) -> Document:
        """Decode a CIF document from BinaryCIF bytes.

        Raises
        ------
        ReadError
            If ``data`` is not valid BinaryCIF.

        Examples
        --------
        >>> with open("1crn.bcif", "rb") as file:
        ...     document = mmcif_io.Document.from_bcif(file.read())
        >>> document.blocks[0].name
        '1CRN'
        """
    @staticmethod
    def open_bcif(path: _Path) -> Document:
        """Read the BinaryCIF document in the file at ``path``.

        Raises
        ------
        OSError
            If the file cannot be read.
        ReadError
            If the contents are not valid BinaryCIF.

        Examples
        --------
        >>> mmcif_io.Document.open_bcif("1crn.bcif").blocks[0].name
        '1CRN'
        """
    @property
    def blocks(self) -> list[Block]:
        """The data blocks of the document, in file order."""
    def block(self, name: str) -> Block | None:
        """Return the data block named ``name``, or ``None``.

        Names are compared without regard to ASCII case.

        Examples
        --------
        >>> document = mmcif_io.Document.open("1crn.cif")
        >>> document.block("1crn").name
        '1CRN'
        >>> document.block("2abc") is None
        True
        """
    def to_cif(self) -> str:
        r"""Return the document written as CIF text.

        Parsing the text gives a document with the same blocks, categories and values.

        Examples
        --------
        >>> document = mmcif_io.Document.parse("data_x\n_cell.length_a 10.0\n")
        >>> print(document.to_cif())
        data_x
        #
        _cell.length_a 10.0
        #
        <BLANKLINE>
        """
    def __len__(self) -> int: ...

@final
class Block:
    """A data block of a CIF document.

    ``block[name]`` returns the category ``name`` and raises ``KeyError`` if there
    is none, and ``name in block`` tests for one. Category names may be given
    with or without the leading underscore.

    Examples
    --------
    >>> block = mmcif_io.Document.open("1crn.cif").blocks[0]
    >>> "atom_site" in block
    True
    >>> block["_entry"]["id"]
    ['1CRN']
    """
    @property
    def name(self) -> str:
        """The name of the block, without the ``data_`` prefix."""
    @property
    def category_names(self) -> list[str]:
        """The names of the categories in the block, in file order, without the leading
        underscore.
        """
    def category(self, name: str) -> Category | None:
        """Return the category ``name``, or ``None``.

        ``name`` may be given with or without the leading underscore.

        Examples
        --------
        >>> block = mmcif_io.Document.open("1crn.cif").blocks[0]
        >>> block.category("cell").name
        'cell'
        >>> block.category("nothing") is None
        True
        """
    def __getitem__(self, key: str, /) -> Category: ...
    def __contains__(self, key: str, /) -> bool: ...

@final
class Category:
    """A category of a data block: a table whose columns are items.

    ``category[item]`` returns the column ``item`` as ``column`` does and raises
    ``KeyError`` if there is none, ``item in category`` tests for one, and
    ``len(category)`` is the number of rows.

    Examples
    --------
    >>> category = mmcif_io.Document.open("1crn.cif").blocks[0]["atom_site"]
    >>> category
    <Category atom_site with 327 rows>
    >>> category["label_comp_id"][:3]
    ['THR', 'THR', 'THR']
    """
    @property
    def name(self) -> str:
        """The name of the category, without the leading underscore."""
    @property
    def items(self) -> list[str]:
        """The names of the items in the category, in file order, without the category
        prefix.
        """
    @property
    def row_count(self) -> int:
        """The number of rows in the category."""
    def column(self, item: str) -> list[str | None] | None:
        """Return the values of ``item`` in row order, or ``None`` if the category has no
        such item.

        Values that are ``?`` or ``.`` in the file are returned as ``None``.

        Examples
        --------
        >>> category = mmcif_io.Document.open("1crn.cif").blocks[0]["atom_site"]
        >>> category.column("label_seq_id")[:3]
        ['1', '1', '1']
        >>> category.column("pdbx_PDB_ins_code")[:3]
        [None, None, None]
        """
    def raw_column(self, item: str) -> list[str] | None:
        """Return the values of ``item`` in row order as written in the file, or
        ``None`` if the category has no such item.

        ``?`` and ``.`` are returned as the strings ``"?"`` and ``"."``.

        Examples
        --------
        >>> category = mmcif_io.Document.open("1crn.cif").blocks[0]["atom_site"]
        >>> category.raw_column("pdbx_PDB_ins_code")[:3]
        ['?', '?', '?']
        """
    def __getitem__(self, key: str, /) -> list[str | None]: ...
    def __contains__(self, key: str, /) -> bool: ...
    def __len__(self) -> int: ...

@final
class CcdIndex:
    """A compiled Chemical Component Dictionary index that returns components by id.

    The index records the name of the dictionary release it was built from.
    ``index[id]`` returns the component ``id`` as ``get`` does and raises
    ``KeyError`` if there is none, ``id in index`` tests for one, and
    ``len(index)`` is the number of components.

    Examples
    --------
    >>> index = mmcif_io.CcdIndex.build("components.cif", "2026-10-20")
    >>> index
    <CcdIndex 2026-10-20 with 8 components>
    >>> index["ALA"].name
    'ALANINE'
    >>> "XYZ" in index
    False
    """
    @staticmethod
    def open(path: _Path) -> CcdIndex:
        """Read the index file at ``path``, as written by ``build``.

        Raises
        ------
        OSError
            If the file cannot be read.
        ReadError
            If the contents are not a component index.

        Examples
        --------
        >>> mmcif_io.CcdIndex.build("components.cif", "2026-10-20", out="ccd.idx")
        <CcdIndex 2026-10-20 with 8 components>
        >>> mmcif_io.CcdIndex.open("ccd.idx")
        <CcdIndex 2026-10-20 with 8 components>
        """
    @staticmethod
    def build(source: _Path, release: str, out: _Path | None = None) -> CcdIndex:
        """Compile the Chemical Component Dictionary file at ``source`` into an index for
        the release named ``release``.

        ``source`` is CIF with one data block per component, as in ``components.cif``.
        If ``out`` is given, the index is also written to that path for ``open`` to
        read.

        Raises
        ------
        OSError
            If ``source`` cannot be read or is not valid UTF-8, or ``out`` cannot be
            written.
        ReadError
            If ``source`` is not valid CIF, a block cannot be read as a component, a
            component id is not 1 to 8 ASCII characters without whitespace, or two
            blocks have the same component id.

        Examples
        --------
        >>> index = mmcif_io.CcdIndex.build("components.cif", "2026-10-20")
        >>> index.release
        '2026-10-20'
        """
    @staticmethod
    def from_bytes(data: bytes) -> CcdIndex:
        """Read an index from ``data``, the contents of a file written by ``build``.

        Raises
        ------
        ReadError
            If ``data`` is not a component index or was written in a format version
            this release does not read.

        Examples
        --------
        >>> mmcif_io.CcdIndex.build("components.cif", "2026-10-20", out="ccd.idx")
        <CcdIndex 2026-10-20 with 8 components>
        >>> with open("ccd.idx", "rb") as file:
        ...     index = mmcif_io.CcdIndex.from_bytes(file.read())
        >>> len(index)
        8
        """
    @property
    def release(self) -> str:
        """The name of the dictionary release given when the index was built."""
    def ids(self) -> list[str]:
        """Return the component ids in the index, in ascending byte order.

        Examples
        --------
        >>> index = mmcif_io.CcdIndex.build("components.cif", "2026-10-20")
        >>> index.ids()[:3]
        ['A1AAA', 'ALA', 'ATP']
        """
    def get(self, id: str) -> Component | None:
        """Return the component with id ``id``, or ``None`` if the index does not hold
        it.

        Ids are matched exactly, including case.

        Raises
        ------
        ReadError
            If the index record for ``id`` is corrupt.

        Examples
        --------
        >>> index = mmcif_io.CcdIndex.build("components.cif", "2026-10-20")
        >>> index.get("MSE").parent
        'MET'
        >>> index.get("mse") is None
        True
        """
    def __getitem__(self, key: str, /) -> Component: ...
    def __contains__(self, key: str, /) -> bool: ...
    def __len__(self) -> int: ...

@final
class Component:
    """A chemical component from the dictionary, with its atoms and the bonds between
    them.

    Examples
    --------
    >>> component = mmcif_io.CcdIndex.build("components.cif", "2026-10-20")["ALA"]
    >>> component
    <Component ALA ALANINE>
    >>> component.kind, component.formula
    ('L-PEPTIDE LINKING', 'C3 H7 N O2')
    >>> component.atoms[1]
    ('CA', 'C', 0, False, False, 'S')
    >>> component.bonds[0]
    (0, 1, 'SING', False, 'N')
    """
    @property
    def id(self) -> str:
        """The component id, such as ``ALA``."""
    @property
    def name(self) -> str:
        """The component name, such as ``ALANINE``, or an empty string if the
        dictionary gives none.
        """
    @property
    def kind(self) -> str:
        """The component type, such as ``L-PEPTIDE LINKING``, or an empty string if the
        dictionary gives none.
        """
    @property
    def one_letter_code(self) -> str | None:
        """The one-letter code, such as ``A`` for alanine, or ``None``."""
    @property
    def parent(self) -> str | None:
        """The id of the parent component, such as ``MET`` for selenomethionine, or
        ``None``.
        """
    @property
    def formula(self) -> str | None:
        """The chemical formula, such as ``C3 H7 N O2``, or ``None``."""
    @property
    def formal_charge(self) -> int:
        """The formal charge of the component, or 0 if the dictionary gives none."""
    @property
    def release_status(self) -> str:
        """The release status, such as ``REL``, or an empty string if the dictionary
        gives none.
        """
    @property
    def atoms(self) -> list[_AtomRow]:
        """The atoms, in the order the dictionary lists them, as
        ``(name, element, charge, is_leaving, is_aromatic, stereo)`` tuples.

        ``element`` is as the dictionary writes it, such as ``C`` or ``FE``.
        ``charge`` is 0 if the dictionary gives none, and ``stereo`` is ``"R"``,
        ``"S"`` or ``"N"`` for none.
        """
    @property
    def bonds(self) -> list[_BondRow]:
        """The bonds, in the order the dictionary lists them, as
        ``(first, second, order, is_aromatic, stereo)`` tuples.

        ``first`` and ``second`` are positions in ``atoms``. ``order`` is the
        ``_chem_comp_bond.value_order`` code: ``SING``, ``DOUB``, ``TRIP``, ``QUAD``,
        ``AROM``, ``DELO``, ``PI`` or ``POLY``. ``stereo`` is ``"E"``, ``"Z"`` or
        ``"N"`` for none.
        """
