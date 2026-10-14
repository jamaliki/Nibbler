"""Pack observed coordinates in native code under a caller-owned slot schema.

ObservationPacker accepts a JSON schema and retains its ligand feature tables.
add_hydrogen_parents stores dictionary hydrogen-to-heavy-atom names per component.
Ligand packing uses observed covalent bonds first, then dictionary parents, and
uses distance only when neither source names a parent. Protein packing uses the
slot schema's named parent before distance-based assignment.
afdb_eligibility scans a read-only A14P buffer without decoding a coordinate
array or retaining the source. It counts complete residues within the same
float32 peptide-chain rule as protein_document.
Its pack method returns (arrays, keys). Each array is a (dtype, shape, buffer)
tuple. The buffers own their storage and support writable NumPy views. Native
packing releases the GIL and never serializes or reparses the source document.
"""

from ._core import ObservationPacker

__all__ = ["ObservationPacker"]
