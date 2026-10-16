"""Pack observed coordinates in native code under a caller-owned slot schema.

ObservationPacker accepts a JSON schema and retains its ligand feature tables.
add_hydrogen_parents stores dictionary hydrogen-to-heavy-atom names per component.
Ligand packing uses observed covalent bonds first, then dictionary parents, and
uses distance only when neither source names a parent. Protein packing uses the
slot schema's named parent before distance-based assignment.
When pack receives chem_data, Rust selects each ligand residue's Reduce3
dictionary. missing_hydrogen_dictionaries lists components that need CCD parents;
it does not send dictionary maps through Python.
afdb_eligibility scans a read-only A14P buffer without decoding a coordinate
array or retaining the source. It counts complete residues within the same
float32 peptide-chain rule as protein_document.
Its pack method returns (arrays, keys). Each array is a (dtype, shape, buffer)
tuple. The buffers own their storage and support writable NumPy views. Native
packing releases the GIL and never serializes or reparses the source document.
prepare_source returns an immutable selected-model snapshot. Its dictionary and
component queries share that source with pack_prepared, which reads the packer's
current feature tables and does not accept a second model argument. Later edits
to the document do not alter the snapshot. kaveh-data creates one per pack call.

product_targets accepts a JSON list of occurrence-specific product records with
parent_key, sequence_id, formal_charge, hybridization, hydrogens, complete, and
modified fields. The packer applies these before token selection. Ligand H and
charge labels require complete product coverage and a matching observed H count;
an empty list clears both known masks. Other modified protein residues retain
observed H coordinates but lose exact state and absent-H labels. The supported
linked Cys states keep their existing rules. Omitting product_targets retains
the caller's component-level feature policy.
"""

from ._core import ObservationPacker, PreparedObservationSource

__all__ = ["ObservationPacker", "PreparedObservationSource"]
