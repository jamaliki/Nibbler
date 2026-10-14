"""Pack observed coordinates in native code under a caller-owned slot schema.

ObservationPacker accepts a JSON schema and retains its ligand feature tables.
add_hydrogen_parents stores CCD hydrogen-to-heavy-atom names once per component.
Packing uses observed covalent bonds first, then CCD parents, and uses distance
only when neither source names a parent.
Its pack method returns (arrays, keys). Each array is a (dtype, shape, buffer)
tuple. The buffers own their storage and support writable NumPy views. Native
packing releases the GIL and never serializes or reparses the source document.
"""

from ._core import ObservationPacker

__all__ = ["ObservationPacker"]
