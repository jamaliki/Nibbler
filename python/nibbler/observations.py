"""Pack observed coordinates in native code under a caller-owned slot schema.

ObservationPacker accepts a JSON schema and retains its ligand feature tables.
Its pack method returns (arrays, keys). Each array is a (dtype, shape, buffer)
tuple. The buffers own their storage and support writable NumPy views. Native
packing releases the GIL and never serializes or reparses the source document.
"""

from ._core import ObservationPacker

__all__ = ["ObservationPacker"]
