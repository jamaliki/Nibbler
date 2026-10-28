"""Ligand token tables and output digests for packer tests."""

import gzip
import hashlib
import json
import os
from contextlib import contextmanager
from pathlib import Path

import nibbler


def token_table(*elements, start=21, other=64, hydrogen=65):
    """Give the elements consecutive tokens from ``start``, as the 66-token contract does."""
    return dict(
        elements={element: start + index for index, element in enumerate(elements)},
        other=other,
        hydrogen=hydrogen,
    )


def digest(value):
    """Packer output with each buffer, and each long list, replaced by its SHA-256 digest."""
    if isinstance(value, dict):
        return {str(name): digest(item) for name, item in value.items()}
    if (
        isinstance(value, tuple)
        and len(value) == 3
        and isinstance(value[2], (bytes, bytearray, memoryview))
    ):
        dtype, shape, buffer = value
        return dict(
            dtype=dtype,
            shape=list(shape),
            sha256=hashlib.sha256(bytes(buffer)).hexdigest(),
        )
    if isinstance(value, (bytes, bytearray, memoryview)):
        return dict(sha256=hashlib.sha256(bytes(value)).hexdigest())
    if isinstance(value, (list, tuple)):
        items = [digest(item) for item in value]
        if len(items) <= 8:
            return items
        encoded = json.dumps(items, sort_keys=True, separators=(",", ":")).encode()
        return dict(items=len(items), sha256=hashlib.sha256(encoded).hexdigest())
    return value


FIXTURES = Path(__file__).parent / "fixtures" / "observations"


def frozen(name):
    """A frozen kaveh-data contract or the recorded digests, from the observation fixtures."""
    with gzip.open(FIXTURES / name, "rt") as file:
        return file.read()


def with_entry(raw):
    """The document with an ``_entry.id``, which the molecular view requires, after its data line."""
    if b"\n_entry.id" in raw:
        return raw
    start = raw.index(b"data_")
    end = raw.index(b"\n", start) + 1
    return raw[:end] + b"_entry.id PARITY\n" + raw[end:]


# Where Reduce3 looks for chem_data when it is given none.
CHEM_DATA_VARIABLES = (
    "REDUCE3_CHEM_DATA",
    "CHEM_DATA",
    "MMTBX_CCP4_MONOMER_LIB",
    "CLIBD_MON",
    "CONDA_PREFIX",
)


@contextmanager
def without_chem_data():
    saved = {
        name: os.environ.pop(name) for name in CHEM_DATA_VARIABLES if name in os.environ
    }
    try:
        yield
    finally:
        os.environ.update(saved)


def packed(packer, raw, chem_data):
    """Digests of the legacy and molecular packs of one document, or the error each raised."""
    document = nibbler.cif.read(raw)

    def outcome(run):
        try:
            return digest(run())
        except Exception as error:
            return dict(error=type(error).__name__, message=str(error))

    def molecular():
        view = packer.molecular_view(
            nibbler.cif.read(with_entry(raw)), chem_data=chem_data
        )
        perceived = [
            bytes(32 * len(graph["atoms"])) for graph in view.perception_graphs()
        ]
        return packer.finish_molecular(view, perceived, complete=True)

    # Given no chem_data, the packer finds one from the environment; these two run without it.
    with without_chem_data():
        result = dict(
            plain=outcome(lambda: packer.pack(document)),
            complete=outcome(
                lambda: packer.pack(
                    document, complete=True, complete_tokens=True, validate=True
                )
            ),
        )
    if chem_data is not None:
        result["chem_data"] = outcome(
            lambda: packer.pack(
                document, complete=True, validate=True, chem_data=chem_data
            )
        )
        result["molecular"] = outcome(molecular)
    return result
