"""Record the packer outputs that the explicit ligand token table must reproduce.

Run once with the Nibbler 9b8041e wheel, which derived ligand tokens from an element list,
and with REDUCE3_CHEM_DATA set to chem_data 2026.9:

    python tests/record_token_parity.py

That packer took the 43 ligand elements and hydrogen as kaveh-data passed them, so element
``i`` carried token ``21 + i``, other elements 64 and a bridging hydride 65. The frozen Atom32
and Atom84 schemas beside the digests are the kaveh-data contracts it read.
Only molecular crop_positions digests were subsequently corrected for the no-bridge bug;
rerunning this historical recorder restores the affected zero-coordinate digests.
"""

import gzip
import json
import os
from pathlib import Path

from nibbler.observations import ObservationPacker

from token_tables import FIXTURES, frozen, packed

ROOT = Path(__file__).parent
SKIPPED = {"fixtures/syntax"}

record = json.loads(frozen("token_table_parity.json.gz"))
packer = ObservationPacker(frozen("atom32_schema.json.gz"), [*record["ligand_elements"], "H"], record["excluded"])
packer.set_molecular_schema(frozen("molecular_schema.json.gz"))
chem_data = os.environ["REDUCE3_CHEM_DATA"]
cases = sorted(str(file.relative_to(ROOT)) for file in ROOT.rglob("*.cif")
               if not any(str(file.relative_to(ROOT)).startswith(skipped) for skipped in SKIPPED))
record["cases"] = {case: packed(packer, (ROOT / case).read_bytes(), chem_data) for case in cases}
with gzip.GzipFile(FIXTURES / "token_table_parity.json.gz", "wb", mtime=0) as file:
    file.write(json.dumps(record, indent=1, sort_keys=True).encode())
