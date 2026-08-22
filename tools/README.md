# Build and qualification workflows

These tools are side-effect boundaries around the offline library. Network access is
limited to explicit schema/corpus fetch commands; runtime parsing and validation never
fetch data.

The PDB corpus fetcher downloads immutable, revision-addressed text/gzip entries from
the wwPDB versioned archive and exact hash-pinned BinaryCIF files from RCSB. It writes
through temporary files and only promotes bytes that match the manifest.

## Lock-pinned schema supply chain

```mermaid
flowchart TB
    LOCK["schemas/locks.toml<br/>URL + version + SHA-256"]:::input --> FETCH["fetch_schemas.py"]:::hot
    FETCH --> CACHE{"Cached file has<br/>locked digest?"}:::proof
    CACHE -- no --> DOWNLOAD["Download to temp<br/>verify SHA-256"]:::hot
    DOWNLOAD --> ATOMIC["Atomic cache replace"]:::proof
    CACHE -- yes --> SOURCE["Verified DDL2 source"]:::data
    ATOMIC --> SOURCE

    SOURCE --> PARSE["Nibbler strict CIF parser"]:::hot
    PARSE --> COMPILE["compile_dictionary<br/>types, categories, items,<br/>aliases, links, inherited rules"]:::hot
    COMPILE --> VERIFY["Validate compiled references"]:::proof
    VERIFY --> ENCODE["Deterministic NIBSCM v2 MessagePack"]:::hot
    ENCODE --> REPLACE["Atomic schemas/compiled/*.nbs replace"]:::proof
    REPLACE --> EMBED["include_bytes! release input"]:::output

    classDef input fill:#eef2ff,stroke:#4f46e5,color:#1e1b4b;
    classDef hot fill:#fff1f2,stroke:#e11d48,color:#881337,stroke-width:2px;
    classDef proof fill:#ecfdf5,stroke:#059669,color:#064e3b;
    classDef data fill:#ecfeff,stroke:#0891b2,color:#164e63;
    classDef output fill:#f5f3ff,stroke:#7c3aed,color:#4c1d95,stroke-width:2px;
```

## Optional native PGO artifact

PGO is host-, compiler-, ABI-, source-, and corpus-specific. It is deliberately separate
from the portable release matrix.

```mermaid
flowchart TB
    START["Current source + lockfile<br/>pinned PDB corpus"]:::input --> CLEAN["Reject inherited profile flags<br/>require native target"]:::proof
    CLEAN --> FINGER["Hash compiler/profile inputs"]:::proof
    FINGER --> INSTR["Build instrumented release wheel"]:::hot
    INSTR --> EXTRACT["Extract to isolated PYTHONPATH"]:::data
    EXTRACT --> TRAIN["Train text, gzip, BinaryCIF<br/>projection + selected full documents"]:::hot
    TRAIN --> ASSERT["Verify corpus hashes<br/>and projected atom counts"]:::proof
    ASSERT --> MERGE["Toolchain-matched llvm-profdata merge"]:::hot
    MERGE --> STABLE{"Inputs and compiler unchanged?"}:::proof
    STABLE -- no --> REFUSE["Refuse stale profile"]:::error
    STABLE -- yes --> OPT["Rebuild with profile-use"]:::hot
    OPT --> FINAL{"Inputs and compiler<br/>still identical?"}:::proof
    FINAL -- no --> REFUSE
    FINAL -- yes --> OUT["Native wheel + .pgo.json fingerprint"]:::output

    classDef input fill:#eef2ff,stroke:#4f46e5,color:#1e1b4b;
    classDef hot fill:#fff1f2,stroke:#e11d48,color:#881337,stroke-width:2px;
    classDef proof fill:#ecfdf5,stroke:#059669,color:#064e3b;
    classDef data fill:#ecfeff,stroke:#0891b2,color:#164e63;
    classDef output fill:#f5f3ff,stroke:#7c3aed,color:#4c1d95,stroke-width:2px;
    classDef error fill:#fef2f2,stroke:#b91c1c,color:#7f1d1d;
```

## Portable release qualification

The release workflow builds 5 CPython versions on 5 platform targets plus one source
distribution. Publication credentials are unavailable until every artifact has passed
the same archive and installed-package boundary.

```mermaid
flowchart TB
    META["Validate Cargo/Python/changelog/tag metadata"]:::proof --> MATRIX["5 CPython ABIs x 5 platforms"]:::input
    META --> SDIST["Build source distribution"]:::input
    MATRIX --> WHEELS["Build locked, stripped,<br/>interpreter-specific wheels"]:::hot

    WHEELS --> QUALIFY["Fresh venv outside checkout"]:::proof
    SDIST --> REBUILD["Fresh venv<br/>rebuild sdist"]:::proof
    QUALIFY --> INSTALL["Install archive + pinned PyArrow"]:::hot
    REBUILD --> INSTALL
    INSTALL --> SMOKE["Text/gzip/BinaryCIF, schema typing,<br/>missing states, writing, Arrow"]:::proof

    SMOKE --> ASSEMBLE["Require 25 unique wheels + 1 sdist<br/>inspect contents and metadata"]:::proof
    ASSEMBLE --> HASH["Write SHA256SUMS"]:::data
    HASH --> EVENT{"Workflow event"}:::decision
    EVENT -- manual --> CANDIDATE["Retain qualified candidate<br/>no publication identity"]:::output
    EVENT -- exact version tag --> PYPI["OIDC publish complete set to PyPI"]:::output
    PYPI --> GH["Create GitHub release<br/>same archives + checksums"]:::output

    classDef input fill:#eef2ff,stroke:#4f46e5,color:#1e1b4b;
    classDef decision fill:#f8fafc,stroke:#64748b,color:#0f172a;
    classDef hot fill:#fff1f2,stroke:#e11d48,color:#881337,stroke-width:2px;
    classDef proof fill:#ecfdf5,stroke:#059669,color:#064e3b;
    classDef data fill:#ecfeff,stroke:#0891b2,color:#164e63;
    classDef output fill:#f5f3ff,stroke:#7c3aed,color:#4c1d95,stroke-width:2px;
```

`qualify_release.py` installs one archive into a temporary isolated environment and
runs `installed_smoke.py` with `python -I`, preventing the checkout from satisfying an
import accidentally.
