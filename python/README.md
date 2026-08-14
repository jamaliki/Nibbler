# Python boundary algorithms

`python/nibbler/` is a typed facade over the native extension. It normalizes Python
inputs, dispatches explicit public operations, converts native errors once, and owns
filesystem transaction semantics. Parsing, projection, validation, and semantic model
construction stay in Rust.

## Single-source calls

```mermaid
flowchart TB
    API["chomp / cif.read<br/>mmcif.read / sniff / dump"]:::input --> NORMAL["Normalize source, schema,<br/>columns, predicates, profile"]:::proof
    NORMAL --> NATIVE["PyO3 call<br/>detach GIL"]:::hot
    NATIVE --> CORE["Native read / validate / model / write logic"]:::hot
    CORE --> RESULT{"Native result"}:::decision
    RESULT -- document --> DOC["Frozen CifDocument"]:::output
    RESULT -- table --> TABLE["Frozen CifTable<br/>Arrow C Stream"]:::output
    RESULT -- semantic --> MODEL["Frozen MmcifModel"]:::output
    RESULT -- report tuple --> REPORT["Immutable ValidationReport"]:::output
    RESULT -- error tuple --> ERROR["Central exception conversion<br/>stable code + structured fields"]:::error

    classDef input fill:#eef2ff,stroke:#4f46e5,color:#1e1b4b;
    classDef decision fill:#f8fafc,stroke:#64748b,color:#0f172a;
    classDef hot fill:#fff1f2,stroke:#e11d48,color:#881337,stroke-width:2px;
    classDef proof fill:#ecfdf5,stroke:#059669,color:#064e3b;
    classDef output fill:#f5f3ff,stroke:#7c3aed,color:#4c1d95,stroke-width:2px;
    classDef error fill:#fef2f2,stroke:#b91c1c,color:#7f1d1d;
```

`chomp` and `feast` are exact aliases of `cif.read` and `cif.scan`. `sniff` and
`dump` dispatch through explicit schema/profile arguments and Nibbler object types.

## Bounded ordered scan

The Python iterator admits at most one source per native worker. Native completion order
may differ from input order; `BTreeMap` is the bounded reorder buffer that restores it.

```mermaid
sequenceDiagram
    autonumber
    participant Caller
    participant ScanResult as Python ScanResult
    participant Jobs as sync_channel(workers)
    participant Pool as Native worker pool
    participant Done as completion channel
    participant Order as BTreeMap reorder buffer

    Caller->>ScanResult: cif.scan(sources, workers=N)
    loop initial admission, at most N
        ScanResult->>Jobs: (source_index, normalized source)
    end
    Pool->>Jobs: receive one job
    Pool->>Pool: execute_read without Python objects
    Pool->>Done: (source_index, result)
    ScanResult->>Done: wait with GIL detached
    ScanResult->>Order: buffer out-of-order completions
    Order-->>ScanResult: next expected source_index
    ScanResult->>Jobs: admit one replacement source
    alt success
        ScanResult-->>Caller: CifDocument or CifTable
    else on_error = collect
        ScanResult->>ScanResult: append ordered BatchError
    else on_error = raise
        ScanResult->>Pool: cancel admission
        ScanResult-->>Caller: raise first input-ordered BatchError
    end
```

Stopping early calls `ScanResult.close()`, closes admission, and joins workers while the
GIL is detached. Tables returned by scan enable source/block/frame provenance columns at
Arrow export.

## Transactional path writing

Streams are written directly and checked for short writes. Filesystem destinations use
a sibling temporary file so validation failure never replaces an existing destination.

```mermaid
flowchart TB
    VALUE["CifDocument or MmcifModel"]:::input --> VALIDATE{"Requested validation"}:::decision
    VALIDATE -- profile --> PROFILE["Validate semantic model"]:::proof
    VALIDATE -- dictionary --> DICT["Validate attached dictionary"]:::proof
    VALIDATE -- syntax/none --> BUILD
    PROFILE --> BUILD["Select source-preserving or<br/>canonical semantic document"]:::hot
    DICT --> BUILD
    BUILD --> FORMAT{"CIF or BinaryCIF<br/>plus optional path gzip"}:::decision
    FORMAT --> BYTES["Complete payload bytes"]:::data
    BYTES --> DEST{"Destination kind"}:::decision
    DEST -- binary stream --> STREAM["write(payload)<br/>reject short write"]:::output
    DEST -- filesystem path --> TEMP["Create sibling temp file"]:::data
    TEMP --> SYNC["write + flush + fsync"]:::proof
    SYNC --> REPARSE{"Validation requested?"}:::decision
    REPARSE -- yes --> CHECK["Reparse temp<br/>optional dictionary validation"]:::proof
    REPARSE -- no --> COMMIT
    CHECK -- valid --> COMMIT["os.replace temp -> destination"]:::output
    CHECK -- invalid --> CLEAN["Delete temp<br/>leave destination unchanged"]:::error

    classDef input fill:#eef2ff,stroke:#4f46e5,color:#1e1b4b;
    classDef decision fill:#f8fafc,stroke:#64748b,color:#0f172a;
    classDef hot fill:#fff1f2,stroke:#e11d48,color:#881337,stroke-width:2px;
    classDef proof fill:#ecfdf5,stroke:#059669,color:#064e3b;
    classDef data fill:#ecfeff,stroke:#0891b2,color:#164e63;
    classDef output fill:#f5f3ff,stroke:#7c3aed,color:#4c1d95,stroke-width:2px;
    classDef error fill:#fef2f2,stroke:#b91c1c,color:#7f1d1d;
```
