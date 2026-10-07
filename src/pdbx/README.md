# PDBx semantic algorithms

This directory builds and validates an immutable coordinate-and-chemistry model above
the generic [`CifDocument`](../cif/README.md). It never reparses text and never fetches
chemical data from the network.

## Model construction

Every category decoder resolves its fixed columns once per scalar or loop occurrence.
Rows then read directly by column index; `_atom_site` does not allocate a map or wrapper
object for each atom.

```mermaid
flowchart TB
    DOC["CifDocument<br/>exactly one data block"]:::input --> VIEW["Borrowed CategoryView"]:::data
    VIEW --> FIELDS["Resolve fixed columns once<br/>per scalar or loop occurrence"]:::proof
    FIELDS --> DECODE["Typed category decoders<br/>entry + audit_conform<br/>entities + sequences + branches + schemes<br/>asymmetric units + atom_site + struct_conn<br/>embedded chem_comp atoms + bonds"]:::hot
    DECODE --> RECORDS["Immutable semantic records"]:::data
    DECODE --> REFS["Referenced and embedded<br/>component definitions"]:::data
    REFS --> RESOLVE["Deterministic component resolution"]:::proof
    LOCAL["Optional immutable local CCD registry"]:::input --> RESOLVE
    BUILTIN["Small built-in amino acid,<br/>water, and ion registry"]:::data --> RESOLVE

    RECORDS --> MODEL["PdbxModel"]:::output
    RESOLVE --> MODEL
    DOC --> SHARED["Shared CifDocument clone<br/>Arc-backed blocks"]:::data --> MODEL

    classDef input fill:#eef2ff,stroke:#4f46e5,color:#1e1b4b;
    classDef hot fill:#fff1f2,stroke:#e11d48,color:#881337,stroke-width:2px;
    classDef proof fill:#ecfdf5,stroke:#059669,color:#064e3b;
    classDef data fill:#ecfeff,stroke:#0891b2,color:#164e63;
    classDef output fill:#f5f3ff,stroke:#7c3aed,color:#4c1d95,stroke-width:2px;
```

Construction checks required identifiers and typed values. Full profile coherence is a
separate validation step, so a model may retain an explicit unresolved component and
report it deterministically.

## Component resolution

Embedded definitions have source priority, but a conflicting caller-supplied definition
is an error rather than a silent override.

```mermaid
flowchart TB
    ID["Referenced component ID"]:::input --> EMB{"Embedded definition?"}:::decision
    EMB -- yes --> BOTH{"Local definition also present?"}:::decision
    BOTH -- yes --> SAME{"Definitions compatible?"}:::proof
    SAME -- no --> CONFLICT["PDBX_COMPONENT_CONFLICT"]:::error
    SAME -- yes --> USE_EMB["Use embedded definition"]:::output
    BOTH -- no --> USE_EMB

    EMB -- no --> LOCAL{"Local CCD definition?"}:::decision
    LOCAL -- yes --> USE_LOCAL["Use LocalCcd definition"]:::output
    LOCAL -- no --> MIN{"Built-in minimal definition?"}:::decision
    MIN -- yes --> USE_MIN["Use MinimalRegistry definition"]:::output
    MIN -- no --> UNRES["Retain Unresolved placeholder"]:::fallback
    UNRES --> PROFILE["Strict profile emits<br/>PDBX_COMPONENT_UNRESOLVED"]:::error

    classDef input fill:#eef2ff,stroke:#4f46e5,color:#1e1b4b;
    classDef decision fill:#f8fafc,stroke:#64748b,color:#0f172a;
    classDef proof fill:#ecfdf5,stroke:#059669,color:#064e3b;
    classDef fallback fill:#fff7ed,stroke:#ea580c,color:#7c2d12,stroke-width:2px,stroke-dasharray:5 3;
    classDef output fill:#f5f3ff,stroke:#7c3aed,color:#4c1d95,stroke-width:2px;
    classDef error fill:#fef2f2,stroke:#b91c1c,color:#7f1d1d;
```

## Validation and canonical output

```mermaid
flowchart TB
    DOC["CifDocument"]:::input --> DDL2["Pinned PDBx 5.416<br/>dictionary validation"]:::proof
    DOC --> BUILD["Build PdbxModel"]:::hot
    BUILD --> INDEX["Case-normalized indices<br/>entities, asyms, components, atoms"]:::data
    INDEX --> SEM["Identifiers + cross-references<br/>schemes + chemistry + connections"]:::hot
    DDL2 --> MERGE["Stable capped diagnostics"]:::proof
    SEM --> MERGE
    BUILD -. structural failure .-> MERGE
    MERGE --> REPORT["ValidationReport<br/>explicit PDBx coverage"]:::output

    BUILD --> CANON["Clone source entries<br/>append resolved component rows"]:::hot
    CANON --> ORDER["Stable profile category order<br/>canonical atom_site columns"]:::proof
    ORDER --> CDOC["Canonical CifDocument"]:::output
    CDOC --> CIFWRITE["Shared CIF/BinaryCIF writer"]:::output

    classDef input fill:#eef2ff,stroke:#4f46e5,color:#1e1b4b;
    classDef hot fill:#fff1f2,stroke:#e11d48,color:#881337,stroke-width:2px;
    classDef proof fill:#ecfdf5,stroke:#059669,color:#064e3b;
    classDef data fill:#ecfeff,stroke:#0891b2,color:#164e63;
    classDef output fill:#f5f3ff,stroke:#7c3aed,color:#4c1d95,stroke-width:2px;
```

The writer preserves all source categories, adds non-embedded resolved component
definitions when required, orders known profile categories, and delegates value quoting
and byte serialization to the generic CIF layer.

## Biological assemblies

`assembly_document` writes one `_pdbx_struct_assembly` out as explicit copies. Operator
expressions expand to combinations (leftmost group outermost, rightmost operator applied
first); each distinct combination is one copy. Copies of an atom that land within 0.2 Å
of each other (an atom on a symmetry axis) are written once; each atom's copies are hashed
into cells of that size, so the check grows with the atoms written. Rows of the per-copy
categories are streamed into text columns, so a large assembly does not allocate a value
object per cell.

```mermaid
flowchart TB
    MODEL["PdbxModel<br/>source document"]:::input --> GEN["_pdbx_struct_assembly_gen<br/>rows of the selected assembly"]:::data
    OPER["_pdbx_struct_oper_list<br/>rotation + translation"]:::data --> COMBO
    GEN --> COMBO["Expand oper_expression<br/>(1,2) (1-60)(61) ..."]:::proof
    COMBO --> COPIES["Copies: distinct combinations<br/>composed transform + asym ids"]:::data
    COPIES --> NAMES{"Renamed chains<br/>X-n unique?"}:::decision
    NAMES -- no --> COLLIDE["PDBX_ASSEMBLY_ID_COLLISION"]:::error
    NAMES -- yes --> MERGE["Copies of an atom within 0.2 Å<br/>written once; absorbed chain copies<br/>left out, connections redirected"]:::proof
    MERGE --> ROWS["Per copy: struct_asym, schemes,<br/>atom_site (transformed), anisotrop (rotated),<br/>same-copy struct_conn"]:::hot
    ROWS --> DOC["New one-block CifDocument<br/>no cell, symmetry, or assembly categories"]:::output

    classDef input fill:#eef2ff,stroke:#4f46e5,color:#1e1b4b;
    classDef hot fill:#fff1f2,stroke:#e11d48,color:#881337,stroke-width:2px;
    classDef proof fill:#ecfdf5,stroke:#059669,color:#064e3b;
    classDef data fill:#ecfeff,stroke:#0891b2,color:#164e63;
    classDef decision fill:#f8fafc,stroke:#64748b,color:#0f172a;
    classDef output fill:#f5f3ff,stroke:#7c3aed,color:#4c1d95,stroke-width:2px;
    classDef error fill:#fef2f2,stroke:#b91c1c,color:#7f1d1d;
```
