# ModelCIF semantic algorithms

ModelCIF extends the shared [`PdbxModel`](../pdbx/README.md); it does not create a second
coordinate graph. This directory owns prediction, provenance, template/alignment, QA,
file/archive, profile-validation, and optional viewer-mirroring semantics.

## Construction and validation

```mermaid
flowchart TB
    DOC["CifDocument"]:::input --> PDBX["Build shared PdbxModel<br/>coordinates + chemistry"]:::hot
    DOC --> VIEW["Borrowed CategoryView"]:::data

    VIEW --> PRED["prediction.rs<br/>targets, models, groups, data flow"]:::hot
    VIEW --> PROV["provenance.rs<br/>software, protocols, files, archives"]:::hot
    VIEW --> TEMPLATE["template.rs<br/>templates, mappings, alignments"]:::hot
    VIEW --> QA["qa.rs<br/>metric definitions + global/local/pairwise values"]:::hot

    PDBX --> MODEL["ModelCifModel"]:::output
    PRED --> MODEL
    PROV --> MODEL
    TEMPLATE --> MODEL
    QA --> MODEL

    MODEL --> INDICES["Stable ID/reference indices"]:::data
    INDICES --> VPRED["prediction checks"]:::proof
    INDICES --> VPROV["provenance/data-flow checks"]:::proof
    INDICES --> VTEMP["template/alignment checks"]:::proof
    INDICES --> VQA["QA definition/value checks"]:::proof
    MODEL --> VPDBX["PDBx semantic checks"]:::proof

    VPRED --> REPORT["ModelCIF ValidationReport"]:::output
    VPROV --> REPORT
    VTEMP --> REPORT
    VQA --> REPORT
    VPDBX --> REPORT

    classDef input fill:#eef2ff,stroke:#4f46e5,color:#1e1b4b;
    classDef hot fill:#fff1f2,stroke:#e11d48,color:#881337,stroke-width:2px;
    classDef proof fill:#ecfdf5,stroke:#059669,color:#064e3b;
    classDef data fill:#ecfeff,stroke:#0891b2,color:#164e63;
    classDef output fill:#f5f3ff,stroke:#7c3aed,color:#4c1d95,stroke-width:2px;
```

The four validation slices share one index-owning validator. Keeping the slices separate
makes each domain's references visible without duplicating indices or diagnostics.

## Canonical writing and explicit QA mirroring

Local prediction confidence is not a crystallographic displacement value. The normal
writer leaves B factors alone. Mirroring is a separate, explicit output view.

```mermaid
flowchart TB
    MODEL["ModelCifModel"]:::input --> MODE{"Output request"}:::decision
    MODE -- normal --> ORDER["PDBx canonical builder<br/>with ModelCIF category order"]:::proof
    ORDER --> DOC["Canonical CifDocument<br/>QA categories intact"]:::output

    MODE -- mirror local metric ID --> FIND{"Metric exists and mode=local?"}:::proof
    FIND -- no --> ERR["Structured mirror error"]:::error
    FIND -- yes --> MAP["Index local QA values by<br/>model, asym, seq, component"]:::data
    MAP --> UNIQUE{"At least one value<br/>and no duplicate identity?"}:::proof
    UNIQUE -- no --> ERR
    UNIQUE -- yes --> ORDER
    ORDER --> ATOMS["Locate atom_site identity columns"]:::hot
    ATOMS --> REWRITE["Copy matching confidence into<br/>B_iso_or_equiv output cells only"]:::hot
    REWRITE --> MIRRORED["Viewer-oriented CifDocument<br/>original QA records preserved"]:::output

    classDef input fill:#eef2ff,stroke:#4f46e5,color:#1e1b4b;
    classDef decision fill:#f8fafc,stroke:#64748b,color:#0f172a;
    classDef hot fill:#fff1f2,stroke:#e11d48,color:#881337,stroke-width:2px;
    classDef proof fill:#ecfdf5,stroke:#059669,color:#064e3b;
    classDef data fill:#ecfeff,stroke:#0891b2,color:#164e63;
    classDef output fill:#f5f3ff,stroke:#7c3aed,color:#4c1d95,stroke-width:2px;
    classDef error fill:#fef2f2,stroke:#b91c1c,color:#7f1d1d;
```

Preserving mode cannot mirror a metric. Pairwise or global metrics cannot be selected,
and a missing `B_iso_or_equiv` column is added only to the constructed output document.
