# CIF core algorithms

This directory owns generic CIF syntax, input decoding, projection, BinaryCIF,
dictionaries, validation, Arrow export, and generic serialization. The diagrams show the
current production paths. Rose nodes are CPU-hot, green nodes are proof or validation
gates, orange dashed nodes are conservative fallbacks, cyan nodes retain data, and
violet nodes are public outputs.

## Read dispatch and shared outputs

Gzip and BinaryCIF are detected from bytes, not filename suffixes. Gzip decoding happens
before the text/BinaryCIF split, so compressed BinaryCIF follows the binary branch.

```mermaid
flowchart TB
    INPUT["Path or owned bytes"]:::input --> LIMIT["Input-size check"]:::proof
    LIMIT --> GZIP{"gzip magic?"}:::decision
    GZIP -- yes --> DEFLATE["zlib-rs decode<br/>size + ratio limits"]:::hot
    GZIP -- no --> BYTES["Owned bytes"]:::data
    DEFLATE --> BYTES
    BYTES --> BINARY{"MessagePack map prefix?"}:::decision

    BINARY -- no --> UTF8["UTF-8 validation<br/>shared SourceBuffer"]:::data
    UTF8 --> TEXTMODE{"category requested?"}:::decision
    TEXTMODE -- no --> PARSE["Lexer + Parser<br/>DocumentSink"]:::hot
    TEXTMODE -- yes --> PROJECT["ProjectionPlan<br/>Lexer + Parser + TableSink"]:::hot

    BINARY -- yes --> BINMODE{"category requested?"}:::decision
    BINMODE -- no --> BINDOC["BinaryCIF shape check<br/>decode all categories"]:::hot
    BINMODE -- yes --> BINPROJ["BinaryCIF shape check<br/>decode selected columns"]:::hot

    PARSE --> DOC["CifDocument"]:::output
    BINDOC --> DOC
    PROJECT --> TABLE["CifTable"]:::output
    BINPROJ --> TABLE

    classDef input fill:#eef2ff,stroke:#4f46e5,color:#1e1b4b;
    classDef decision fill:#f8fafc,stroke:#64748b,color:#0f172a;
    classDef hot fill:#fff1f2,stroke:#e11d48,color:#881337,stroke-width:2px;
    classDef proof fill:#ecfdf5,stroke:#059669,color:#064e3b;
    classDef data fill:#ecfeff,stroke:#0891b2,color:#164e63;
    classDef output fill:#f5f3ff,stroke:#7c3aed,color:#4c1d95,stroke-width:2px;
```

Code owners: `input.rs`, `source.rs`, `parser.rs`, `projection.rs`, and `binary/`.

## One text grammar, two sinks

`Parser` is the sole owner of block, frame, scalar, loop, duplicate-tag, row-completeness,
control-word, and resource-limit rules. `ParseSink` receives events only after the parser
accepts the corresponding grammar transition.

```mermaid
flowchart TB
    SOURCE["SourceBuffer"]:::data --> LEXER["Lexer::next_token<br/>span + kind + quote style"]:::hot
    LEXER --> LOOK["Parser one-token lookahead"]:::hot
    LOOK --> GRAMMAR{"Current grammar scope"}:::decision

    subgraph OWNED["Parser-owned grammar"]
        DOCSCOPE["Document<br/>data_ or global_ only"]:::proof
        BLOCK["Block<br/>item, loop, frame, next block"]:::proof
        FRAME["Save frame<br/>item, loop, save_ end"]:::proof
        LOOP["Loop<br/>tags, values, optional stop_"]:::proof
        DOCSCOPE --> BLOCK
        BLOCK --> FRAME
        BLOCK --> LOOP
        FRAME --> LOOP
    end

    GRAMMAR --> OWNED
    OWNED --> CHECKS["Charge limits before retention<br/>reject duplicates and incomplete rows"]:::proof
    CHECKS --> EVENTS["ParseSink events"]:::input

    EVENTS --> DS["DocumentSink"]:::hot
    DS --> SCALAR["Scalar CifValue<br/>shared source slice"]:::data
    DS --> CELLS["Loop SourceCell chunks<br/>8 bytes per value"]:::data
    SCALAR --> DOCUMENT["CifDocument"]:::output
    CELLS --> DOCUMENT

    EVENTS --> TS["TableSink"]:::hot
    TS --> SELECT["Selected + predicate columns only"]:::data
    SELECT --> TABLE["Segmented CifTable"]:::output

    FAIL["Located ParseError<br/>stable code + earliest span"]:::error
    LEXER -. lexical failure .-> FAIL
    CHECKS -. grammar or limit failure .-> FAIL

    classDef input fill:#eef2ff,stroke:#4f46e5,color:#1e1b4b;
    classDef decision fill:#f8fafc,stroke:#64748b,color:#0f172a;
    classDef hot fill:#fff1f2,stroke:#e11d48,color:#881337,stroke-width:2px;
    classDef proof fill:#ecfdf5,stroke:#059669,color:#064e3b;
    classDef data fill:#ecfeff,stroke:#0891b2,color:#164e63;
    classDef output fill:#f5f3ff,stroke:#7c3aed,color:#4c1d95,stroke-width:2px;
    classDef error fill:#fef2f2,stroke:#b91c1c,color:#7f1d1d;
```

The ordinary unquoted-token path fuses delimiter search and forbidden-character
validation in one scalar pass. Quoted values, comments, semicolon text fields, and
control words remain in the same lexer.

## Large-loop parallel kernel

This is the main text-CIF throughput path. It is an optimization under the parser, not
a second grammar. Default/custom resource limits that could alter error ordering keep a
loop serial.

### Admission and safe ranges

```mermaid
flowchart TB
    FIRST["Parser has loop tags<br/>and first value offset"]:::input --> SAFE_LIMITS{"Parallel limits<br/>provably unreachable?"}:::proof
    SAFE_LIMITS -- no --> SERIAL["Serial parser loop"]:::fallback
    SAFE_LIMITS -- yes --> SIZE{"At least 8 MiB remains?"}:::proof
    SIZE -- no --> SERIAL

    SIZE -- yes --> MODE{"Sink kernel"}:::decision
    MODE -- retain document --> PROBE{"Lightweight probe proves<br/>this loop reaches 8 MiB?"}:::proof
    PROBE -- no --> SERIAL
    PROBE -- yes --> LEASE
    MODE -- project table --> LEASE["Atomic process-wide worker lease<br/>ceil(bytes / 2 MiB), max host CPUs"]:::proof
    LEASE --> WORKERS{"At least two workers granted?"}:::proof
    WORKERS -- no --> SERIAL
    WORKERS -- yes --> BOUNDS["safe_chunks<br/>boundaries outside text fields"]:::hot
    BOUNDS --> READY["Worker lease + RawChunk ranges"]:::output

    classDef input fill:#eef2ff,stroke:#4f46e5,color:#1e1b4b;
    classDef decision fill:#f8fafc,stroke:#64748b,color:#0f172a;
    classDef hot fill:#fff1f2,stroke:#e11d48,color:#881337,stroke-width:2px;
    classDef proof fill:#ecfdf5,stroke:#059669,color:#064e3b;
    classDef fallback fill:#fff7ed,stroke:#ea580c,color:#7c2d12,stroke-width:2px,stroke-dasharray:5 3;
    classDef output fill:#f5f3ff,stroke:#7c3aed,color:#4c1d95,stroke-width:2px;
```

`parallel.rs` owns worker leases, validation, prefix counts, error arbitration, segmented
retention, speculative projection, and the row-aligned fallback. `parallel_probe.rs`
owns the allocation-free document eligibility proof.

### Parallel scan and error arbitration

```mermaid
flowchart TB
    READY["Worker lease + safe RawChunk ranges"]:::input --> MODE{"Sink kernel"}:::decision
    MODE -- retain --> RETAIN["Production lexer<br/>validate + count + retain<br/>64K-cell segments"]:::hot
    MODE -- project --> SPEC["Production lexer<br/>project with local phase 0"]:::hot

    RETAIN --> RESULT["Per chunk:<br/>count, last value, control/error, cells"]:::data
    SPEC --> PRESULT["Per chunk:<br/>count, last value, control/error, table"]:::data
    RESULT --> CONTROL["Select earliest control token<br/>as global loop end"]:::proof
    PRESULT --> CONTROL
    CONTROL --> ERROR["Select earliest lexical error<br/>strictly before loop end"]:::proof
    ERROR --> PREFIX["Prefix counts in source order<br/>retain global last value"]:::hot
    PREFIX --> LOOP["ParallelLoop<br/>chunks + total + control + lease"]:::output

    classDef input fill:#eef2ff,stroke:#4f46e5,color:#1e1b4b;
    classDef decision fill:#f8fafc,stroke:#64748b,color:#0f172a;
    classDef hot fill:#fff1f2,stroke:#e11d48,color:#881337,stroke-width:2px;
    classDef proof fill:#ecfdf5,stroke:#059669,color:#064e3b;
    classDef data fill:#ecfeff,stroke:#0891b2,color:#164e63;
    classDef output fill:#f5f3ff,stroke:#7c3aed,color:#4c1d95,stroke-width:2px;
```

### Proof and commit

```mermaid
flowchart TB
    LOOP["ParallelLoop + optional<br/>speculative tables"]:::input --> ROWS{"Total values nonzero and<br/>divisible by tag count?"}:::proof
    ROWS -- no --> ERROR["Located missing/incomplete-row error"]:::error
    ROWS -- yes --> KIND{"Sink kernel"}:::decision

    KIND -- retain --> CELLS["Commit source-ordered<br/>SourceCell segments"]:::output
    KIND -- project --> DIRECT{"Every chunk prefix mod<br/>tag count = 0?"}:::proof
    DIRECT -- yes --> TABLES["Commit speculative<br/>CifTable segments"]:::output
    DIRECT -- no --> DISCARD["Discard speculative tables<br/>and typed errors"]:::fallback
    DISCARD --> GENERAL["General validation/count scan"]:::fallback
    GENERAL --> ALIGNED["Advance misphased starts<br/>to row boundaries"]:::proof
    ALIGNED --> ENOUGH{"At least two aligned chunks?"}:::proof
    ENOUGH -- no --> SERIAL["Serial parser loop"]:::fallback
    ENOUGH -- yes --> REPROJECT["Parallel row-aligned projection"]:::hot
    REPROJECT --> TABLES

    CELLS --> HANDOFF["Advance parser lexer to<br/>stop_, next control, or EOF"]:::proof
    TABLES --> HANDOFF
    SERIAL --> HANDOFF

    classDef input fill:#eef2ff,stroke:#4f46e5,color:#1e1b4b;
    classDef decision fill:#f8fafc,stroke:#64748b,color:#0f172a;
    classDef hot fill:#fff1f2,stroke:#e11d48,color:#881337,stroke-width:2px;
    classDef proof fill:#ecfdf5,stroke:#059669,color:#064e3b;
    classDef fallback fill:#fff7ed,stroke:#ea580c,color:#7c2d12,stroke-width:2px,stroke-dasharray:5 3;
    classDef output fill:#f5f3ff,stroke:#7c3aed,color:#4c1d95,stroke-width:2px;
    classDef error fill:#fef2f2,stroke:#b91c1c,color:#7f1d1d;
```

## Safe boundary index

Quoted values cannot cross a line and comments end at a line ending. Only semicolon text
fields can cross a candidate line boundary, so `parallel_index.rs` summarizes the exact
two-state text-field machine instead of lexing every line serially.

```mermaid
flowchart TB
    TARGETS["Nominal byte-balanced<br/>worker targets"]:::input --> ALIGN["Move each target to<br/>the next line start"]:::proof
    ALIGN --> BLOCKS["Boundary scan blocks"]:::data

    subgraph PAR["Parallel memchr(';') summaries"]
        B0["Block 0"]:::hot
        B1["Block 1"]:::hot
        BN["Block N"]:::hot
    end
    BLOCKS --> B0
    BLOCKS --> B1
    BLOCKS --> BN

    B0 --> S0["end_state[outside, inside]<br/>first_safe[outside, inside]"]:::data
    B1 --> S1["same two-state summary"]:::data
    BN --> SN["same two-state summary"]:::data

    S0 --> PREFIX["Compose incoming text-field state<br/>in source order"]:::proof
    S1 --> PREFIX
    SN --> PREFIX
    PREFIX --> PICK["For each target, choose first<br/>safe boundary under true state"]:::proof
    PICK --> RANGES["Non-overlapping RawChunk ranges"]:::output

    classDef input fill:#eef2ff,stroke:#4f46e5,color:#1e1b4b;
    classDef hot fill:#fff1f2,stroke:#e11d48,color:#881337,stroke-width:2px;
    classDef proof fill:#ecfdf5,stroke:#059669,color:#064e3b;
    classDef data fill:#ecfeff,stroke:#0891b2,color:#164e63;
    classDef output fill:#f5f3ff,stroke:#7c3aed,color:#4c1d95,stroke-width:2px;
```

Each summary is computed for both possible incoming states. The coordinator therefore
composes exact results without guessing whether an earlier block opened a text field.

## Document and table storage

The shared logical API hides three loop representations. Each representation matches
its source rather than forcing every format through a row-major object graph.

```mermaid
flowchart TB
    subgraph TEXT["Parsed text CIF"]
        SRC["SourceBuffer<br/>Arc-owned UTF-8 + name"]:::data
        LAZY["OnceLock line index<br/>built only for diagnostics"]:::data
        SCALAR["Scalar TextValue<br/>source + content range"]:::data
        CELL["SourceCell<br/>u32 start + u32 end"]:::hot
        SEG["SourceCells<br/>64K-value chunks + offsets"]:::data
        SRC --> LAZY
        SRC --> SCALAR
        SRC --> CELL --> SEG
    end

    subgraph LOOPS["CifLoop storage variants"]
        SOURCE["Source<br/>shared text + SourceCells"]:::data
        COLUMNS["Columns<br/>typed vectors + masks"]:::data
        OWNED["Owned<br/>Vec<CifValue>"]:::data
    end
    SRC --> SOURCE
    SEG --> SOURCE
    BCIF["BinaryCIF decoder"]:::input --> COLUMNS
    BUILD["Constructed/canonical document"]:::input --> OWNED

    SOURCE --> LOOP["CifLoop logical rows/values"]:::output
    COLUMNS --> LOOP
    OWNED --> LOOP
    LOOP --> BLOCKS["Arc<[CifBlock]><br/>O(1) document clone"]:::output

    subgraph TABLE["Projected CifTable"]
        ARRAYS["Segmented Arrow value arrays"]:::data
        KINDS["Segmented UInt8 missing-kind arrays"]:::data
        PROV["Run-length source/block/frame provenance"]:::data
    end
    ARRAYS --> CT["CifTable"]:::output
    KINDS --> CT
    PROV --> CT

    classDef input fill:#eef2ff,stroke:#4f46e5,color:#1e1b4b;
    classDef hot fill:#fff1f2,stroke:#e11d48,color:#881337,stroke-width:2px;
    classDef data fill:#ecfeff,stroke:#0891b2,color:#164e63;
    classDef output fill:#f5f3ff,stroke:#7c3aed,color:#4c1d95,stroke-width:2px;
```

Text loop values recover quote style, missing kind, and content bounds lazily from the
two stored offsets. Sequential iteration walks chunks directly; random access uses the
chunk-offset index.

## Projection and Arrow export

Selectors are normalized before parsing. Predicate columns are assigned indexes once
per matching category occurrence, so the row loop does not search names or call Python.

```mermaid
flowchart TB
    REQUEST["category + columns + where + schema"]:::input --> PLAN["ProjectionPlan<br/>normalize category/items"]:::proof
    SCHEMA["Optional compiled dictionary"]:::data --> TYPES["Resolve UTF-8 / int64 / float64"]:::proof
    TYPES --> PLAN
    PLAN --> OCC["Matching scalar or loop occurrence"]:::hot
    OCC --> PREP["Check layout + required items<br/>compile target/predicate indexes"]:::proof
    PREP --> ROWS["Validate every token<br/>decode selected cells only"]:::hot
    ROWS --> PRED{"Predicates match?"}:::decision
    PRED -- no --> DROP["Drop row"]:::fallback
    PRED -- yes --> BUILD["Append typed value array<br/>and exact missing-kind byte"]:::hot
    BUILD --> SEGMENTS["Source-ordered Arrow chunks<br/>no final concatenation"]:::data

    SEGMENTS --> POLICY{"Arrow missing policy"}:::decision
    POLICY -- collapse --> NULLS["Arrow nulls<br/>counts in field metadata"]:::output
    POLICY -- columns --> COMP["Value + non-null<br/>__missing_kind column"]:::output
    POLICY -- extension --> EXT["struct<value, kind><br/>nibbler.cif_missing"]:::output
    SEGMENTS --> PROVENANCE["Scan only:<br/>source/block/frame columns"]:::output

    classDef input fill:#eef2ff,stroke:#4f46e5,color:#1e1b4b;
    classDef decision fill:#f8fafc,stroke:#64748b,color:#0f172a;
    classDef hot fill:#fff1f2,stroke:#e11d48,color:#881337,stroke-width:2px;
    classDef proof fill:#ecfdf5,stroke:#059669,color:#064e3b;
    classDef fallback fill:#fff7ed,stroke:#ea580c,color:#7c2d12,stroke-dasharray:5 3;
    classDef data fill:#ecfeff,stroke:#0891b2,color:#164e63;
    classDef output fill:#f5f3ff,stroke:#7c3aed,color:#4c1d95,stroke-width:2px;
```

`arrow.rs` emits one record batch per retained chunk through the Arrow C Stream release
callback. PyArrow and Polars take ownership of those arrays without Nibbler joining the
payload buffers.

## BinaryCIF 0.3

BinaryCIF has a separate MessagePack codec but converges on the same document, table,
missing-state, schema, and writer contracts.

```mermaid
flowchart TB
    RAW["MessagePack bytes"]:::input --> DECODE["Deserialize BinaryFile"]:::hot
    DECODE --> SHAPE["Require 0.3.x + data blocks<br/>non-empty names with a text CIF spelling<br/>no case-folded repeats per scope<br/>columns for every category with rows<br/>checked counts + <=1000x expansion"]:::proof
    SHAPE --> ROWS["Drop zero-row categories<br/>text CIF has no empty loop"]:::proof
    ROWS --> MODE{"Requested output"}:::decision

    MODE -- document --> ALL["Reverse every encoding chain<br/>validate row counts and masks"]:::hot
    ALL --> COLS["ColumnValues<br/>int64 / float64 / string dictionary"]:::data
    COLS --> LOOPS["CifLoop::Columns"]:::data
    LOOPS --> DOC["CifDocument"]:::output

    MODE -- projection --> MATCH["Visit matching categories only"]:::hot
    MATCH --> SELECT["Decode output + predicate columns only"]:::hot
    SELECT --> FILTER["Shared ProjectionPlan predicates<br/>shared TableBuilder"]:::proof
    FILTER --> TABLE["CifTable"]:::output

    DOC --> WRITE["Reject global blocks, frames,<br/>empty categories, lossy values"]:::proof
    WRITE --> ENCODE["Int32, float64, or StringArray<br/>plus 0/1/2 mask"]:::hot
    ENCODE --> OUT["Deterministic BinaryCIF 0.3"]:::output

    classDef input fill:#eef2ff,stroke:#4f46e5,color:#1e1b4b;
    classDef decision fill:#f8fafc,stroke:#64748b,color:#0f172a;
    classDef hot fill:#fff1f2,stroke:#e11d48,color:#881337,stroke-width:2px;
    classDef proof fill:#ecfdf5,stroke:#059669,color:#064e3b;
    classDef data fill:#ecfeff,stroke:#0891b2,color:#164e63;
    classDef output fill:#f5f3ff,stroke:#7c3aed,color:#4c1d95,stroke-width:2px;
```

## Dictionaries and validation

Built-in dictionaries are deterministic artifacts generated from lock-pinned DDL2
sources. Runtime loading is local, lazy, and immutable.

```mermaid
flowchart TB
    ART["Embedded .nbs artifact"]:::input --> ENVELOPE["Magic + format v2<br/>16 MiB limit + no trailing bytes"]:::proof
    ENVELOPE --> LOCK["Schema name + version + source SHA<br/>match compiled lock"]:::proof
    LOCK --> REGEX["Compile anchored DDL2 type regexes"]:::hot
    REGEX --> ONCE["OnceLock<LoadedSchema>"]:::data

    DOC["CifDocument"]:::input --> VIEW["Borrowed CategoryView<br/>scalar + distinct loop occurrences"]:::data
    VIEW --> RESOLVE["Resolve aliases and item definitions"]:::proof
    ONCE --> RESOLVE
    RESOLVE --> INDEX["Index parent values<br/>without copying columns"]:::hot
    INDEX --> RULES["Mandatory categories/items<br/>types, enums, ranges, keys, parents"]:::hot
    RULES --> CAP["Stable source order<br/>10,000 diagnostics + truncation"]:::proof
    CAP --> REPORT["ValidationReport"]:::output

    ONCE --> PROJECTION["Projection typing<br/>and selector validation"]:::output

    classDef input fill:#eef2ff,stroke:#4f46e5,color:#1e1b4b;
    classDef hot fill:#fff1f2,stroke:#e11d48,color:#881337,stroke-width:2px;
    classDef proof fill:#ecfdf5,stroke:#059669,color:#064e3b;
    classDef data fill:#ecfeff,stroke:#0891b2,color:#164e63;
    classDef output fill:#f5f3ff,stroke:#7c3aed,color:#4c1d95,stroke-width:2px;
```

The same `CategoryView` feeds generic validation and the semantic decoders in
[`../pdbx/`](../pdbx/README.md) and [`../modelcif/`](../modelcif/README.md).
