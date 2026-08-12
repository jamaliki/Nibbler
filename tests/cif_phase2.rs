//! Projection, Arrow-storage, gzip, and resource-limit tests for Phase 2.
//!
//! Assertions intentionally use panic-based helpers, as permitted for tests by
//! `ENGINEERING.md`; production code retains the workspace-wide deny policy.

#![allow(clippy::expect_used)]

use std::io::Write;

use _core::cif::{
    CifCellRef, InputErrorCode, InputLimits, Limits, MissingKind, ParseOptions, Predicate,
    ProjectionErrorCode, ProjectionPlan, SourceBuffer, decode_source, parse, parse_with_options,
    project, project_source, project_with_options,
};
use flate2::Compression;
use flate2::write::GzEncoder;

const CHEMISTRY: &[u8] = include_bytes!("fixtures/chemistry/ligand_ion_water.cif");
const MISSING_VALUES: &[u8] = include_bytes!("fixtures/syntax/missing_values.cif");
const MALFORMED_LOOP: &[u8] = include_bytes!("fixtures/syntax/malformed_loop.cif");

#[test]
fn projects_requested_columns_and_filters_before_retaining_rows() {
    let plan = ProjectionPlan::new("atom_site")
        .and_then(|plan| plan.with_columns(["label_comp_id", "Cartn_x"]))
        .and_then(|plan| {
            plan.with_predicate(Predicate::Equal {
                column: "pdbx_PDB_model_num".to_owned(),
                value: "1".to_owned(),
            })
        })
        .and_then(|plan| {
            plan.with_predicate(Predicate::In {
                column: "group_PDB".to_owned(),
                values: vec!["ATOM".to_owned()],
            })
        })
        .expect("projection plan is valid");

    let table = project(CHEMISTRY, plan).expect("fixture must project");
    assert_eq!(table.category(), "atom_site");
    assert_eq!(
        table
            .columns()
            .iter()
            .map(|column| column.name())
            .collect::<Vec<_>>(),
        ["label_comp_id", "Cartn_x"]
    );
    assert_eq!(table.row_count(), 2);
    assert_eq!(
        table
            .column("LABEL_COMP_ID")
            .expect("column exists")
            .values()
            .collect::<Vec<_>>(),
        [CifCellRef::Text("ALA"), CifCellRef::Text("MSE")]
    );
    assert!(
        table
            .columns()
            .iter()
            .all(|column| column.values().len() == table.row_count())
    );
}

#[test]
fn implicit_projection_preserves_source_order_and_both_missing_kinds() {
    let table = project(
        MISSING_VALUES,
        ProjectionPlan::new("_nibbler_missing").expect("category is valid"),
    )
    .expect("fixture must project");

    assert_eq!(table.row_count(), 2);
    assert_eq!(table.columns()[0].name(), "id");
    assert_eq!(
        table
            .column("unknown_value")
            .expect("column exists")
            .values()
            .collect::<Vec<_>>(),
        [CifCellRef::Unknown, CifCellRef::NotApplicable]
    );
    assert_eq!(
        table
            .column("not_applicable_value")
            .expect("column exists")
            .values()
            .collect::<Vec<_>>(),
        [CifCellRef::NotApplicable, CifCellRef::Unknown]
    );
    assert_eq!(
        table
            .column("present_empty_value")
            .expect("column exists")
            .values()
            .collect::<Vec<_>>(),
        [CifCellRef::Text(""), CifCellRef::Text("")]
    );
}

#[test]
fn missing_predicates_distinguish_unknown_from_not_applicable() {
    let plan = ProjectionPlan::new("nibbler_missing")
        .and_then(|plan| plan.with_columns(["id"]))
        .and_then(|plan| {
            plan.with_predicate(Predicate::Missing {
                column: "unknown_value".to_owned(),
                kind: Some(MissingKind::Unknown),
            })
        })
        .expect("projection plan is valid");
    let table = project(MISSING_VALUES, plan).expect("fixture must project");
    assert_eq!(
        table
            .column("id")
            .expect("column exists")
            .values()
            .collect::<Vec<_>>(),
        [CifCellRef::Text("1")]
    );
}

#[test]
fn projection_keeps_source_block_and_frame_provenance() {
    let source = SourceBuffer::from_text(
        "blocks.cif",
        "data_one\nloop_\n_x.id\n1\n2\ndata_two\nsave_frame\n_x.id 3\nsave_\n",
    );
    let table = project_source(source, ProjectionPlan::new("x").expect("category is valid"))
        .expect("source must project");

    assert_eq!(table.row_count(), 3);
    assert_eq!(
        table
            .row_provenance(0)
            .expect("row provenance exists")
            .source_name(),
        "blocks.cif"
    );
    assert_eq!(
        table
            .row_provenance(0)
            .expect("row provenance exists")
            .block_code(),
        Some("one")
    );
    assert_eq!(
        table
            .row_provenance(2)
            .expect("row provenance exists")
            .block_code(),
        Some("two")
    );
    assert_eq!(
        table
            .row_provenance(2)
            .expect("row provenance exists")
            .frame_code(),
        Some("frame")
    );
    assert_eq!(table.provenance().len(), 2);
}

#[test]
fn required_columns_never_turn_into_fabricated_nulls() {
    let plan = ProjectionPlan::new("nibbler_missing")
        .and_then(|plan| plan.with_columns(["id", "absent"]))
        .expect("projection plan is valid");
    let error = project(MISSING_VALUES, plan).expect_err("missing column must fail");
    assert_eq!(error.code(), ProjectionErrorCode::MissingColumn);
}

#[test]
fn syntax_errors_remain_located_parser_errors() {
    let error = project(
        MALFORMED_LOOP,
        ProjectionPlan::new("bad").expect("category is valid"),
    )
    .expect_err("malformed source must fail");
    assert_eq!(error.code(), ProjectionErrorCode::Parse);
    let parse_error = error
        .parse_error()
        .expect("located parse error is retained");
    assert_eq!(parse_error.code_str(), "CIF_LOOP_VALUE_COUNT");
    assert_eq!(parse_error.span().line(), 6);
}

#[test]
fn absent_categories_return_a_well_typed_empty_projection() {
    let plan = ProjectionPlan::new("absent")
        .and_then(|plan| plan.with_columns(["id", "value"]))
        .expect("projection plan is valid");
    let table = project(MISSING_VALUES, plan).expect("absence is not malformed input");
    assert_eq!(table.row_count(), 0);
    assert_eq!(table.columns().len(), 2);
}

#[test]
fn parallel_projection_matches_forced_serial_on_adversarial_loop() {
    let mut source = String::from("data_parallel\nloop_\n_x.id\n_x.payload\n_x.flag\n");
    source.push_str("0\n;\n");
    while source.len() < 3 * 1024 * 1024 {
        source.push_str("a semicolon text field may cross a nominal worker boundary\n");
    }
    source.push_str(";\nyes\n");

    let payload = "x".repeat(1_024);
    let mut row_count = 1;
    let mut retained = 1;
    while source.len() < 9 * 1024 * 1024 {
        let flag = if row_count % 3 == 0 { "yes" } else { "no" };
        let value = if row_count == 7 { "'stop_'" } else { &payload };
        source.push_str(&format!("{row_count}\n{value}\n{flag}\n"));
        retained += usize::from(flag == "yes");
        row_count += 1;
    }
    source.push_str("stop_\n_y.after parallel-control-handoff\n");

    let plan = ProjectionPlan::new("x")
        .and_then(|plan| plan.with_columns(["id"]))
        .and_then(|plan| {
            plan.with_predicate(Predicate::Equal {
                column: "flag".to_owned(),
                value: "yes".to_owned(),
            })
        })
        .expect("projection plan is valid");
    let parallel = project(source.as_bytes(), plan.clone()).expect("parallel source must project");

    let limits = Limits {
        values: row_count * 3 + 1,
        rows: row_count + 1,
        ..Limits::default()
    };
    let serial = project_with_options(source.as_bytes(), plan, ParseOptions { limits })
        .expect("forced-serial source must project");

    assert_eq!(parallel.row_count(), retained);
    assert_eq!(parallel.row_count(), serial.row_count());
    assert_eq!(parallel.provenance(), serial.provenance());
    assert_eq!(
        parallel
            .column("id")
            .expect("parallel column exists")
            .values()
            .collect::<Vec<_>>(),
        serial
            .column("id")
            .expect("serial column exists")
            .values()
            .collect::<Vec<_>>()
    );

    let parallel_document = parse(source.as_bytes()).expect("parallel source must parse");
    let serial_document = parse_with_options(source.as_bytes(), ParseOptions { limits })
        .expect("forced-serial source must parse");
    assert_eq!(parallel_document, serial_document);
}

#[test]
fn gzip_is_detected_by_content_and_decompressed_with_limits() {
    let mut encoder = GzEncoder::new(Vec::new(), Compression::fast());
    encoder
        .write_all(MISSING_VALUES)
        .expect("in-memory compression must succeed");
    let compressed = encoder.finish().expect("gzip trailer must be written");

    let source = decode_source("no-gzip-suffix.cif", &compressed, InputLimits::default())
        .expect("gzip magic must be detected");
    assert_eq!(source.as_str().as_bytes(), MISSING_VALUES);

    let raw = decode_source("advisory.cif.gz", MISSING_VALUES, InputLimits::default())
        .expect("suffix must not force decompression");
    assert_eq!(raw.as_str().as_bytes(), MISSING_VALUES);
}

#[test]
fn gzip_failures_and_expansion_limits_are_structured() {
    let truncated = [0x1f, 0x8b, 0x08, 0x00];
    let error = decode_source("truncated", &truncated, InputLimits::default())
        .expect_err("truncated gzip must fail");
    assert_eq!(error.code(), InputErrorCode::InvalidGzip);

    let mut encoder = GzEncoder::new(Vec::new(), Compression::fast());
    encoder
        .write_all(&vec![b'x'; 16_384])
        .expect("in-memory compression must succeed");
    let compressed = encoder.finish().expect("gzip trailer must be written");
    let error = decode_source(
        "limited",
        &compressed,
        InputLimits {
            input_bytes: compressed.len(),
            decompressed_bytes: 1_024,
            decompression_ratio: 1_000,
        },
    )
    .expect_err("decompressed limit must be enforced while reading");
    assert_eq!(error.code(), InputErrorCode::DecompressedLimit);
}
