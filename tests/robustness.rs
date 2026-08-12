#![allow(clippy::expect_used)]

//! Bounded property and regression tests for untrusted CIF input.

use std::io::Write;

use _core::cif::{
    CifDocument, CifEntry, CifTable, InputErrorCode, InputLimits, Limits, ParseOptions, Predicate,
    ProjectionPlan, SchemaName, decode_binary, decode_source, encode_binary, exercise_lexer,
    format_text, parse, parse_serial_reference, parse_with_options, project_binary,
    project_serial_reference, project_with_options, write_canonical, write_preserving,
};
use flate2::Compression;
use flate2::write::GzEncoder;
use proptest::prelude::*;
use proptest::test_runner::FileFailurePersistence;

const MINIMAL: &[u8] = include_bytes!("fixtures/syntax/minimal.cif");
const MISSING_VALUES: &[u8] = include_bytes!("fixtures/syntax/missing_values.cif");
const MALFORMED_LOOP: &[u8] = include_bytes!("fixtures/syntax/malformed_loop.cif");
const CHEMISTRY: &[u8] = include_bytes!("fixtures/chemistry/ligand_ion_water.cif");
const GLYCAN: &[u8] = include_bytes!("fixtures/chemistry/branched_glycan.cif");
const PREDICTION: &[u8] = include_bytes!("fixtures/modelcif/prediction_with_qa.cif");

const MAX_INPUT_BYTES: usize = 256 * 1024;
const MAX_TOKEN_BYTES: usize = 64 * 1024;

#[test]
fn regression_corpus_satisfies_every_property() {
    for source in [
        MINIMAL,
        MISSING_VALUES,
        MALFORMED_LOOP,
        CHEMISTRY,
        GLYCAN,
        PREDICTION,
    ] {
        check_input(source, 9);
    }

    let source = generated_loop(1_000, 0xDEAD_BEEF);
    check_input(source.as_bytes(), 9);
    check_filtered_projection(source.as_bytes());
}

#[test]
fn semantic_documents_rebuild_from_canonical_output() {
    check_semantics(CHEMISTRY);
    check_semantics(GLYCAN);
    check_semantics(PREDICTION);
}

#[test]
fn schema_typed_projection_matches_serial_reference() {
    let plan = ProjectionPlan::new("atom_site")
        .and_then(|plan| plan.with_schema(SchemaName::Pdbx))
        .and_then(|plan| plan.with_columns(["id", "Cartn_x", "occupancy"]))
        .expect("the embedded PDBx projection is valid");
    assert_projection_results(
        project_with_options(CHEMISTRY, plan.clone(), bounded_options()),
        project_serial_reference(CHEMISTRY, plan, bounded_options()),
        true,
    );
}

proptest! {
    #![proptest_config(ProptestConfig {
        cases: 128,
        failure_persistence: Some(Box::new(FileFailurePersistence::Direct(
            "tests/robustness.proptest-regressions",
        ))),
        max_shrink_iters: 4_096,
        ..ProptestConfig::default()
    })]

    #[test]
    fn arbitrary_bytes_obey_the_bounded_contract(
        bytes in prop::collection::vec(any::<u8>(), 0..=8 * 1024),
        selector in any::<u8>(),
    ) {
        check_input(&bytes, selector);
    }

    #[test]
    fn nearby_fixture_mutations_obey_the_bounded_contract(
        fixture in 0_usize..6,
        mutations in prop::collection::vec((any::<usize>(), any::<u8>()), 1..=8),
        selector in any::<u8>(),
    ) {
        let seeds = [MINIMAL, MISSING_VALUES, MALFORMED_LOOP, CHEMISTRY, GLYCAN, PREDICTION];
        let mut bytes = seeds[fixture].to_vec();
        for (offset, byte) in mutations {
            if !bytes.is_empty() {
                let index = offset % bytes.len();
                bytes[index] = byte;
            }
        }
        check_input(&bytes, selector);
    }

    #[test]
    fn nearby_binary_mutations_remain_bounded(
        fixture in 0_usize..5,
        mutations in prop::collection::vec((any::<usize>(), any::<u8>()), 1..=8),
    ) {
        let seeds = [MINIMAL, MISSING_VALUES, CHEMISTRY, GLYCAN, PREDICTION];
        let document = parse(seeds[fixture]).expect("binary seed must parse");
        let Ok(mut bytes) = encode_binary(&document) else {
            return Ok(());
        };
        for (offset, byte) in mutations {
            let index = offset % bytes.len();
            bytes[index] = byte;
        }
        check_binary_roundtrip(&bytes);
    }

    #[test]
    fn generated_valid_loops_match_every_encoding(
        rows in 1_usize..=600,
        salt in any::<u64>(),
    ) {
        let source = generated_loop(rows, salt);
        check_input(source.as_bytes(), 9);
        check_filtered_projection(source.as_bytes());
    }

    #[test]
    fn representable_text_survives_formatting(value in any::<String>()) {
        check_formatted_text(&value);
    }
}

fn check_input(bytes: &[u8], selector: u8) {
    if bytes.len() > MAX_INPUT_BYTES {
        return;
    }
    exercise_lexer(bytes);
    check_text_roundtrip(bytes);
    check_parallel_equivalence(bytes, selected_options(selector));
    check_binary_roundtrip(bytes);
    check_compression(bytes);
    check_semantics(bytes);
}

fn check_formatted_text(text: &str) {
    if text.len() > MAX_TOKEN_BYTES {
        return;
    }
    let Ok(formatted) = format_text(text) else {
        return;
    };
    let source = if formatted.starts_with(';') {
        format!("data_property\n_property.value\n{formatted}\n")
    } else {
        format!("data_property\n_property.value {formatted}\n")
    };
    let document = parse(source.as_bytes()).expect("formatter output must parse");
    assert_eq!(scalar_text(&document), Some(text));
}

fn check_text_roundtrip(bytes: &[u8]) {
    let Ok(document) = parse_with_options(bytes, bounded_options()) else {
        return;
    };
    let canonical = write_canonical(&document).expect("parsed document must be writable");
    let reparsed = parse(canonical.as_bytes()).expect("canonical output must parse");
    assert_eq!(reparsed, document);
    assert_eq!(
        write_canonical(&reparsed).expect("reparsed document must be writable"),
        canonical
    );

    let preserving = write_preserving(&document).expect("parsed document must preserve-write");
    assert_eq!(
        parse(preserving.as_bytes()).expect("preserving output must parse"),
        document
    );
}

fn check_parallel_equivalence(bytes: &[u8], options: ParseOptions) {
    assert_eq!(
        parse_with_options(bytes, options),
        parse_serial_reference(bytes, options)
    );
    for category in ["atom_site", "property", "case", "example"] {
        let plan = ProjectionPlan::new(category).expect("static category is valid");
        assert_projection_results(
            project_with_options(bytes, plan.clone(), options),
            project_serial_reference(bytes, plan, options),
            true,
        );
    }
}

fn check_filtered_projection(bytes: &[u8]) {
    let plan = ProjectionPlan::new("case")
        .and_then(|plan| plan.with_columns(["id", "value"]))
        .and_then(|plan| {
            plan.with_predicate(Predicate::Equal {
                column: "kind".to_owned(),
                value: "keep".to_owned(),
            })
        })
        .and_then(|plan| {
            plan.with_predicate(Predicate::Missing {
                column: "value".to_owned(),
                kind: None,
            })
        })
        .expect("static projection is valid");
    assert_projection_results(
        project_with_options(bytes, plan.clone(), bounded_options()),
        project_serial_reference(bytes, plan, bounded_options()),
        true,
    );
}

fn check_binary_roundtrip(bytes: &[u8]) {
    if let Ok(document) = decode_binary(bytes) {
        check_binary_document(&document, true);
    }
    if let Ok(document) = parse_with_options(bytes, bounded_options()) {
        check_binary_document(&document, false);
    }
}

fn check_binary_document(document: &CifDocument, structurally_round_trips: bool) {
    let Ok(binary) = encode_binary(document) else {
        return;
    };
    let decoded = decode_binary(&binary).expect("encoded BinaryCIF must decode");
    if structurally_round_trips {
        assert_eq!(decoded, *document);
    }
    let Ok(text) = write_canonical(document) else {
        return;
    };
    for category in document_categories(document) {
        let Ok(plan) = ProjectionPlan::new(category) else {
            continue;
        };
        assert_cross_encoding_results(
            project_with_options(text.as_bytes(), plan.clone(), bounded_options()),
            project_binary(&binary, plan),
        );
    }
}

fn document_categories(document: &CifDocument) -> Vec<&str> {
    let mut categories = Vec::new();
    for block in document.blocks() {
        for entry in block.entries() {
            let tags: Box<dyn Iterator<Item = &str>> = match entry {
                CifEntry::Item(item) => Box::new(std::iter::once(item.tag())),
                CifEntry::Loop(loop_) => Box::new(loop_.tags().iter().map(String::as_str)),
                CifEntry::Frame(_) => continue,
            };
            for tag in tags {
                let Some((category, _)) = tag.strip_prefix('_').and_then(|tag| tag.split_once('.'))
                else {
                    continue;
                };
                if !categories.contains(&category) {
                    categories.push(category);
                }
            }
        }
    }
    categories
}

fn check_compression(bytes: &[u8]) {
    if bytes.len() > MAX_TOKEN_BYTES {
        return;
    }
    let limits = InputLimits {
        input_bytes: MAX_INPUT_BYTES,
        decompressed_bytes: MAX_INPUT_BYTES,
        decompression_ratio: 128,
    };
    let _ = decode_source("<arbitrary>", bytes, limits);

    let mut encoder = GzEncoder::new(Vec::new(), Compression::fast());
    encoder
        .write_all(bytes)
        .expect("writing to a byte vector cannot fail");
    let compressed = encoder
        .finish()
        .expect("finishing a byte vector cannot fail");
    let permissive = InputLimits {
        input_bytes: MAX_INPUT_BYTES,
        decompressed_bytes: MAX_INPUT_BYTES,
        decompression_ratio: usize::MAX,
    };
    let decoded = decode_source("<compressed>", &compressed, permissive);
    if let Ok(text) = std::str::from_utf8(bytes) {
        assert_eq!(
            decoded.as_ref().ok().map(|source| source.as_str()),
            Some(text),
            "valid UTF-8 compression round trip failed: {decoded:?}"
        );
    } else {
        assert_eq!(
            decoded.err().map(|error| error.code()),
            Some(InputErrorCode::InvalidUtf8)
        );
    }

    for end in [
        0,
        1,
        compressed.len() / 2,
        compressed.len().saturating_sub(1),
    ] {
        let _ = decode_source("<truncated>", &compressed[..end], permissive);
    }
    for index in [0, compressed.len() / 2, compressed.len().saturating_sub(1)] {
        let mut altered = compressed.clone();
        altered[index] ^= 0xFF;
        let _ = decode_source("<altered>", &altered, permissive);
    }
    if !bytes.is_empty() {
        assert_eq!(
            decode_source(
                "<limited>",
                &compressed,
                InputLimits {
                    decompressed_bytes: bytes.len() - 1,
                    ..permissive
                },
            )
            .err()
            .map(|error| error.code()),
            Some(InputErrorCode::DecompressedLimit)
        );
        assert_eq!(
            decode_source(
                "<ratio-limited>",
                &compressed,
                InputLimits {
                    decompression_ratio: 0,
                    ..permissive
                },
            )
            .err()
            .map(|error| error.code()),
            Some(InputErrorCode::DecompressionRatio)
        );
    }
}

fn check_semantics(bytes: &[u8]) {
    let Ok(document) = parse(bytes) else {
        return;
    };
    if let Ok(model) = _core::pdbx::build_model(&document) {
        let _ = _core::pdbx::validate_model(&model);
        let rebuilt_document = _core::pdbx::canonical_document(&model);
        let encoded = write_canonical(&rebuilt_document).expect("PDBx model must be writable");
        let reparsed = parse(encoded.as_bytes()).expect("PDBx output must parse");
        let rebuilt = _core::pdbx::build_model(&reparsed).expect("PDBx output must rebuild");
        assert_pdbx_summary(&model, &rebuilt);
    }
    if let Ok(model) = _core::modelcif::build_model(&document) {
        let _ = _core::modelcif::validate_model(&model);
        check_modelcif_document(&model, _core::modelcif::canonical_document(&model));
        for metric_id in 0..=3 {
            if let Ok(mirrored) =
                _core::modelcif::canonical_document_with_local_metric(&model, metric_id)
            {
                check_modelcif_document(&model, mirrored);
            }
        }
    }
}

fn assert_projection_results(
    left: Result<CifTable, _core::cif::ProjectionError>,
    right: Result<CifTable, _core::cif::ProjectionError>,
    compare_source_name: bool,
) {
    assert_eq!(
        left.is_ok(),
        right.is_ok(),
        "projection result mismatch: {left:?} {right:?}"
    );
    if let (Ok(left), Ok(right)) = (&left, &right) {
        assert_tables_equal(left, right, compare_source_name);
    }
    if let (Err(left), Err(right)) = (left, right) {
        assert_eq!(left, right);
    }
}

fn assert_cross_encoding_results<E1: std::fmt::Debug, E2: std::fmt::Debug>(
    left: Result<CifTable, E1>,
    right: Result<CifTable, E2>,
) {
    assert_eq!(
        left.is_ok(),
        right.is_ok(),
        "encoding projection mismatch: {left:?} {right:?}"
    );
    if let (Ok(left), Ok(right)) = (&left, &right) {
        assert_tables_equal(left, right, false);
    }
}

fn assert_tables_equal(left: &CifTable, right: &CifTable, compare_source_name: bool) {
    assert_eq!(left.category(), right.category());
    assert_eq!(left.row_count(), right.row_count());
    assert_eq!(left.columns().len(), right.columns().len());
    for (left, right) in left.columns().iter().zip(right.columns()) {
        assert_eq!(left.name(), right.name());
        assert!(left.values().eq(right.values()));
    }
    assert_eq!(left.provenance().len(), right.provenance().len());
    for (left, right) in left.provenance().iter().zip(right.provenance()) {
        if compare_source_name {
            assert_eq!(left.source_name(), right.source_name());
        }
        assert_eq!(left.block_code(), right.block_code());
        assert_eq!(left.frame_code(), right.frame_code());
        assert_eq!(left.row_count(), right.row_count());
    }
}

fn assert_pdbx_summary(left: &_core::pdbx::PdbxModel, right: &_core::pdbx::PdbxModel) {
    assert_eq!(left.entry_id(), right.entry_id());
    assert_eq!(left.entity_count(), right.entity_count());
    assert_eq!(left.asym_unit_count(), right.asym_unit_count());
    assert_eq!(left.component_count(), right.component_count());
    assert_eq!(left.atom_site_count(), right.atom_site_count());
    assert_eq!(left.connection_count(), right.connection_count());
    assert!(left.entity_kinds().eq(right.entity_kinds()));
    assert!(left.component_ids().eq(right.component_ids()));
}

fn check_modelcif_document(model: &_core::modelcif::ModelCifModel, document: CifDocument) {
    let encoded = write_canonical(&document).expect("ModelCIF model must be writable");
    let reparsed = parse(encoded.as_bytes()).expect("ModelCIF output must parse");
    let rebuilt = _core::modelcif::build_model(&reparsed).expect("ModelCIF output must rebuild");
    assert_eq!(
        model.coordinates().entry_id(),
        rebuilt.coordinates().entry_id()
    );
    assert_eq!(
        model.prediction_model_count(),
        rebuilt.prediction_model_count()
    );
    assert_eq!(model.target_entity_count(), rebuilt.target_entity_count());
    assert_eq!(model.template_count(), rebuilt.template_count());
    assert_eq!(model.qa_metric_count(), rebuilt.qa_metric_count());
    assert_eq!(model.qa_value_count(), rebuilt.qa_value_count());
    assert!(model.software_names().eq(rebuilt.software_names()));
    assert!(model.qa_metric_names().eq(rebuilt.qa_metric_names()));
    assert!(model.qa_metric_modes().eq(rebuilt.qa_metric_modes()));
}

fn scalar_text(document: &CifDocument) -> Option<&str> {
    let CifEntry::Item(item) = document.blocks().first()?.entries().first()? else {
        return None;
    };
    item.value().as_ref().as_text()
}

fn bounded_options() -> ParseOptions {
    ParseOptions {
        limits: Limits {
            source_bytes: MAX_INPUT_BYTES,
            token_bytes: MAX_TOKEN_BYTES,
            blocks: 256,
            frames: 1_024,
            loops: MAX_INPUT_BYTES,
            tags: MAX_INPUT_BYTES,
            loop_columns: 4_096,
            rows: MAX_INPUT_BYTES,
            values: MAX_INPUT_BYTES,
        },
    }
}

fn selected_options(selector: u8) -> ParseOptions {
    let mut options = bounded_options();
    match selector % 10 {
        0 => options.limits.token_bytes = usize::from(selector / 10) + 1,
        1 => options.limits.blocks = 0,
        2 => options.limits.frames = 0,
        3 => options.limits.loops = 0,
        4 => options.limits.tags = 0,
        5 => options.limits.loop_columns = 1,
        6 => options.limits.rows = 1,
        7 => options.limits.values = 1,
        8 | 9 => {}
        _ => {}
    }
    options
}

fn generated_loop(rows: usize, salt: u64) -> String {
    let mut source = String::with_capacity(rows * 40 + 64);
    source.push_str("data_generated\nloop_\n_case.id\n_case.kind\n_case.value\n");
    for row in 0..rows {
        let kind = if row % 2 == 0 { "keep" } else { "drop" };
        let value = match row % 11 {
            0 => ".".to_owned(),
            1 => "?".to_owned(),
            2 => "'quoted value'".to_owned(),
            _ => format!("v{:016x}", salt.rotate_left((row % 64) as u32) ^ row as u64),
        };
        source.push_str(&format!("{row} {kind} {value}\n"));
    }
    source.push_str("stop_\n");
    source
}
