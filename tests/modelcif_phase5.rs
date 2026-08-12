//! Phase 5 ModelCIF prediction-semantics integration tests.

#![allow(clippy::expect_used)]

use _core::cif::{CifDocument, CifEntry, CifValueRef, parse, write_canonical};
use _core::modelcif::{
    build_model, canonical_document, canonical_document_with_local_metric, validate_document,
    validate_model,
};

const PREDICTION: &[u8] = include_bytes!("fixtures/modelcif/prediction_with_qa.cif");

#[test]
fn source_builds_one_coordinate_graph_with_typed_prediction_metadata() {
    let document = parse(PREDICTION).expect("fixture syntax must parse");
    let model = build_model(&document).expect("fixture semantics must build");

    assert_eq!(model.coordinates().entry_id(), "NIBBLER_LIGAND_ION_WATER");
    assert_eq!(model.prediction_model_count(), 1);
    assert_eq!(model.target_entity_count(), 1);
    assert_eq!(model.template_count(), 1);
    assert_eq!(model.qa_metric_count(), 3);
    assert_eq!(model.qa_value_count(), 4);
    assert_eq!(model.software_names().collect::<Vec<_>>(), ["Rosetta"]);
    assert_eq!(
        model.qa_metric_names().collect::<Vec<_>>(),
        ["pLDDT", "pTM", "PAE"]
    );
    assert_eq!(
        model.qa_metric_modes().collect::<Vec<_>>(),
        ["local", "global", "local-pairwise"]
    );
    assert!(validate_model(&model).is_valid());
}

#[test]
fn confidence_is_not_an_atom_displacement_value() {
    let document = parse(PREDICTION).expect("fixture syntax must parse");
    let model = build_model(&document).expect("fixture semantics must build");
    let canonical = canonical_document(&model);
    let b_values = loop_floats(&canonical, "_atom_site.B_iso_or_equiv");
    let confidence = loop_floats(&canonical, "_ma_qa_metric_local.metric_value");

    assert_eq!(&b_values[..2], [Some(10.0), Some(12.0)]);
    assert_eq!(confidence, [Some(91.0), Some(73.0)]);
}

#[test]
fn strict_profile_accepts_the_complete_prediction_graph() {
    let document = parse(PREDICTION).expect("fixture syntax must parse");
    let report = validate_document(&document).expect("dictionary must load");

    assert_eq!(report.schema_name(), "modelcif");
    assert_eq!(report.dictionary_version(), "1.4.9");
    assert!(
        report
            .coverage()
            .iter()
            .any(|value| value == "confidence-b-factor-separation")
    );
    assert!(report.is_valid(), "unexpected: {:?}", report.diagnostics());
}

#[test]
fn strict_profile_rejects_missing_software_version() {
    let source = String::from_utf8(PREDICTION.to_vec())
        .expect("fixture is UTF-8")
        .replace(
            "1 Rosetta 'model building' 2.1.0",
            "1 Rosetta 'model building' .",
        );
    let document = parse(source.as_bytes()).expect("modified fixture syntax must parse");
    let model =
        build_model(&document).expect("missing optional source field remains representable");

    assert!(
        validate_model(&model)
            .diagnostics()
            .iter()
            .any(|diagnostic| diagnostic.code() == "MODELCIF_SOFTWARE_VERSION")
    );
}

#[test]
fn strict_profile_rejects_dangling_qa_sites_and_data_flow() {
    let source = String::from_utf8(PREDICTION.to_vec())
        .expect("fixture is UTF-8")
        .replace("2 1 A 2 MSE 1 73.0", "2 1 A 99 MSE 1 73.0")
        .replace(
            "1 1 1 modeling prediction 'coordinates generated from supplied target sequence' 1 1 2",
            "1 1 1 modeling prediction 'coordinates generated from supplied target sequence' 1 99 2",
        );
    let document = parse(source.as_bytes()).expect("modified fixture syntax must parse");
    let model = build_model(&document).expect("dangling references remain inspectable");
    let report = validate_model(&model);
    let codes = report
        .diagnostics()
        .iter()
        .map(|diagnostic| diagnostic.code())
        .collect::<Vec<_>>();

    assert!(codes.contains(&"MODELCIF_PROTOCOL_INPUT"));
    assert!(codes.contains(&"MODELCIF_QA_SITE"));
}

#[test]
fn canonical_writer_preserves_qa_and_uses_profile_order() {
    let document = parse(PREDICTION).expect("fixture syntax must parse");
    let model = build_model(&document).expect("fixture semantics must build");
    let canonical = canonical_document(&model);
    let output = write_canonical(&canonical).expect("canonical document must serialize");
    let reparsed = parse(output.as_bytes()).expect("canonical output must reparse");

    assert_eq!(reparsed, canonical);
    assert!(
        output.find("_ma_target_entity.entity_id").expect("target")
            < output.find("_ma_model_list.ordinal_id").expect("model")
    );
    assert!(
        output.find("_atom_site.group_PDB").expect("coordinates")
            < output.find("_ma_qa_metric.id").expect("QA definitions")
    );
    assert!(
        validate_document(&reparsed)
            .expect("dictionary must load")
            .is_valid()
    );
}

#[test]
fn explicit_local_metric_mirroring_changes_only_the_output_view() {
    let document = parse(PREDICTION).expect("fixture syntax must parse");
    let model = build_model(&document).expect("fixture semantics must build");
    let mirrored =
        canonical_document_with_local_metric(&model, 1).expect("pLDDT is a local metric");
    let rebuilt = build_model(&mirrored).expect("mirrored output must remain ModelCIF");
    let b_values = loop_floats(&mirrored, "_atom_site.B_iso_or_equiv");

    assert_eq!(&b_values[..2], [Some(91.0), Some(73.0)]);
    assert_eq!(
        loop_floats(&mirrored, "_ma_qa_metric_local.metric_value"),
        loop_floats(
            &canonical_document(&model),
            "_ma_qa_metric_local.metric_value"
        ),
        "mirroring must preserve ModelCIF confidence values"
    );
    assert_eq!(
        &loop_floats(&canonical_document(&model), "_atom_site.B_iso_or_equiv")[..2],
        [Some(10.0), Some(12.0)],
        "mirroring must not mutate the source-backed model"
    );
    assert!(validate_model(&rebuilt).is_valid());
}

#[test]
fn mirroring_rejects_a_global_metric() {
    let document = parse(PREDICTION).expect("fixture syntax must parse");
    let model = build_model(&document).expect("fixture semantics must build");
    let error =
        canonical_document_with_local_metric(&model, 2).expect_err("pTM is global, not local");

    assert_eq!(error.code(), "MODELCIF_MIRROR_MODE");
}

fn loop_floats(document: &CifDocument, tag: &str) -> Vec<Option<f64>> {
    document
        .blocks()
        .iter()
        .flat_map(|block| block.entries())
        .find_map(|entry| {
            let CifEntry::Loop(cif_loop) = entry else {
                return None;
            };
            let column = cif_loop
                .tags()
                .iter()
                .position(|candidate| candidate.eq_ignore_ascii_case(tag))?;
            Some(
                (0..cif_loop.row_count())
                    .map(|row| match cif_loop.value(row, column) {
                        Some(CifValueRef::Float(value, ..)) => Some(value),
                        Some(CifValueRef::Integer(value, _)) => Some(value as f64),
                        Some(CifValueRef::Text(value)) => value.as_str().parse().ok(),
                        _ => None,
                    })
                    .collect(),
            )
        })
        .unwrap_or_default()
}
