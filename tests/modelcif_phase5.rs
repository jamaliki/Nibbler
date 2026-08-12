//! Phase 5 ModelCIF prediction-semantics integration tests.

#![allow(clippy::expect_used)]

use _core::cif::{parse, write_canonical};
use _core::modelcif::{
    MirrorPolicy, QaValue, build_model, canonical_document, canonical_document_with_mirror,
    validate_document, validate_model,
};

const PREDICTION: &[u8] = include_bytes!("fixtures/modelcif/prediction_with_qa.cif");

#[test]
fn source_builds_one_coordinate_graph_with_typed_prediction_metadata() {
    let document = parse(PREDICTION).expect("fixture syntax must parse");
    let model = build_model(&document).expect("fixture semantics must build");

    assert_eq!(model.coordinates().entry_id(), "NIBBLER_LIGAND_ION_WATER");
    assert_eq!(model.data().len(), 5);
    assert_eq!(model.data_groups().len(), 2);
    assert_eq!(model.models().len(), 1);
    assert_eq!(model.model_groups().len(), 1);
    assert_eq!(model.model_group_links().len(), 1);
    assert_eq!(model.representatives().len(), 1);
    assert_eq!(model.targets().len(), 1);
    assert_eq!(model.target_instances().len(), 1);
    assert_eq!(model.software().len(), 1);
    assert_eq!(model.software_groups().len(), 1);
    assert_eq!(model.protocol_steps().len(), 1);
    assert_eq!(model.templates().len(), 1);
    assert_eq!(model.template_segments().len(), 1);
    assert_eq!(model.template_mappings().len(), 1);
    assert_eq!(model.alignments().len(), 1);
    assert_eq!(model.alignment_details().len(), 1);
    assert_eq!(model.alignment_sequences().len(), 2);
    assert_eq!(model.qa_metrics().len(), 3);
    assert_eq!(model.qa_values().len(), 4);
    assert_eq!(model.associated_files().len(), 1);
    assert_eq!(model.archive_members().len(), 1);
}

#[test]
fn confidence_is_not_an_atom_displacement_value() {
    let document = parse(PREDICTION).expect("fixture syntax must parse");
    let model = build_model(&document).expect("fixture semantics must build");
    let b_values = model
        .coordinates()
        .atom_sites()
        .iter()
        .take(2)
        .map(|atom| atom.occupancy_and_b_iso().1)
        .collect::<Vec<_>>();
    let confidence = model
        .qa_values()
        .iter()
        .filter_map(|value| match value {
            QaValue::Local {
                metric_id: 1,
                value,
                ..
            } => Some(*value),
            _ => None,
        })
        .collect::<Vec<_>>();

    assert_eq!(b_values, [Some(10.0), Some(12.0)]);
    assert_eq!(confidence, [91.0, 73.0]);
}

#[test]
fn strict_profile_accepts_the_complete_prediction_graph() {
    let document = parse(PREDICTION).expect("fixture syntax must parse");
    let report = validate_document(&document).expect("dictionary must load");

    assert_eq!(report.dictionary_version(), "1.4.9");
    assert!(
        report
            .coverage()
            .contains(&"confidence-b-factor-separation")
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
    let mirrored = canonical_document_with_mirror(&model, MirrorPolicy::LocalMetric(1))
        .expect("pLDDT is a local metric");
    let rebuilt = build_model(&mirrored).expect("mirrored output must remain ModelCIF");
    let b_values = rebuilt
        .coordinates()
        .atom_sites()
        .iter()
        .take(2)
        .map(|atom| atom.occupancy_and_b_iso().1)
        .collect::<Vec<_>>();

    assert_eq!(b_values, [Some(91.0), Some(73.0)]);
    assert_eq!(rebuilt.qa_values(), model.qa_values());
    assert_eq!(
        model.coordinates().atom_sites()[0].occupancy_and_b_iso().1,
        Some(10.0),
        "mirroring must not mutate the source-backed model"
    );
}

#[test]
fn mirroring_rejects_a_global_metric() {
    let document = parse(PREDICTION).expect("fixture syntax must parse");
    let model = build_model(&document).expect("fixture semantics must build");
    let error = canonical_document_with_mirror(&model, MirrorPolicy::LocalMetric(2))
        .expect_err("pTM is global, not local");

    assert_eq!(error.code(), "MODELCIF_MIRROR_MODE");
}
