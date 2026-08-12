//! BinaryCIF and preserving-writer acceptance tests for Phase 6.

#![allow(clippy::expect_used)]

use _core::cif::{
    BinaryCifErrorCode, CifCellRef, ProjectionPlan, decode_binary, encode_binary, parse, project,
    project_binary, write_preserving,
};

const SOURCE: &[u8] = b"data_binary\n\
_entry.id binary\n\
loop_\n\
_atom_site.label_comp_id\n\
_atom_site.Cartn_x\n\
ALA 1.25\n\
? .\n";

#[test]
fn binary_projection_matches_text_and_preserves_missing_kinds() {
    let document = parse(SOURCE).expect("fixture must parse");
    let binary = encode_binary(&document).expect("fixture must encode");
    let plan = || {
        ProjectionPlan::new("atom_site")
            .and_then(|plan| plan.with_columns(["label_comp_id", "Cartn_x"]))
            .expect("projection must compile")
    };
    let text = project(SOURCE, plan()).expect("text fixture must project");
    let decoded = project_binary(&binary, plan()).expect("binary fixture must project");

    for column_name in ["label_comp_id", "Cartn_x"] {
        let text_values = text
            .column(column_name)
            .expect("text column")
            .values()
            .collect::<Vec<_>>();
        let binary_values = decoded
            .column(column_name)
            .expect("binary column")
            .values()
            .collect::<Vec<_>>();
        assert_eq!(binary_values, text_values);
    }
    assert_eq!(
        decoded
            .column("label_comp_id")
            .expect("label column")
            .values()
            .nth(1),
        Some(CifCellRef::Unknown)
    );
    assert_eq!(
        decoded
            .column("Cartn_x")
            .expect("coordinate column")
            .values()
            .nth(1),
        Some(CifCellRef::NotApplicable)
    );
}

#[test]
fn binary_encoding_is_deterministic_and_decodes_to_the_shared_model() {
    let document = parse(SOURCE).expect("fixture must parse");
    let first = encode_binary(&document).expect("fixture must encode");
    let second = encode_binary(&document).expect("fixture must encode repeatedly");
    assert_eq!(first, second);

    let decoded = decode_binary(&first).expect("encoded fixture must decode");
    assert_eq!(decoded.blocks().len(), 1);
    assert_eq!(decoded.blocks()[0].code(), Some("binary"));
}

#[test]
fn binary_container_failures_are_structured() {
    let error = decode_binary(&[0x83]).expect_err("truncated MessagePack must fail");
    assert_eq!(error.code(), BinaryCifErrorCode::Container);
}

#[test]
fn preserving_output_reuses_source_quote_styles_and_reparses() {
    let source = b"data_quotes\n_a.one 'alpha beta'\n_a.two \"gamma\"\n";
    let document = parse(source).expect("fixture must parse");
    let output = write_preserving(&document).expect("fixture must serialize");

    assert!(output.contains("_a.one 'alpha beta'\n"));
    assert!(output.contains("_a.two \"gamma\"\n"));
    assert_eq!(parse(output.as_bytes()), Ok(document));
}
