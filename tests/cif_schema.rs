//! Pinned-schema projection and DDL2 validation integration tests.

#![allow(clippy::expect_used)]

use _core::cif::{
    CifCellRef, ProjectionErrorCode, ProjectionPlan, SchemaName, parse, project, validate_document,
};

const CHEMISTRY: &[u8] = include_bytes!("fixtures/chemistry/ligand_ion_water.cif");

#[test]
fn pinned_pdbx_schema_validates_the_qualified_fixture() {
    let document = parse(CHEMISTRY).expect("fixture must parse");
    let report = validate_document(&document, SchemaName::Pdbx).expect("schema must load");

    assert_eq!(report.schema_name(), "pdbx");
    assert_eq!(report.dictionary_version(), "5.416");
    assert!(report.coverage().contains(&"category-keys".to_owned()));
    assert!(report.is_valid(), "unexpected: {:?}", report.diagnostics());
}

#[test]
fn schema_projection_builds_numeric_arrow_storage() {
    let plan = ProjectionPlan::new("atom_site")
        .and_then(|plan| plan.with_columns(["pdbx_PDB_model_num", "Cartn_x"]))
        .and_then(|plan| plan.with_schema(SchemaName::Pdbx))
        .expect("known PDBx projection");
    let table = project(CHEMISTRY, plan).expect("fixture must project");

    assert_eq!(
        table
            .column("pdbx_PDB_model_num")
            .expect("model-number column")
            .values()
            .next(),
        Some(CifCellRef::Integer(1))
    );
    assert_eq!(
        table
            .column("Cartn_x")
            .expect("coordinate column")
            .values()
            .next(),
        Some(CifCellRef::Float(0.0))
    );
}

#[test]
fn schema_projection_rejects_bad_numeric_values() {
    let plan = ProjectionPlan::new("atom_site")
        .and_then(|plan| plan.with_columns(["Cartn_x"]))
        .and_then(|plan| plan.with_schema(SchemaName::Pdbx))
        .expect("known PDBx projection");
    let error = project(b"data_bad\nloop_\n_atom_site.Cartn_x\nnot-a-number\n", plan)
        .expect_err("typed projection must reject invalid numerics");

    assert_eq!(error.code(), ProjectionErrorCode::SchemaType);
}

#[test]
fn schema_projection_plan_order_does_not_change_typing_or_checks() {
    let plan = ProjectionPlan::new("atom_site")
        .and_then(|plan| plan.with_schema(SchemaName::Pdbx))
        .and_then(|plan| plan.with_columns(["pdbx_PDB_model_num"]))
        .expect("columns may follow schema selection");
    let table = project(CHEMISTRY, plan).expect("fixture must project");
    assert_eq!(
        table
            .column("pdbx_PDB_model_num")
            .expect("model-number column")
            .values()
            .next(),
        Some(CifCellRef::Integer(1))
    );

    let error = ProjectionPlan::new("atom_site")
        .and_then(|plan| plan.with_schema(SchemaName::Pdbx))
        .and_then(|plan| plan.with_columns(["not_a_dictionary_item"]))
        .expect_err("unknown columns must fail in either builder order");
    assert_eq!(error.code(), ProjectionErrorCode::SchemaItem);
}

#[test]
fn validator_reports_supported_constraint_families() {
    let document = parse(
        b"data_bad\n\
          loop_\n\
          _chem_comp.id\n\
          _chem_comp.formula_weight\n\
          DUP 0\n\
          DUP 0\n\
          _unknown_category.value x\n",
    )
    .expect("invalid dictionary content remains valid CIF syntax");
    let report = validate_document(&document, SchemaName::Pdbx).expect("schema must load");
    let codes = report
        .diagnostics()
        .iter()
        .map(|diagnostic| diagnostic.code())
        .collect::<Vec<_>>();

    assert!(!report.is_valid());
    assert!(codes.contains(&"CIF_SCHEMA_RANGE"));
    assert!(codes.contains(&"CIF_SCHEMA_KEY_DUPLICATE"));
    assert!(codes.contains(&"CIF_SCHEMA_ITEM_UNKNOWN"));
}

#[test]
fn validator_requires_a_declared_parent_category() {
    let document = parse(b"data_bad\n_atom_site.type_symbol C\n")
        .expect("missing dictionary parents remain valid CIF syntax");
    let report = validate_document(&document, SchemaName::Pdbx).expect("schema must load");

    assert!(
        report
            .diagnostics()
            .iter()
            .any(|diagnostic| diagnostic.code() == "CIF_SCHEMA_PARENT_MISSING")
    );
}

#[test]
fn validator_compares_uchar_keys_and_links_case_insensitively() {
    let duplicate = parse(b"data_bad\nloop_\n_chem_comp.id\nala\nALA\n")
        .expect("case variants remain valid CIF syntax");
    let duplicate_report =
        validate_document(&duplicate, SchemaName::Pdbx).expect("schema must load");
    assert!(
        duplicate_report
            .diagnostics()
            .iter()
            .any(|diagnostic| diagnostic.code() == "CIF_SCHEMA_KEY_DUPLICATE")
    );

    let linked = parse(b"data_linked\n_chem_comp.id zn\n_atom_site.label_comp_id ZN\n")
        .expect("case variants remain valid CIF syntax");
    let linked_report = validate_document(&linked, SchemaName::Pdbx).expect("schema must load");
    assert!(!linked_report.diagnostics().iter().any(|diagnostic| {
        diagnostic.code() == "CIF_SCHEMA_PARENT_MISSING"
            && diagnostic
                .context()
                .iter()
                .any(|value| value == "item=_atom_site.label_comp_id")
    }));
}
