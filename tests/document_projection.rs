//! Parsed-document projection shares the text and Arrow contracts.

use _core::cif::{
    CifCellRef, Predicate, ProjectionErrorCode, ProjectionPlan, SchemaName, decode_binary,
    encode_binary, parse, project, project_document,
};

fn cells(table: &_core::cif::CifTable) -> Vec<Vec<CifCellRef<'_>>> {
    table
        .columns()
        .iter()
        .map(|column| column.values().collect())
        .collect()
}

#[test]
fn document_projection_matches_text_projection() -> Result<(), Box<dyn std::error::Error>> {
    let raw = b"data_a\n_entry.id a\nloop_\n_atom_site.id\n_atom_site.Cartn_x\n_atom_site.label_alt_id\n1 1.25 ?\n2 -2 .\ndata_b\nloop_\n_atom_site.id\n_atom_site.Cartn_x\n_atom_site.label_alt_id\n3 0 'x'\n";
    let plan = ProjectionPlan::new("_ATOM_SITE")?
        .with_columns(["label_alt_id", "id", "Cartn_x"])?
        .with_schema(SchemaName::Pdbx)?;
    let from_document = project_document(&parse(raw)?, plan.clone())?;
    let from_text = project(raw, plan)?;
    assert_eq!(cells(&from_document), cells(&from_text));
    assert_eq!(from_document.row_count(), 3);
    assert_eq!(from_document.provenance()[0].source_name(), "<document>");
    assert_eq!(from_document.provenance()[0].block_code(), Some("a"));
    assert_eq!(from_document.provenance()[1].block_code(), Some("b"));
    Ok(())
}

#[test]
fn scalar_missing_values_and_save_frames_match() -> Result<(), Box<dyn std::error::Error>> {
    let raw =
        b"data_a\n_entry.id ?\n_entry.note '.'\nsave_part\n_entry.id .\n_entry.note '?'\nsave_\n";
    let plan = ProjectionPlan::new("entry")?;
    let from_document = project_document(&parse(raw)?, plan.clone())?;
    assert_eq!(cells(&from_document), cells(&project(raw, plan)?));
    assert_eq!(from_document.row_count(), 2);
    assert_eq!(from_document.provenance()[0].frame_code(), None);
    assert_eq!(from_document.provenance()[1].frame_code(), Some("part"));
    Ok(())
}

#[test]
fn predicates_use_unreturned_columns() -> Result<(), Box<dyn std::error::Error>> {
    let raw = b"data_a\nloop_\n_atom_site.id\n_atom_site.group_PDB\n1 ATOM\n2 HETATM\n3 ATOM\n";
    let plan = ProjectionPlan::new("atom_site")?
        .with_columns(["id"])?
        .with_predicate(Predicate::Equal {
            column: "group_PDB".into(),
            value: "ATOM".into(),
        })?;
    let from_document = project_document(&parse(raw)?, plan.clone())?;
    assert_eq!(cells(&from_document), cells(&project(raw, plan)?));
    assert_eq!(from_document.row_count(), 2);
    Ok(())
}

#[test]
fn binary_typed_values_project_without_text_round_trip() -> Result<(), Box<dyn std::error::Error>> {
    let raw = b"data_a\nloop_\n_atom_site.id\n_atom_site.Cartn_x\n1 1.25\n2 ?\n";
    let document = decode_binary(&encode_binary(&parse(raw)?)?)?;
    let plan = ProjectionPlan::new("atom_site")?.with_schema(SchemaName::Pdbx)?;
    assert_eq!(
        cells(&project_document(&document, plan.clone())?),
        cells(&project(raw, plan)?)
    );
    Ok(())
}

#[test]
fn absent_category_keeps_selected_column_layout() -> Result<(), Box<dyn std::error::Error>> {
    let plan = ProjectionPlan::new("atom_site")?.with_columns(["id", "Cartn_x"])?;
    let table = project_document(&parse(b"data_a\n_entry.id a\n")?, plan)?;
    assert_eq!(table.row_count(), 0);
    assert_eq!(table.columns().len(), 2);
    Ok(())
}

#[test]
fn missing_columns_and_bad_types_keep_error_codes() -> Result<(), Box<dyn std::error::Error>> {
    let document = parse(b"data_a\n_atom_site.Cartn_x nope\n")?;
    let missing = project_document(
        &document,
        ProjectionPlan::new("atom_site")?.with_columns(["id"])?,
    );
    assert!(
        matches!(missing, Err(ref error) if error.code() == ProjectionErrorCode::MissingColumn)
    );
    let bad_type = project_document(
        &document,
        ProjectionPlan::new("atom_site")?.with_schema(SchemaName::Pdbx)?,
    );
    assert!(matches!(bad_type, Err(ref error) if error.code() == ProjectionErrorCode::SchemaType));
    Ok(())
}
