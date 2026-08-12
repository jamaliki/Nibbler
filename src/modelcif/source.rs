//! Construction of typed ModelCIF metadata over one PDBx coordinate model.

use crate::cif::CifDocument;
use crate::pdbx::category::CategoryIndex;
use crate::pdbx::{ComponentRegistry, SemanticError, build_model_with_registry as build_pdbx};

use super::aggregate::ModelCifModel;
use super::{source_prediction, source_provenance, source_qa, source_template};

/// Build ModelCIF semantics without a caller-selected component registry.
///
/// # Errors
///
/// Returns an error when the coordinate graph or a present ModelCIF row is structurally invalid.
pub fn build_model(document: &CifDocument) -> Result<ModelCifModel, SemanticError> {
    build_model_with_registry(document, None)
}

/// Build ModelCIF semantics with an optional immutable local CCD registry.
///
/// # Errors
///
/// Returns an error when the coordinate graph or a present ModelCIF row is structurally invalid.
pub fn build_model_with_registry(
    document: &CifDocument,
    registry: Option<&ComponentRegistry>,
) -> Result<ModelCifModel, SemanticError> {
    let coordinates = build_pdbx(document, registry)?;
    let categories = CategoryIndex::new(&document.blocks()[0]);

    Ok(ModelCifModel {
        coordinates,
        audit_conform: source_prediction::audit_conform(&categories)?,
        data: source_prediction::data(&categories)?,
        data_groups: source_prediction::data_groups(&categories)?,
        targets: source_prediction::targets(&categories)?,
        target_instances: source_prediction::target_instances(&categories)?,
        models: source_prediction::models(&categories)?,
        model_groups: source_prediction::model_groups(&categories)?,
        model_group_links: source_prediction::model_group_links(&categories)?,
        representatives: source_prediction::representatives(&categories)?,
        software: source_provenance::software(&categories)?,
        software_groups: source_provenance::software_groups(&categories)?,
        protocol_steps: source_provenance::protocol_steps(&categories)?,
        templates: source_template::templates(&categories)?,
        template_segments: source_template::template_segments(&categories)?,
        template_mappings: source_template::template_mappings(&categories)?,
        alignments: source_template::alignments(&categories)?,
        alignment_details: source_template::alignment_details(&categories)?,
        alignment_sequences: source_template::alignment_sequences(&categories)?,
        qa_metrics: source_qa::qa_metrics(&categories)?,
        qa_values: source_qa::qa_values(&categories)?,
        associated_files: source_provenance::associated_files(&categories)?,
        archive_members: source_provenance::archive_members(&categories)?,
    })
}
