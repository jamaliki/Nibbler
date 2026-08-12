//! Construction of typed ModelCIF metadata over one PDBx coordinate model.

use std::error::Error;
use std::fmt::{self, Display, Formatter};

use crate::cif::CifDocument;
use crate::pdbx::category::CategoryIndex;
use crate::pdbx::{ComponentRegistry, build_model_with_registry as build_pdbx};

use super::aggregate::ModelCifModel;
use super::{source_prediction, source_provenance, source_qa, source_template};

/// A structural or typed failure while building a ModelCIF semantic model.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ModelCifError {
    code: &'static str,
    message: String,
    context: Vec<String>,
}

impl ModelCifError {
    pub(super) fn new(
        code: &'static str,
        message: impl Into<String>,
        context: Vec<String>,
    ) -> Self {
        Self {
            code,
            message: message.into(),
            context,
        }
    }

    /// Return the stable machine-readable error code.
    #[must_use]
    pub const fn code(&self) -> &'static str {
        self.code
    }

    /// Return the human-readable failure detail.
    #[must_use]
    pub fn message(&self) -> &str {
        &self.message
    }

    /// Return category, item, and row context when available.
    #[must_use]
    pub fn context(&self) -> &[String] {
        &self.context
    }
}

impl Display for ModelCifError {
    fn fmt(&self, formatter: &mut Formatter<'_>) -> fmt::Result {
        write!(formatter, "{}: {}", self.code, self.message)
    }
}

impl Error for ModelCifError {}

/// Build ModelCIF semantics without a caller-selected component registry.
///
/// # Errors
///
/// Returns an error when the coordinate graph or a present ModelCIF row is structurally invalid.
pub fn build_model(document: &CifDocument) -> Result<ModelCifModel, ModelCifError> {
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
) -> Result<ModelCifModel, ModelCifError> {
    let coordinates = build_pdbx(document, registry).map_err(|error| {
        ModelCifError::new(error.code(), error.message(), error.context().to_vec())
    })?;
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
