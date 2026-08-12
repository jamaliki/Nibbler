//! Source-backed ModelCIF prediction model.

use crate::pdbx::PdbxModel;

use super::prediction::{
    Data, ModelGroup, ModelGroupLink, ModelRepresentative, PredictionModel, TargetEntity,
    TargetInstance,
};
use super::provenance::{
    ArchiveMember, AssociatedFile, DataGroupMember, ProtocolStep, Software, SoftwareGroupMember,
};
use super::qa::{QaMetric, QaValue};
use super::template::{
    Alignment, AlignmentDetail, AlignmentSequence, Template, TemplateMapping, TemplateSegment,
};

/// An immutable ModelCIF prediction model backed by one generic CIF document.
#[derive(Clone, Debug, PartialEq)]
pub struct ModelCifModel {
    pub(super) coordinates: PdbxModel,
    pub(super) audit_conform: Vec<(String, String)>,
    pub(super) data: Vec<Data>,
    pub(super) data_groups: Vec<DataGroupMember>,
    pub(super) targets: Vec<TargetEntity>,
    pub(super) target_instances: Vec<TargetInstance>,
    pub(super) models: Vec<PredictionModel>,
    pub(super) model_groups: Vec<ModelGroup>,
    pub(super) model_group_links: Vec<ModelGroupLink>,
    pub(super) representatives: Vec<ModelRepresentative>,
    pub(super) software: Vec<Software>,
    pub(super) software_groups: Vec<SoftwareGroupMember>,
    pub(super) protocol_steps: Vec<ProtocolStep>,
    pub(super) templates: Vec<Template>,
    pub(super) template_segments: Vec<TemplateSegment>,
    pub(super) template_mappings: Vec<TemplateMapping>,
    pub(super) alignments: Vec<Alignment>,
    pub(super) alignment_details: Vec<AlignmentDetail>,
    pub(super) alignment_sequences: Vec<AlignmentSequence>,
    pub(super) qa_metrics: Vec<QaMetric>,
    pub(super) qa_values: Vec<QaValue>,
    pub(super) associated_files: Vec<AssociatedFile>,
    pub(super) archive_members: Vec<ArchiveMember>,
}

impl ModelCifModel {
    /// Return the shared PDBx coordinate and chemistry model.
    #[must_use]
    pub const fn coordinates(&self) -> &PdbxModel {
        &self.coordinates
    }

    /// Return the number of deposited prediction models.
    #[must_use]
    pub fn prediction_model_count(&self) -> usize {
        self.models.len()
    }

    /// Return the number of modeled target entities.
    #[must_use]
    pub fn target_entity_count(&self) -> usize {
        self.targets.len()
    }

    /// Return the number of structural templates.
    #[must_use]
    pub fn template_count(&self) -> usize {
        self.templates.len()
    }

    /// Return the number of QA metric definitions.
    #[must_use]
    pub fn qa_metric_count(&self) -> usize {
        self.qa_metrics.len()
    }

    /// Return the number of global, local, and pairwise QA values.
    #[must_use]
    pub fn qa_value_count(&self) -> usize {
        self.qa_values.len()
    }

    /// Iterate over deposited software names in source order.
    pub fn software_names(&self) -> impl ExactSizeIterator<Item = &str> {
        self.software.iter().map(|software| software.name.as_str())
    }

    /// Iterate over QA metric names in source order.
    pub fn qa_metric_names(&self) -> impl ExactSizeIterator<Item = &str> {
        self.qa_metrics.iter().map(|metric| metric.name.as_str())
    }

    /// Iterate over QA metric modes in source order.
    pub fn qa_metric_modes(&self) -> impl ExactSizeIterator<Item = &str> {
        self.qa_metrics.iter().map(|metric| metric.mode.as_str())
    }

    pub(crate) fn audit_conform(&self) -> &[(String, String)] {
        &self.audit_conform
    }
}
