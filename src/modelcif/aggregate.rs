//! Aggregate source-backed ModelCIF model.

use crate::pdbx::PdbxModel;

use super::model::{
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
    /// Return the shared PDBx coordinate and chemistry graph.
    #[must_use]
    pub const fn coordinates(&self) -> &PdbxModel {
        &self.coordinates
    }

    /// Return the typed modeling data records.
    #[must_use]
    pub fn data(&self) -> &[Data] {
        &self.data
    }

    /// Return membership edges in protocol data groups.
    #[must_use]
    pub fn data_groups(&self) -> &[DataGroupMember] {
        &self.data_groups
    }

    /// Return modeled target entities.
    #[must_use]
    pub fn targets(&self) -> &[TargetEntity] {
        &self.targets
    }

    /// Return label-space target instances.
    #[must_use]
    pub fn target_instances(&self) -> &[TargetInstance] {
        &self.target_instances
    }

    /// Return deposited prediction models.
    #[must_use]
    pub fn models(&self) -> &[PredictionModel] {
        &self.models
    }

    /// Return named groups of prediction models.
    #[must_use]
    pub fn model_groups(&self) -> &[ModelGroup] {
        &self.model_groups
    }

    /// Return explicit model-to-group membership edges.
    #[must_use]
    pub fn model_group_links(&self) -> &[ModelGroupLink] {
        &self.model_group_links
    }

    /// Return declared representatives of model groups.
    #[must_use]
    pub fn representatives(&self) -> &[ModelRepresentative] {
        &self.representatives
    }

    /// Return software records referenced by prediction provenance.
    #[must_use]
    pub fn software(&self) -> &[Software] {
        &self.software
    }

    /// Return membership edges in ModelCIF software groups.
    #[must_use]
    pub fn software_groups(&self) -> &[SoftwareGroupMember] {
        &self.software_groups
    }

    /// Return ordered modeling-protocol steps.
    #[must_use]
    pub fn protocol_steps(&self) -> &[ProtocolStep] {
        &self.protocol_steps
    }

    /// Return structural templates used by the prediction.
    #[must_use]
    pub fn templates(&self) -> &[Template] {
        &self.templates
    }

    /// Return contiguous polymer segments declared on templates.
    #[must_use]
    pub fn template_segments(&self) -> &[TemplateSegment] {
        &self.template_segments
    }

    /// Return target-to-template segment mappings.
    #[must_use]
    pub fn template_mappings(&self) -> &[TemplateMapping] {
        &self.template_mappings
    }

    /// Return target-template alignment declarations.
    #[must_use]
    pub fn alignments(&self) -> &[Alignment] {
        &self.alignments
    }

    /// Return target-template participant and score records for alignments.
    #[must_use]
    pub fn alignment_details(&self) -> &[AlignmentDetail] {
        &self.alignment_details
    }

    /// Return target and template sequences participating in alignments.
    #[must_use]
    pub fn alignment_sequences(&self) -> &[AlignmentSequence] {
        &self.alignment_sequences
    }

    /// Return QA metric definitions.
    #[must_use]
    pub fn qa_metrics(&self) -> &[QaMetric] {
        &self.qa_metrics
    }

    /// Return global, local, and pairwise QA values.
    #[must_use]
    pub fn qa_values(&self) -> &[QaValue] {
        &self.qa_values
    }

    /// Return files associated with the ModelCIF entry.
    #[must_use]
    pub fn associated_files(&self) -> &[AssociatedFile] {
        &self.associated_files
    }

    /// Return files contained in associated archives.
    #[must_use]
    pub fn archive_members(&self) -> &[ArchiveMember] {
        &self.archive_members
    }

    pub(crate) fn audit_conform(&self) -> &[(String, String)] {
        &self.audit_conform
    }
}
