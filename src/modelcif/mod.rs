//! Immutable ModelCIF prediction semantics over the shared PDBx coordinate graph.

mod aggregate;
mod fields;
mod mirror;
mod model;
mod provenance;
mod qa;
mod source;
mod source_prediction;
mod source_provenance;
mod source_qa;
mod source_template;
mod template;
mod validation;
mod validation_prediction;
mod validation_provenance;
mod validation_qa;
mod validation_template;
mod writer;

pub use aggregate::ModelCifModel;
pub use mirror::{MirrorPolicy, WriteError, canonical_document_with_mirror};
pub use model::{
    Data, ModelGroup, ModelGroupLink, ModelRepresentative, PredictionModel, TargetEntity,
    TargetInstance,
};
pub use provenance::{
    ArchiveMember, AssociatedFile, DataGroupMember, ProtocolStep, Software, SoftwareGroupMember,
};
pub use qa::{QaMetric, QaValue, ResidueSite};
pub use source::{build_model, build_model_with_registry};
pub use template::{
    Alignment, AlignmentDetail, AlignmentSequence, Template, TemplateMapping, TemplateSegment,
};
pub use validation::{validate_document, validate_model};
pub use writer::canonical_document;
