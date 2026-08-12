//! Immutable ModelCIF prediction semantics over the shared PDBx coordinate graph.

mod aggregate;
mod fields;
mod mirror;
mod prediction;
mod provenance;
mod qa;
mod source;
mod template;
mod validation;
mod validation_prediction;
mod validation_provenance;
mod validation_qa;
mod validation_template;
mod writer;

pub use aggregate::ModelCifModel;
pub use mirror::{WriteError, canonical_document_with_local_metric};
pub use source::{build_model, build_model_with_registry};
pub use validation::{validate_document, validate_model};
pub use writer::canonical_document;
