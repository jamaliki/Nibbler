//! Immutable PDBx coordinate-model semantics over the shared generic CIF document.

mod assembly;
mod component_source;
mod error;
pub(crate) mod fields;
mod model;
mod registry;
mod source;
mod validation;
mod validation_atom;
mod validation_entity;
pub(crate) mod writer;

pub use assembly::{MAX_ASSEMBLY_ATOM_SITES, assembly_document};
pub use error::SemanticError;
pub use model::PdbxModel;
pub use registry::ComponentRegistry;
pub use source::{build_component_registry, build_model, build_model_with_registry};
pub use validation::{validate_document, validate_model};
pub use writer::canonical_document;
