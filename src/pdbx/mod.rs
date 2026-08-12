//! Immutable PDBx coordinate-model semantics over the shared generic CIF document.

pub(crate) mod category;
mod component_source;
mod error;
mod fields;
mod model;
pub(crate) mod profile;
mod registry;
mod source;
mod validation;
mod validation_atom;
mod validation_entity;
pub(crate) mod writer;

pub use error::SemanticError;
pub use model::{
    AsymUnit, AtomSite, BranchedLink, BranchedNode, ComponentAtom, ComponentBond,
    ComponentDefinition, ComponentResolution, Connection, ConnectionEndpoint, Entity, EntityKind,
    PdbxModel, PolymerMonomer,
};
pub use profile::{ProfileDiagnostic, ProfileSeverity, ProfileValidationReport};
pub use registry::ComponentRegistry;
pub use source::{build_component_registry, build_model, build_model_with_registry};
pub use validation::{validate_document, validate_model};
pub use writer::canonical_document;
