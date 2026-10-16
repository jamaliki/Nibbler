//! Immutable PDBx coordinate-model semantics over the shared generic CIF document.

mod assembly;
mod assembly_overlap;
mod chemistry;
mod component_source;
mod conformer;
mod error;
pub(crate) mod fields;
mod model;
mod registry;
mod selected_document;
mod source;
mod validation;
mod validation_atom;
mod validation_entity;
pub(crate) mod writer;

#[cfg(feature = "python")]
pub(crate) use assembly::assembly_mapping_names;

pub use assembly::{MAX_ASSEMBLY_ATOM_SITES, assembly_document, assembly_document_with_mapping};
pub use assembly_overlap::{AssemblyDeduplication, deduplicate_assembly_copies};
pub use chemistry::{
    BondSource, ChemicalEdit, PreparedChemistry, ProductAtom, ProductBond, prepare_chemistry,
};
pub(crate) use chemistry::library::registry as chemical_library_registry;
pub(crate) use chemistry::prepare_product;
pub use conformer::{
    ConformerAtom, ConformerSelection, select_atom_conformers, select_conformer,
    select_conformer_document,
};
pub use error::SemanticError;
pub use model::PdbxModel;
pub use registry::ComponentRegistry;
pub use selected_document::selected_atom_document;
pub use source::{build_component_registry, build_model, build_model_with_registry};
pub use validation::{validate_document, validate_model};
pub use writer::canonical_document;
