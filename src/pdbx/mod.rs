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
mod sasa;
mod selected_document;
mod short_text;
mod source;
mod validation;
mod validation_atom;
mod validation_entity;
pub(crate) mod writer;

#[cfg(feature = "python")]
pub(crate) use assembly::{assembly_mapping_names, assembly_with_model};
#[cfg(feature = "python")]
pub(crate) use model::{AtomSite, only_atom_site_loop};
#[cfg(feature = "python")]
pub(crate) use model::{SiteArray, atom_site_arrays};

pub use assembly::{MAX_ASSEMBLY_ATOM_SITES, assembly_document, assembly_document_with_mapping};
pub use assembly_overlap::{AssemblyDeduplication, deduplicate_assembly_copies};
pub(crate) use chemistry::library::registry as chemical_library_registry;
pub(crate) use chemistry::{hydrogen as hydrogen_element, prepare_product, product_rows_in, product_rows_with};
pub use short_text::ShortText;
pub use chemistry::{
    BondSource, ChemicalEdit, PreparedChemistry, ProductAtom, ProductBond, prepare_chemistry,
};
pub use conformer::{
    ConformerAtom, ConformerPolicy, ConformerSelection, select_atom_conformers,
    select_atom_conformers_with_policy, select_conformer, select_conformer_document,
    select_conformer_document_with_policy, select_conformer_with_policy,
};
pub use error::SemanticError;
pub use model::PdbxModel;
pub use registry::ComponentRegistry;
pub use sasa::{ligand_sasa, model_ligand_sasa};
pub use selected_document::selected_atom_document;
pub use source::{build_component_registry, build_model, build_model_with_registry};
#[cfg(feature = "python")]
pub(crate) use source::build_model_with_copies;
pub use validation::{validate_document, validate_model};
pub use writer::canonical_document;
