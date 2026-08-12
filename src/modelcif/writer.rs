//! Deterministic ModelCIF category ordering.

use crate::cif::CifDocument;
use crate::pdbx::writer::canonical_document_ordered;

use super::aggregate::ModelCifModel;

const CATEGORY_ORDER: &[&str] = &[
    "entry",
    "audit_conform",
    "entity",
    "entity_poly",
    "entity_poly_seq",
    "pdbx_entity_nonpoly",
    "pdbx_entity_branch",
    "pdbx_entity_branch_list",
    "pdbx_entity_branch_link",
    "chem_comp",
    "chem_comp_atom",
    "chem_comp_bond",
    "struct_asym",
    "ma_target_entity",
    "ma_target_entity_instance",
    "ma_struct_assembly_details",
    "ma_struct_assembly",
    "pdbx_poly_seq_scheme",
    "pdbx_nonpoly_scheme",
    "pdbx_branch_scheme",
    "software",
    "ma_software_parameter",
    "ma_software_group",
    "ma_data",
    "ma_data_group",
    "ma_protocol_step",
    "ma_model_group",
    "ma_model_list",
    "ma_model_group_link",
    "ma_model_representative",
    "ma_template_trans_matrix",
    "ma_template_details",
    "ma_template_poly",
    "ma_template_non_poly",
    "ma_template_branched",
    "ma_template_poly_segment",
    "ma_target_template_poly_mapping",
    "ma_alignment_info",
    "ma_alignment_details",
    "ma_alignment",
    "atom_type",
    "atom_site",
    "struct_conn_type",
    "struct_conn",
    "ma_qa_metric",
    "ma_qa_metric_global",
    "ma_qa_metric_local",
    "ma_qa_metric_local_pairwise",
    "ma_entry_associated_files",
    "ma_associated_archive_file_details",
];

/// Construct canonical ModelCIF category order without confidence mirroring.
#[must_use]
pub fn canonical_document(model: &ModelCifModel) -> CifDocument {
    canonical_document_ordered(model.coordinates(), CATEGORY_ORDER)
}
