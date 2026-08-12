//! Deterministic ModelCIF category ordering.

use crate::cif::{BlockKind, CifBlock, CifDocument, CifEntry};
use crate::pdbx;

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
    reorder_document(pdbx::canonical_document(model.coordinates()))
}

fn reorder_document(document: CifDocument) -> CifDocument {
    CifDocument::new(document.blocks().iter().map(reorder_block).collect())
}

fn reorder_block(block: &CifBlock) -> CifBlock {
    let mut entries = block
        .entries()
        .iter()
        .cloned()
        .enumerate()
        .map(|(source_index, entry)| (category_rank(&entry), source_index, entry))
        .collect::<Vec<_>>();
    entries.sort_by_key(|(rank, source_index, _)| (*rank, *source_index));
    let entries = entries.into_iter().map(|(_, _, entry)| entry).collect();
    match block.kind() {
        BlockKind::Data => CifBlock::data(block.code().unwrap_or_default().to_owned(), entries),
        BlockKind::Global => CifBlock::global(entries),
    }
}

fn category_rank(entry: &CifEntry) -> usize {
    let tag = match entry {
        CifEntry::Item(item) => item.tag(),
        CifEntry::Loop(cif_loop) => cif_loop.tags().first().map_or("", String::as_str),
        CifEntry::Frame(_) => return CATEGORY_ORDER.len() + 1,
    };
    let name = category(tag).unwrap_or("");
    CATEGORY_ORDER
        .iter()
        .position(|expected| name.eq_ignore_ascii_case(expected))
        .unwrap_or(CATEGORY_ORDER.len())
}

fn category(tag: &str) -> Option<&str> {
    tag.strip_prefix('_')?.split_once('.').map(|(name, _)| name)
}
