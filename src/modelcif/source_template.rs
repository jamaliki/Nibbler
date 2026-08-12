//! Structural-template and alignment row decoding.

use crate::pdbx::SemanticError;
use crate::pdbx::category::CategoryIndex;

use super::fields::{
    optional_float, optional_integer, optional_text, required_integer, required_text,
};
use super::template::{
    Alignment, AlignmentDetail, AlignmentSequence, Template, TemplateMapping, TemplateSegment,
};

pub(super) fn templates(categories: &CategoryIndex<'_>) -> Result<Vec<Template>, SemanticError> {
    categories
        .rows("ma_template_details")
        .iter()
        .enumerate()
        .map(|(row_index, row)| {
            Ok(Template {
                ordinal_id: required_integer(row, "ma_template_details", "ordinal_id", row_index)?,
                id: required_integer(row, "ma_template_details", "template_id", row_index)?,
                origin: required_text(row, "ma_template_details", "template_origin", row_index)?,
                entity_type: required_text(
                    row,
                    "ma_template_details",
                    "template_entity_type",
                    row_index,
                )?,
                data_id: required_integer(
                    row,
                    "ma_template_details",
                    "template_data_id",
                    row_index,
                )?,
                target_asym_id: required_text(
                    row,
                    "ma_template_details",
                    "target_asym_id",
                    row_index,
                )?,
                auth_asym_id: required_text(
                    row,
                    "ma_template_details",
                    "template_auth_asym_id",
                    row_index,
                )?,
                label_asym_id: optional_text(row, "template_label_asym_id"),
                label_entity_id: optional_text(row, "template_label_entity_id"),
                model_number: required_integer(
                    row,
                    "ma_template_details",
                    "template_model_num",
                    row_index,
                )?,
                transform_id: required_integer(
                    row,
                    "ma_template_details",
                    "template_trans_matrix_id",
                    row_index,
                )?,
                name: optional_text(row, "template_name"),
            })
        })
        .collect()
}

pub(super) fn template_segments(
    categories: &CategoryIndex<'_>,
) -> Result<Vec<TemplateSegment>, SemanticError> {
    categories
        .rows("ma_template_poly_segment")
        .iter()
        .enumerate()
        .map(|(row_index, row)| {
            Ok(TemplateSegment {
                id: required_integer(row, "ma_template_poly_segment", "id", row_index)?,
                template_id: required_integer(
                    row,
                    "ma_template_poly_segment",
                    "template_id",
                    row_index,
                )?,
                sequence_begin: optional_integer(
                    row,
                    "ma_template_poly_segment",
                    "residue_number_begin",
                    row_index,
                )?,
                sequence_end: optional_integer(
                    row,
                    "ma_template_poly_segment",
                    "residue_number_end",
                    row_index,
                )?,
            })
        })
        .collect()
}

pub(super) fn template_mappings(
    categories: &CategoryIndex<'_>,
) -> Result<Vec<TemplateMapping>, SemanticError> {
    categories
        .rows("ma_target_template_poly_mapping")
        .iter()
        .enumerate()
        .map(|(row_index, row)| {
            Ok(TemplateMapping {
                id: required_integer(row, "ma_target_template_poly_mapping", "id", row_index)?,
                template_segment_id: required_integer(
                    row,
                    "ma_target_template_poly_mapping",
                    "template_segment_id",
                    row_index,
                )?,
                target_asym_id: required_text(
                    row,
                    "ma_target_template_poly_mapping",
                    "target_asym_id",
                    row_index,
                )?,
                target_begin: optional_integer(
                    row,
                    "ma_target_template_poly_mapping",
                    "target_seq_id_begin",
                    row_index,
                )?,
                target_end: optional_integer(
                    row,
                    "ma_target_template_poly_mapping",
                    "target_seq_id_end",
                    row_index,
                )?,
            })
        })
        .collect()
}

pub(super) fn alignments(categories: &CategoryIndex<'_>) -> Result<Vec<Alignment>, SemanticError> {
    categories
        .rows("ma_alignment_info")
        .iter()
        .enumerate()
        .map(|(row_index, row)| {
            Ok(Alignment {
                id: required_integer(row, "ma_alignment_info", "alignment_id", row_index)?,
                data_id: required_integer(row, "ma_alignment_info", "data_id", row_index)?,
                software_group_id: optional_integer(
                    row,
                    "ma_alignment_info",
                    "software_group_id",
                    row_index,
                )?,
                length: optional_integer(row, "ma_alignment_info", "alignment_length", row_index)?,
                alignment_type: optional_text(row, "alignment_type"),
                mode: optional_text(row, "alignment_mode"),
            })
        })
        .collect()
}

pub(super) fn alignment_details(
    categories: &CategoryIndex<'_>,
) -> Result<Vec<AlignmentDetail>, SemanticError> {
    categories
        .rows("ma_alignment_details")
        .iter()
        .enumerate()
        .map(|(row_index, row)| {
            Ok(AlignmentDetail {
                ordinal_id: required_integer(row, "ma_alignment_details", "ordinal_id", row_index)?,
                alignment_id: required_integer(
                    row,
                    "ma_alignment_details",
                    "alignment_id",
                    row_index,
                )?,
                template_segment_id: required_integer(
                    row,
                    "ma_alignment_details",
                    "template_segment_id",
                    row_index,
                )?,
                target_asym_id: required_text(
                    row,
                    "ma_alignment_details",
                    "target_asym_id",
                    row_index,
                )?,
                score_type: optional_text(row, "score_type"),
                score_value: optional_float(row, "ma_alignment_details", "score_value", row_index)?,
                sequence_identity: optional_float(
                    row,
                    "ma_alignment_details",
                    "percent_sequence_identity",
                    row_index,
                )?,
            })
        })
        .collect()
}

pub(super) fn alignment_sequences(
    categories: &CategoryIndex<'_>,
) -> Result<Vec<AlignmentSequence>, SemanticError> {
    categories
        .rows("ma_alignment")
        .iter()
        .enumerate()
        .map(|(row_index, row)| {
            Ok(AlignmentSequence {
                ordinal_id: required_integer(row, "ma_alignment", "ordinal_id", row_index)?,
                alignment_id: required_integer(row, "ma_alignment", "alignment_id", row_index)?,
                target_template_flag: required_text(
                    row,
                    "ma_alignment",
                    "target_template_flag",
                    row_index,
                )?,
                sequence: required_text(row, "ma_alignment", "sequence", row_index)?,
            })
        })
        .collect()
}
