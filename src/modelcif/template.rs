//! Template and alignment records and row decoding.

use crate::cif::CategoryView;
use crate::pdbx::SemanticError;
use crate::pdbx::fields::category_rows;

use super::fields::{
    optional_float, optional_integer, optional_text, required_integer, required_text,
};

#[derive(Clone, Debug, Eq, PartialEq)]
pub(super) struct Template {
    pub ordinal_id: i64,
    pub id: i64,
    pub origin: String,
    pub entity_type: String,
    pub data_id: i64,
    pub target_asym_id: String,
    pub auth_asym_id: String,
    pub label_asym_id: Option<String>,
    pub label_entity_id: Option<String>,
    pub model_number: i64,
    pub transform_id: i64,
    pub name: Option<String>,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(super) struct TemplateSegment {
    pub id: i64,
    pub template_id: i64,
    pub sequence_begin: Option<i64>,
    pub sequence_end: Option<i64>,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(super) struct TemplateMapping {
    pub id: i64,
    pub template_segment_id: i64,
    pub target_asym_id: String,
    pub target_begin: Option<i64>,
    pub target_end: Option<i64>,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(super) struct Alignment {
    pub id: i64,
    pub data_id: i64,
    pub software_group_id: Option<i64>,
    pub length: Option<i64>,
    pub alignment_type: Option<String>,
    pub mode: Option<String>,
}

#[derive(Clone, Debug, PartialEq)]
pub(super) struct AlignmentDetail {
    pub ordinal_id: i64,
    pub alignment_id: i64,
    pub template_segment_id: i64,
    pub target_asym_id: String,
    pub score_type: Option<String>,
    pub score_value: Option<f64>,
    pub sequence_identity: Option<f64>,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(super) struct AlignmentSequence {
    pub ordinal_id: i64,
    pub alignment_id: i64,
    pub target_template_flag: String,
    pub sequence: String,
}

pub(super) fn templates(categories: CategoryView<'_>) -> Result<Vec<Template>, SemanticError> {
    category_rows(
        categories,
        "ma_template_details",
        [
            "ordinal_id",
            "template_id",
            "template_origin",
            "template_entity_type",
            "template_data_id",
            "target_asym_id",
            "template_auth_asym_id",
            "template_label_asym_id",
            "template_label_entity_id",
            "template_model_num",
            "template_trans_matrix_id",
            "template_name",
        ],
    )
    .map(|row| {
        Ok(Template {
            ordinal_id: required_integer(&row, "ma_template_details", "ordinal_id", 0)?,
            id: required_integer(&row, "ma_template_details", "template_id", 1)?,
            origin: required_text(&row, "ma_template_details", "template_origin", 2)?,
            entity_type: required_text(&row, "ma_template_details", "template_entity_type", 3)?,
            data_id: required_integer(&row, "ma_template_details", "template_data_id", 4)?,
            target_asym_id: required_text(&row, "ma_template_details", "target_asym_id", 5)?,
            auth_asym_id: required_text(&row, "ma_template_details", "template_auth_asym_id", 6)?,
            label_asym_id: optional_text(&row, 7),
            label_entity_id: optional_text(&row, 8),
            model_number: required_integer(&row, "ma_template_details", "template_model_num", 9)?,
            transform_id: required_integer(
                &row,
                "ma_template_details",
                "template_trans_matrix_id",
                10,
            )?,
            name: optional_text(&row, 11),
        })
    })
    .collect()
}

pub(super) fn template_segments(
    categories: CategoryView<'_>,
) -> Result<Vec<TemplateSegment>, SemanticError> {
    category_rows(
        categories,
        "ma_template_poly_segment",
        [
            "id",
            "template_id",
            "residue_number_begin",
            "residue_number_end",
        ],
    )
    .map(|row| {
        Ok(TemplateSegment {
            id: required_integer(&row, "ma_template_poly_segment", "id", 0)?,
            template_id: required_integer(&row, "ma_template_poly_segment", "template_id", 1)?,
            sequence_begin: optional_integer(
                &row,
                "ma_template_poly_segment",
                "residue_number_begin",
                2,
            )?,
            sequence_end: optional_integer(
                &row,
                "ma_template_poly_segment",
                "residue_number_end",
                3,
            )?,
        })
    })
    .collect()
}

pub(super) fn template_mappings(
    categories: CategoryView<'_>,
) -> Result<Vec<TemplateMapping>, SemanticError> {
    category_rows(
        categories,
        "ma_target_template_poly_mapping",
        [
            "id",
            "template_segment_id",
            "target_asym_id",
            "target_seq_id_begin",
            "target_seq_id_end",
        ],
    )
    .map(|row| {
        Ok(TemplateMapping {
            id: required_integer(&row, "ma_target_template_poly_mapping", "id", 0)?,
            template_segment_id: required_integer(
                &row,
                "ma_target_template_poly_mapping",
                "template_segment_id",
                1,
            )?,
            target_asym_id: required_text(
                &row,
                "ma_target_template_poly_mapping",
                "target_asym_id",
                2,
            )?,
            target_begin: optional_integer(
                &row,
                "ma_target_template_poly_mapping",
                "target_seq_id_begin",
                3,
            )?,
            target_end: optional_integer(
                &row,
                "ma_target_template_poly_mapping",
                "target_seq_id_end",
                4,
            )?,
        })
    })
    .collect()
}

pub(super) fn alignments(categories: CategoryView<'_>) -> Result<Vec<Alignment>, SemanticError> {
    category_rows(
        categories,
        "ma_alignment_info",
        [
            "alignment_id",
            "data_id",
            "software_group_id",
            "alignment_length",
            "alignment_type",
            "alignment_mode",
        ],
    )
    .map(|row| {
        Ok(Alignment {
            id: required_integer(&row, "ma_alignment_info", "alignment_id", 0)?,
            data_id: required_integer(&row, "ma_alignment_info", "data_id", 1)?,
            software_group_id: optional_integer(&row, "ma_alignment_info", "software_group_id", 2)?,
            length: optional_integer(&row, "ma_alignment_info", "alignment_length", 3)?,
            alignment_type: optional_text(&row, 4),
            mode: optional_text(&row, 5),
        })
    })
    .collect()
}

pub(super) fn alignment_details(
    categories: CategoryView<'_>,
) -> Result<Vec<AlignmentDetail>, SemanticError> {
    category_rows(
        categories,
        "ma_alignment_details",
        [
            "ordinal_id",
            "alignment_id",
            "template_segment_id",
            "target_asym_id",
            "score_type",
            "score_value",
            "percent_sequence_identity",
        ],
    )
    .map(|row| {
        Ok(AlignmentDetail {
            ordinal_id: required_integer(&row, "ma_alignment_details", "ordinal_id", 0)?,
            alignment_id: required_integer(&row, "ma_alignment_details", "alignment_id", 1)?,
            template_segment_id: required_integer(
                &row,
                "ma_alignment_details",
                "template_segment_id",
                2,
            )?,
            target_asym_id: required_text(&row, "ma_alignment_details", "target_asym_id", 3)?,
            score_type: optional_text(&row, 4),
            score_value: optional_float(&row, "ma_alignment_details", "score_value", 5)?,
            sequence_identity: optional_float(
                &row,
                "ma_alignment_details",
                "percent_sequence_identity",
                6,
            )?,
        })
    })
    .collect()
}

pub(super) fn alignment_sequences(
    categories: CategoryView<'_>,
) -> Result<Vec<AlignmentSequence>, SemanticError> {
    category_rows(
        categories,
        "ma_alignment",
        [
            "ordinal_id",
            "alignment_id",
            "target_template_flag",
            "sequence",
        ],
    )
    .map(|row| {
        Ok(AlignmentSequence {
            ordinal_id: required_integer(&row, "ma_alignment", "ordinal_id", 0)?,
            alignment_id: required_integer(&row, "ma_alignment", "alignment_id", 1)?,
            target_template_flag: required_text(&row, "ma_alignment", "target_template_flag", 2)?,
            sequence: required_text(&row, "ma_alignment", "sequence", 3)?,
        })
    })
    .collect()
}
