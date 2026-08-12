//! Target, data, and prediction-model row decoding.

use crate::pdbx::SemanticError;
use crate::pdbx::category::CategoryIndex;

use super::fields::{optional_integer, optional_text, required_integer, required_text};
use super::model::{
    Data, ModelGroup, ModelGroupLink, ModelRepresentative, PredictionModel, TargetEntity,
    TargetInstance,
};
use super::provenance::DataGroupMember;

pub(super) fn audit_conform(
    categories: &CategoryIndex<'_>,
) -> Result<Vec<(String, String)>, SemanticError> {
    categories
        .rows("audit_conform")
        .iter()
        .enumerate()
        .map(|(row_index, row)| {
            Ok((
                required_text(row, "audit_conform", "dict_name", row_index)?,
                required_text(row, "audit_conform", "dict_version", row_index)?,
            ))
        })
        .collect()
}

pub(super) fn data(categories: &CategoryIndex<'_>) -> Result<Vec<Data>, SemanticError> {
    categories
        .rows("ma_data")
        .iter()
        .enumerate()
        .map(|(row_index, row)| {
            Ok(Data {
                id: required_integer(row, "ma_data", "id", row_index)?,
                name: required_text(row, "ma_data", "name", row_index)?,
                content_type: required_text(row, "ma_data", "content_type", row_index)?,
            })
        })
        .collect()
}

pub(super) fn data_groups(
    categories: &CategoryIndex<'_>,
) -> Result<Vec<DataGroupMember>, SemanticError> {
    categories
        .rows("ma_data_group")
        .iter()
        .enumerate()
        .map(|(row_index, row)| {
            Ok(DataGroupMember {
                ordinal_id: required_integer(row, "ma_data_group", "ordinal_id", row_index)?,
                group_id: required_integer(row, "ma_data_group", "group_id", row_index)?,
                data_id: required_integer(row, "ma_data_group", "data_id", row_index)?,
            })
        })
        .collect()
}

pub(super) fn targets(categories: &CategoryIndex<'_>) -> Result<Vec<TargetEntity>, SemanticError> {
    categories
        .rows("ma_target_entity")
        .iter()
        .enumerate()
        .map(|(row_index, row)| {
            Ok(TargetEntity {
                entity_id: required_text(row, "ma_target_entity", "entity_id", row_index)?,
                data_id: required_integer(row, "ma_target_entity", "data_id", row_index)?,
                origin: required_text(row, "ma_target_entity", "origin", row_index)?,
            })
        })
        .collect()
}

pub(super) fn target_instances(
    categories: &CategoryIndex<'_>,
) -> Result<Vec<TargetInstance>, SemanticError> {
    categories
        .rows("ma_target_entity_instance")
        .iter()
        .enumerate()
        .map(|(row_index, row)| {
            Ok(TargetInstance {
                asym_id: required_text(row, "ma_target_entity_instance", "asym_id", row_index)?,
                entity_id: required_text(row, "ma_target_entity_instance", "entity_id", row_index)?,
                details: optional_text(row, "details"),
            })
        })
        .collect()
}

pub(super) fn models(
    categories: &CategoryIndex<'_>,
) -> Result<Vec<PredictionModel>, SemanticError> {
    categories
        .rows("ma_model_list")
        .iter()
        .enumerate()
        .map(|(row_index, row)| {
            Ok(PredictionModel {
                id: required_integer(row, "ma_model_list", "ordinal_id", row_index)?,
                name: optional_text(row, "model_name"),
                assembly_id: optional_integer(row, "ma_model_list", "assembly_id", row_index)?,
                model_type: required_text(row, "ma_model_list", "model_type", row_index)?,
                type_details: optional_text(row, "model_type_other_details"),
                data_id: required_integer(row, "ma_model_list", "data_id", row_index)?,
            })
        })
        .collect()
}

pub(super) fn model_groups(
    categories: &CategoryIndex<'_>,
) -> Result<Vec<ModelGroup>, SemanticError> {
    categories
        .rows("ma_model_group")
        .iter()
        .enumerate()
        .map(|(row_index, row)| {
            Ok(ModelGroup {
                id: required_integer(row, "ma_model_group", "id", row_index)?,
                name: optional_text(row, "name"),
                details: optional_text(row, "details"),
            })
        })
        .collect()
}

pub(super) fn model_group_links(
    categories: &CategoryIndex<'_>,
) -> Result<Vec<ModelGroupLink>, SemanticError> {
    categories
        .rows("ma_model_group_link")
        .iter()
        .enumerate()
        .map(|(row_index, row)| {
            Ok(ModelGroupLink {
                model_id: required_integer(row, "ma_model_group_link", "model_id", row_index)?,
                group_id: required_integer(row, "ma_model_group_link", "group_id", row_index)?,
            })
        })
        .collect()
}

pub(super) fn representatives(
    categories: &CategoryIndex<'_>,
) -> Result<Vec<ModelRepresentative>, SemanticError> {
    categories
        .rows("ma_model_representative")
        .iter()
        .enumerate()
        .map(|(row_index, row)| {
            Ok(ModelRepresentative {
                id: required_integer(row, "ma_model_representative", "id", row_index)?,
                group_id: required_integer(
                    row,
                    "ma_model_representative",
                    "model_group_id",
                    row_index,
                )?,
                model_id: required_integer(row, "ma_model_representative", "model_id", row_index)?,
                selection_criteria: required_text(
                    row,
                    "ma_model_representative",
                    "selection_criteria",
                    row_index,
                )?,
            })
        })
        .collect()
}
