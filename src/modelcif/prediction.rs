//! Prediction graph records and row decoding.

use crate::cif::CategoryView;
use crate::pdbx::SemanticError;
use crate::pdbx::fields::category_rows;

use super::fields::{optional_integer, optional_text, required_integer, required_text};
use super::provenance::DataGroupMember;

#[derive(Clone, Debug, Eq, PartialEq)]
pub(super) struct Data {
    pub id: i64,
    pub name: String,
    pub content_type: String,
}
#[derive(Clone, Debug, Eq, PartialEq)]
pub(super) struct TargetEntity {
    pub entity_id: String,
    pub data_id: i64,
    pub origin: String,
}
#[derive(Clone, Debug, Eq, PartialEq)]
pub(super) struct TargetInstance {
    pub asym_id: String,
    pub entity_id: String,
    pub details: Option<String>,
}
#[derive(Clone, Debug, Eq, PartialEq)]
pub(super) struct PredictionModel {
    pub id: i64,
    pub name: Option<String>,
    pub assembly_id: Option<i64>,
    pub model_type: String,
    pub type_details: Option<String>,
    pub data_id: i64,
}
#[derive(Clone, Debug, Eq, PartialEq)]
pub(super) struct ModelGroup {
    pub id: i64,
    pub name: Option<String>,
    pub details: Option<String>,
}
#[derive(Clone, Debug, Eq, PartialEq)]
pub(super) struct ModelGroupLink {
    pub model_id: i64,
    pub group_id: i64,
}
#[derive(Clone, Debug, Eq, PartialEq)]
pub(super) struct ModelRepresentative {
    pub id: i64,
    pub group_id: i64,
    pub model_id: i64,
    pub selection_criteria: String,
}

pub(super) fn audit_conform(
    categories: CategoryView<'_>,
) -> Result<Vec<(String, String)>, SemanticError> {
    category_rows(categories, "audit_conform", ["dict_name", "dict_version"])
        .map(|row| {
            Ok((
                required_text(&row, "audit_conform", "dict_name", 0)?,
                required_text(&row, "audit_conform", "dict_version", 1)?,
            ))
        })
        .collect()
}

pub(super) fn data(categories: CategoryView<'_>) -> Result<Vec<Data>, SemanticError> {
    category_rows(categories, "ma_data", ["id", "name", "content_type"])
        .map(|row| {
            Ok(Data {
                id: required_integer(&row, "ma_data", "id", 0)?,
                name: required_text(&row, "ma_data", "name", 1)?,
                content_type: required_text(&row, "ma_data", "content_type", 2)?,
            })
        })
        .collect()
}

pub(super) fn data_groups(
    categories: CategoryView<'_>,
) -> Result<Vec<DataGroupMember>, SemanticError> {
    category_rows(
        categories,
        "ma_data_group",
        ["ordinal_id", "group_id", "data_id"],
    )
    .map(|row| {
        Ok(DataGroupMember {
            ordinal_id: required_integer(&row, "ma_data_group", "ordinal_id", 0)?,
            group_id: required_integer(&row, "ma_data_group", "group_id", 1)?,
            data_id: required_integer(&row, "ma_data_group", "data_id", 2)?,
        })
    })
    .collect()
}

pub(super) fn targets(categories: CategoryView<'_>) -> Result<Vec<TargetEntity>, SemanticError> {
    category_rows(
        categories,
        "ma_target_entity",
        ["entity_id", "data_id", "origin"],
    )
    .map(|row| {
        Ok(TargetEntity {
            entity_id: required_text(&row, "ma_target_entity", "entity_id", 0)?,
            data_id: required_integer(&row, "ma_target_entity", "data_id", 1)?,
            origin: required_text(&row, "ma_target_entity", "origin", 2)?,
        })
    })
    .collect()
}

pub(super) fn target_instances(
    categories: CategoryView<'_>,
) -> Result<Vec<TargetInstance>, SemanticError> {
    category_rows(
        categories,
        "ma_target_entity_instance",
        ["asym_id", "entity_id", "details"],
    )
    .map(|row| {
        Ok(TargetInstance {
            asym_id: required_text(&row, "ma_target_entity_instance", "asym_id", 0)?,
            entity_id: required_text(&row, "ma_target_entity_instance", "entity_id", 1)?,
            details: optional_text(&row, 2),
        })
    })
    .collect()
}

pub(super) fn models(categories: CategoryView<'_>) -> Result<Vec<PredictionModel>, SemanticError> {
    category_rows(
        categories,
        "ma_model_list",
        [
            "ordinal_id",
            "model_name",
            "assembly_id",
            "model_type",
            "model_type_other_details",
            "data_id",
        ],
    )
    .map(|row| {
        Ok(PredictionModel {
            id: required_integer(&row, "ma_model_list", "ordinal_id", 0)?,
            name: optional_text(&row, 1),
            assembly_id: optional_integer(&row, "ma_model_list", "assembly_id", 2)?,
            model_type: required_text(&row, "ma_model_list", "model_type", 3)?,
            type_details: optional_text(&row, 4),
            data_id: required_integer(&row, "ma_model_list", "data_id", 5)?,
        })
    })
    .collect()
}

pub(super) fn model_groups(categories: CategoryView<'_>) -> Result<Vec<ModelGroup>, SemanticError> {
    category_rows(categories, "ma_model_group", ["id", "name", "details"])
        .map(|row| {
            Ok(ModelGroup {
                id: required_integer(&row, "ma_model_group", "id", 0)?,
                name: optional_text(&row, 1),
                details: optional_text(&row, 2),
            })
        })
        .collect()
}

pub(super) fn model_group_links(
    categories: CategoryView<'_>,
) -> Result<Vec<ModelGroupLink>, SemanticError> {
    category_rows(categories, "ma_model_group_link", ["model_id", "group_id"])
        .map(|row| {
            Ok(ModelGroupLink {
                model_id: required_integer(&row, "ma_model_group_link", "model_id", 0)?,
                group_id: required_integer(&row, "ma_model_group_link", "group_id", 1)?,
            })
        })
        .collect()
}

pub(super) fn representatives(
    categories: CategoryView<'_>,
) -> Result<Vec<ModelRepresentative>, SemanticError> {
    category_rows(
        categories,
        "ma_model_representative",
        ["id", "model_group_id", "model_id", "selection_criteria"],
    )
    .map(|row| {
        Ok(ModelRepresentative {
            id: required_integer(&row, "ma_model_representative", "id", 0)?,
            group_id: required_integer(&row, "ma_model_representative", "model_group_id", 1)?,
            model_id: required_integer(&row, "ma_model_representative", "model_id", 2)?,
            selection_criteria: required_text(
                &row,
                "ma_model_representative",
                "selection_criteria",
                3,
            )?,
        })
    })
    .collect()
}
