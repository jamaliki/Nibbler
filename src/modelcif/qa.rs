//! Quality metric records and row decoding.

use crate::cif::CategoryView;
use crate::pdbx::SemanticError;
use crate::pdbx::fields::category_rows;

use super::fields::{
    optional_integer, optional_text, required_float, required_integer, required_text,
};

#[derive(Clone, Debug, Eq, PartialEq)]
pub(super) struct QaMetric {
    pub id: i64,
    pub name: String,
    pub description: Option<String>,
    pub metric_type: String,
    pub mode: String,
    pub software_group_id: Option<i64>,
    pub data_id: Option<i64>,
}

#[derive(Clone, Debug, Eq, PartialEq, Ord, PartialOrd)]
pub(super) struct ResidueSite {
    pub asym_id: String,
    pub sequence_id: i64,
    pub component_id: String,
}

#[derive(Clone, Debug, PartialEq)]
pub(super) enum QaValue {
    Global {
        ordinal_id: i64,
        model_id: i64,
        metric_id: i64,
        value: f64,
    },
    Local {
        ordinal_id: i64,
        model_id: i64,
        site: ResidueSite,
        metric_id: i64,
        value: f64,
    },
    Pairwise {
        ordinal_id: i64,
        model_id: i64,
        first: ResidueSite,
        second: ResidueSite,
        metric_id: i64,
        value: f64,
    },
}

impl QaValue {
    pub(super) const fn model_id(&self) -> i64 {
        match self {
            Self::Global { model_id, .. }
            | Self::Local { model_id, .. }
            | Self::Pairwise { model_id, .. } => *model_id,
        }
    }

    pub(super) const fn metric_id(&self) -> i64 {
        match self {
            Self::Global { metric_id, .. }
            | Self::Local { metric_id, .. }
            | Self::Pairwise { metric_id, .. } => *metric_id,
        }
    }
}

pub(super) fn qa_metrics(categories: CategoryView<'_>) -> Result<Vec<QaMetric>, SemanticError> {
    category_rows(
        categories,
        "ma_qa_metric",
        [
            "id",
            "name",
            "description",
            "type",
            "mode",
            "software_group_id",
            "data_id",
        ],
    )
    .map(|row| {
        Ok(QaMetric {
            id: required_integer(&row, "ma_qa_metric", "id", 0)?,
            name: required_text(&row, "ma_qa_metric", "name", 1)?,
            description: optional_text(&row, 2),
            metric_type: required_text(&row, "ma_qa_metric", "type", 3)?,
            mode: required_text(&row, "ma_qa_metric", "mode", 4)?,
            software_group_id: optional_integer(&row, "ma_qa_metric", "software_group_id", 5)?,
            data_id: optional_integer(&row, "ma_qa_metric", "data_id", 6)?,
        })
    })
    .collect()
}

pub(super) fn qa_values(categories: CategoryView<'_>) -> Result<Vec<QaValue>, SemanticError> {
    let mut values = global_qa(categories)?;
    values.extend(local_qa(categories)?);
    values.extend(pairwise_qa(categories)?);
    Ok(values)
}

fn global_qa(categories: CategoryView<'_>) -> Result<Vec<QaValue>, SemanticError> {
    category_rows(
        categories,
        "ma_qa_metric_global",
        ["ordinal_id", "model_id", "metric_id", "metric_value"],
    )
    .map(|row| {
        Ok(QaValue::Global {
            ordinal_id: required_integer(&row, "ma_qa_metric_global", "ordinal_id", 0)?,
            model_id: required_integer(&row, "ma_qa_metric_global", "model_id", 1)?,
            metric_id: required_integer(&row, "ma_qa_metric_global", "metric_id", 2)?,
            value: required_float(&row, "ma_qa_metric_global", "metric_value", 3)?,
        })
    })
    .collect()
}

fn local_qa(categories: CategoryView<'_>) -> Result<Vec<QaValue>, SemanticError> {
    let category = "ma_qa_metric_local";
    category_rows(
        categories,
        category,
        [
            "ordinal_id",
            "model_id",
            "label_asym_id",
            "label_seq_id",
            "label_comp_id",
            "metric_id",
            "metric_value",
        ],
    )
    .map(|row| {
        Ok(QaValue::Local {
            ordinal_id: required_integer(&row, category, "ordinal_id", 0)?,
            model_id: required_integer(&row, category, "model_id", 1)?,
            site: ResidueSite {
                asym_id: required_text(&row, category, "label_asym_id", 2)?,
                sequence_id: required_integer(&row, category, "label_seq_id", 3)?,
                component_id: required_text(&row, category, "label_comp_id", 4)?,
            },
            metric_id: required_integer(&row, category, "metric_id", 5)?,
            value: required_float(&row, category, "metric_value", 6)?,
        })
    })
    .collect()
}

fn pairwise_qa(categories: CategoryView<'_>) -> Result<Vec<QaValue>, SemanticError> {
    let category = "ma_qa_metric_local_pairwise";
    category_rows(
        categories,
        category,
        [
            "ordinal_id",
            "model_id",
            "label_asym_id_1",
            "label_seq_id_1",
            "label_comp_id_1",
            "label_asym_id_2",
            "label_seq_id_2",
            "label_comp_id_2",
            "metric_id",
            "metric_value",
        ],
    )
    .map(|row| {
        Ok(QaValue::Pairwise {
            ordinal_id: required_integer(&row, category, "ordinal_id", 0)?,
            model_id: required_integer(&row, category, "model_id", 1)?,
            first: ResidueSite {
                asym_id: required_text(&row, category, "label_asym_id_1", 2)?,
                sequence_id: required_integer(&row, category, "label_seq_id_1", 3)?,
                component_id: required_text(&row, category, "label_comp_id_1", 4)?,
            },
            second: ResidueSite {
                asym_id: required_text(&row, category, "label_asym_id_2", 5)?,
                sequence_id: required_integer(&row, category, "label_seq_id_2", 6)?,
                component_id: required_text(&row, category, "label_comp_id_2", 7)?,
            },
            metric_id: required_integer(&row, category, "metric_id", 8)?,
            value: required_float(&row, category, "metric_value", 9)?,
        })
    })
    .collect()
}
