//! Global, local, and pairwise QA row decoding.

use crate::pdbx::category::{CategoryIndex, Row};

use super::fields::{
    optional_integer, optional_text, required_float, required_integer, required_text,
};
use super::qa::{QaMetric, QaValue, ResidueSite};
use super::source::ModelCifError;

pub(super) fn qa_metrics(categories: &CategoryIndex<'_>) -> Result<Vec<QaMetric>, ModelCifError> {
    categories
        .rows("ma_qa_metric")
        .iter()
        .enumerate()
        .map(|(row_index, row)| {
            Ok(QaMetric {
                id: required_integer(row, "ma_qa_metric", "id", row_index)?,
                name: required_text(row, "ma_qa_metric", "name", row_index)?,
                description: optional_text(row, "description"),
                metric_type: required_text(row, "ma_qa_metric", "type", row_index)?,
                mode: required_text(row, "ma_qa_metric", "mode", row_index)?,
                software_group_id: optional_integer(
                    row,
                    "ma_qa_metric",
                    "software_group_id",
                    row_index,
                )?,
                data_id: optional_integer(row, "ma_qa_metric", "data_id", row_index)?,
            })
        })
        .collect()
}

pub(super) fn qa_values(categories: &CategoryIndex<'_>) -> Result<Vec<QaValue>, ModelCifError> {
    let mut values = global_qa(categories)?;
    values.extend(local_qa(categories)?);
    values.extend(pairwise_qa(categories)?);
    Ok(values)
}

fn global_qa(categories: &CategoryIndex<'_>) -> Result<Vec<QaValue>, ModelCifError> {
    categories
        .rows("ma_qa_metric_global")
        .iter()
        .enumerate()
        .map(|(row_index, row)| {
            Ok(QaValue::Global {
                ordinal_id: required_integer(row, "ma_qa_metric_global", "ordinal_id", row_index)?,
                model_id: required_integer(row, "ma_qa_metric_global", "model_id", row_index)?,
                metric_id: required_integer(row, "ma_qa_metric_global", "metric_id", row_index)?,
                value: required_float(row, "ma_qa_metric_global", "metric_value", row_index)?,
            })
        })
        .collect()
}

fn local_qa(categories: &CategoryIndex<'_>) -> Result<Vec<QaValue>, ModelCifError> {
    let category = "ma_qa_metric_local";
    categories
        .rows(category)
        .iter()
        .enumerate()
        .map(|(row_index, row)| {
            Ok(QaValue::Local {
                ordinal_id: required_integer(row, category, "ordinal_id", row_index)?,
                model_id: required_integer(row, category, "model_id", row_index)?,
                site: residue_site(row, category, "", row_index)?,
                metric_id: required_integer(row, category, "metric_id", row_index)?,
                value: required_float(row, category, "metric_value", row_index)?,
            })
        })
        .collect()
}

fn pairwise_qa(categories: &CategoryIndex<'_>) -> Result<Vec<QaValue>, ModelCifError> {
    let category = "ma_qa_metric_local_pairwise";
    categories
        .rows(category)
        .iter()
        .enumerate()
        .map(|(row_index, row)| {
            Ok(QaValue::Pairwise {
                ordinal_id: required_integer(row, category, "ordinal_id", row_index)?,
                model_id: required_integer(row, category, "model_id", row_index)?,
                first: residue_site(row, category, "_1", row_index)?,
                second: residue_site(row, category, "_2", row_index)?,
                metric_id: required_integer(row, category, "metric_id", row_index)?,
                value: required_float(row, category, "metric_value", row_index)?,
            })
        })
        .collect()
}

fn residue_site(
    row: &Row<'_>,
    category: &str,
    suffix: &str,
    row_index: usize,
) -> Result<ResidueSite, ModelCifError> {
    Ok(ResidueSite {
        asym_id: required_text(row, category, &format!("label_asym_id{suffix}"), row_index)?,
        sequence_id: required_integer(row, category, &format!("label_seq_id{suffix}"), row_index)?,
        component_id: required_text(row, category, &format!("label_comp_id{suffix}"), row_index)?,
    })
}
