//! Direct BinaryCIF category projection into shared Arrow builders.

use std::collections::HashMap;

use crate::cif::document::ColumnValues;
use crate::cif::projection::{ColumnSpec, ProjectionPlan, predicate_matches};
use crate::cif::table::{CifCell, CifTable, ColumnBuilder, MissingKind, TableBuilder};

use super::codec::decode_data;
use super::document::{decode_file, decode_mask, normalize_category_name, shape_error};
use super::error::BinaryCifError;
use super::model::BinaryCategory;

struct ProjectedColumn {
    values: ColumnValues,
    mask: Vec<u8>,
}

/// Decode one BinaryCIF category directly into the shared columnar table model.
///
/// Only output and predicate columns are decoded. Category and column selectors use
/// the same [`ProjectionPlan`] contract as text CIF projection.
///
/// # Errors
///
/// Returns a structured BinaryCIF error for malformed containers, names without a text
/// CIF spelling, block headers, categories within a block, or columns within a category
/// repeated under ASCII case folding, malformed encoding chains, missing selected
/// columns, a category repeated across blocks with incompatible columns, or invalid
/// typed values.
pub fn project_binary(bytes: &[u8], plan: ProjectionPlan) -> Result<CifTable, BinaryCifError> {
    project_binary_named(bytes, plan, "<binary>")
}

pub(crate) fn project_binary_named(
    bytes: &[u8],
    plan: ProjectionPlan,
    source_name: &str,
) -> Result<CifTable, BinaryCifError> {
    let file = decode_file(bytes)?;
    let mut table: Option<TableBuilder> = None;
    let mut expected_columns: Option<Vec<ColumnSpec>> = None;
    for block in file.data_blocks {
        for category in block.categories {
            let category_name = normalize_category_name(&category.name)?;
            if !category_name.eq_ignore_ascii_case(&plan.category) {
                continue;
            }
            let output_specs = output_specs(&plan, &category)?;
            if let Some(expected) = &expected_columns {
                if !same_specs(expected, &output_specs) {
                    return Err(shape_error(
                        "matching BinaryCIF categories have incompatible column layouts",
                    ));
                }
            } else {
                table = Some(TableBuilder {
                    category: plan.category.clone(),
                    columns: output_specs
                        .iter()
                        .map(|spec| ColumnBuilder::new(spec.name.clone(), spec.column_type))
                        .collect(),
                    provenance: Vec::new(),
                });
                expected_columns = Some(output_specs.clone());
            }
            let Some(output_table) = table.as_mut() else {
                return Err(shape_error(
                    "BinaryCIF projection table was not initialized",
                ));
            };
            project_category(
                category,
                &plan,
                &output_specs,
                output_table,
                source_name,
                &block.header,
            )?;
        }
    }
    Ok(table
        .unwrap_or_else(|| TableBuilder {
            category: plan.category,
            columns: plan.columns.map_or_else(Vec::new, |columns| {
                columns
                    .into_iter()
                    .map(|spec| ColumnBuilder::new(spec.name, spec.column_type))
                    .collect()
            }),
            provenance: Vec::new(),
        })
        .finish())
}

fn output_specs(
    plan: &ProjectionPlan,
    category: &BinaryCategory<'_>,
) -> Result<Vec<ColumnSpec>, BinaryCifError> {
    if let Some(columns) = &plan.columns {
        return Ok(columns.clone());
    }
    category
        .columns
        .iter()
        .map(|column| {
            let column_type = plan.column_type(&column.name).ok_or_else(|| {
                shape_error(format!(
                    "item _{}.{} is not defined by the selected schema",
                    plan.category, column.name
                ))
            })?;
            Ok(ColumnSpec::new(&column.name, column_type))
        })
        .collect()
}

fn project_category(
    category: BinaryCategory<'_>,
    plan: &ProjectionPlan,
    output_specs: &[ColumnSpec],
    table: &mut TableBuilder,
    source_name: &str,
    block_code: &str,
) -> Result<(), BinaryCifError> {
    let prepared = plan.prepare_columns(output_specs.to_vec());
    let needed = &prepared.columns;
    let mut available = category
        .columns
        .into_iter()
        .map(|column| (column.name.to_ascii_lowercase(), column))
        .collect::<HashMap<_, _>>();
    let mut decoded = Vec::with_capacity(needed.len());
    for spec in needed {
        let column = available.remove(&spec.key).ok_or_else(|| {
            shape_error(format!(
                "category {:?} does not contain required column {:?}",
                plan.category, spec.name
            ))
        })?;
        let values = decode_data(column.data)?;
        if values.len() != category.row_count {
            return Err(shape_error(format!(
                "_{}.{} decodes {} rows, expected {}",
                plan.category,
                spec.name,
                values.len(),
                category.row_count
            )));
        }
        let mask = column
            .mask
            .map(decode_mask)
            .transpose()?
            .unwrap_or_else(|| vec![0; category.row_count]);
        if mask.len() != category.row_count {
            return Err(shape_error(format!(
                "_{}.{} mask decodes {} rows, expected {}",
                plan.category,
                spec.name,
                mask.len(),
                category.row_count
            )));
        }
        decoded.push(ProjectedColumn { values, mask });
    }
    if prepared.predicate_columns.is_empty() {
        for row_index in 0..category.row_count {
            for (builder, column) in table
                .columns
                .iter_mut()
                .zip(decoded.iter().take(output_specs.len()))
            {
                append_at(builder, column, row_index)?;
            }
        }
        table.append_provenance_run(source_name, Some(block_code), None, category.row_count);
        return Ok(());
    }
    let mut row = Vec::with_capacity(decoded.len());
    for row_index in 0..category.row_count {
        row.clear();
        for column in &decoded {
            row.push(cell_at(column, row_index)?);
        }
        if !plan
            .predicates
            .iter()
            .zip(&prepared.predicate_columns)
            .all(|(predicate, &column)| predicate_matches(&predicate.predicate, &row[column]))
        {
            continue;
        }
        for (column, value) in table
            .columns
            .iter_mut()
            .zip(row.drain(..output_specs.len()))
        {
            column.append(value).map_err(shape_error)?;
        }
        table.append_provenance(source_name, Some(block_code), None);
    }
    Ok(())
}

fn append_at(
    builder: &mut ColumnBuilder,
    column: &ProjectedColumn,
    row_index: usize,
) -> Result<(), BinaryCifError> {
    match column.mask[row_index] {
        1 => {
            builder.append_missing(MissingKind::NotApplicable);
            return Ok(());
        }
        2 => {
            builder.append_missing(MissingKind::Unknown);
            return Ok(());
        }
        0 => {}
        value => return Err(shape_error(format!("invalid BinaryCIF mask value {value}"))),
    }
    match &column.values {
        ColumnValues::Integers(values) => builder.append_integer_value(values[row_index]),
        ColumnValues::Floats(values) => {
            let number = values[row_index];
            if !number.is_finite() {
                return Err(shape_error(
                    "BinaryCIF present floating-point values must be finite",
                ));
            }
            builder.append_float_value(number).map_err(shape_error)?;
        }
        ColumnValues::Strings(values) => builder
            .append_text_value(values.value(row_index))
            .map_err(shape_error)?,
    }
    Ok(())
}

fn cell_at(column: &ProjectedColumn, row_index: usize) -> Result<CifCell, BinaryCifError> {
    match column.mask[row_index] {
        1 => return Ok(CifCell::NotApplicable),
        2 => return Ok(CifCell::Unknown),
        0 => {}
        value => return Err(shape_error(format!("invalid BinaryCIF mask value {value}"))),
    }
    Ok(match &column.values {
        ColumnValues::Integers(values) => CifCell::Text(values[row_index].to_string()),
        ColumnValues::Floats(values) => {
            let value = values[row_index];
            if !value.is_finite() {
                return Err(shape_error(
                    "BinaryCIF present floating-point values must be finite",
                ));
            }
            CifCell::Text(value.to_string())
        }
        ColumnValues::Strings(values) => CifCell::Text(values.value(row_index).to_owned()),
    })
}

fn same_specs(left: &[ColumnSpec], right: &[ColumnSpec]) -> bool {
    left.len() == right.len()
        && left
            .iter()
            .zip(right)
            .all(|(left, right)| left.key == right.key && left.column_type == right.column_type)
}
