//! The `CifTable` class, projection plans built from Python arguments, and Arrow export.

use std::sync::Arc;

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::{PyAny, PyCapsule, PyModule};

use super::document::read::ReadFailure;
use crate::cif::{
    ArrowMissingPolicy, CifTable, MissingKind, Predicate, ProjectionError, ProjectionPlan,
    SchemaName, export_arrow_stream,
};

pub(crate) type PredicateSpec = (String, String, Vec<String>);

/// A projected CIF category with lossless internal missing states.
#[pyclass(name = "CifTable", module = "nibbler", frozen)]
pub(crate) struct PyCifTable {
    table: Arc<CifTable>,
    missing: ArrowMissingPolicy,
    include_provenance: bool,
}

impl PyCifTable {
    /// Wrap a projected table with the default collapsed missing-value representation.
    pub(crate) const fn new(table: Arc<CifTable>, include_provenance: bool) -> Self {
        Self {
            table,
            missing: ArrowMissingPolicy::Collapse,
            include_provenance,
        }
    }
}

#[pymethods]
impl PyCifTable {
    #[getter]
    /// Return the category name without a leading underscore.
    fn category(&self) -> &str {
        self.table.category()
    }

    #[getter]
    /// Return projected item names in output order.
    fn columns(&self) -> Vec<String> {
        self.table
            .columns()
            .iter()
            .map(|column| column.name().to_owned())
            .collect()
    }

    fn __len__(&self) -> usize {
        self.table.row_count()
    }

    /// Return a view with the selected missing-value Arrow representation.
    fn with_missing(&self, missing: &str) -> PyResult<Self> {
        Ok(Self {
            table: Arc::clone(&self.table),
            missing: parse_missing_policy(missing)?,
            include_provenance: self.include_provenance,
        })
    }

    #[pyo3(signature = (requested_schema=None))]
    /// Export this table through the zero-copy Arrow C Stream protocol.
    fn __arrow_c_stream__<'py>(
        &self,
        py: Python<'py>,
        requested_schema: Option<&Bound<'py, PyAny>>,
    ) -> PyResult<Bound<'py, PyCapsule>> {
        if requested_schema.is_some_and(|schema| !schema.is_none()) {
            return Err(PyValueError::new_err(
                "requested Arrow schema negotiation is not supported",
            ));
        }
        let stream = export_arrow_stream(&self.table, self.missing, self.include_provenance)
            .map_err(|error| PyValueError::new_err(error.to_string()))?;
        PyCapsule::new_with_value(py, stream, c"arrow_array_stream")
    }

    fn __repr__(&self) -> String {
        format!(
            "CifTable(category={:?}, rows={}, columns={})",
            self.table.category(),
            self.table.row_count(),
            self.table.columns().len()
        )
    }
}

pub(crate) fn build_optional_plan(
    category: Option<String>,
    columns: Option<Vec<String>>,
    predicates: Vec<PredicateSpec>,
    schema: Option<String>,
) -> Result<Option<ProjectionPlan>, ReadFailure> {
    let Some(category) = category else {
        return Ok(None);
    };
    let mut plan = ProjectionPlan::new(category).map_err(ReadFailure::Projection)?;
    if let Some(columns) = columns {
        plan = plan
            .with_columns(columns)
            .map_err(ReadFailure::Projection)?;
    }
    for (column, operator, operands) in predicates {
        let predicate = predicate_from_spec(column, &operator, operands)?;
        plan = plan
            .with_predicate(predicate)
            .map_err(ReadFailure::Projection)?;
    }
    if let Some(schema) = schema {
        let schema = schema.parse::<SchemaName>().map_err(|error| {
            ReadFailure::Projection(ProjectionError::schema(error.code(), error.message()))
        })?;
        plan = plan.with_schema(schema).map_err(ReadFailure::Projection)?;
    }
    Ok(Some(plan))
}

fn predicate_from_spec(
    column: String,
    operator: &str,
    operands: Vec<String>,
) -> Result<Predicate, ReadFailure> {
    let invalid = || {
        ReadFailure::Projection(ProjectionError::invalid_predicate(format!(
            "invalid predicate operator or operands: {operator:?}"
        )))
    };
    match (operator, operands.as_slice()) {
        ("eq", [value]) => Ok(Predicate::Equal {
            column,
            value: value.clone(),
        }),
        ("ne", [value]) => Ok(Predicate::NotEqual {
            column,
            value: value.clone(),
        }),
        ("in", _) => Ok(Predicate::In {
            column,
            values: operands,
        }),
        ("missing", []) => Ok(Predicate::Missing { column, kind: None }),
        ("missing", [value]) if value == "?" => Ok(Predicate::Missing {
            column,
            kind: Some(MissingKind::Unknown),
        }),
        ("missing", [value]) if value == "." => Ok(Predicate::Missing {
            column,
            kind: Some(MissingKind::NotApplicable),
        }),
        _ => Err(invalid()),
    }
}

fn parse_missing_policy(policy: &str) -> PyResult<ArrowMissingPolicy> {
    match policy {
        "collapse" => Ok(ArrowMissingPolicy::Collapse),
        "columns" => Ok(ArrowMissingPolicy::Columns),
        "extension" => Ok(ArrowMissingPolicy::Extension),
        _ => Err(PyValueError::new_err(format!(
            "missing must be 'collapse', 'columns', or 'extension', not {policy:?}"
        ))),
    }
}

pub(super) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyCifTable>()
}
