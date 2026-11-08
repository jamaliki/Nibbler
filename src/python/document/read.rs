//! Reading text CIF or BinaryCIF from a file or bytes into a document or projected table,
//! and the Python form of read failures.

use std::sync::Arc;

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::{PyAny, PyModule};

use super::PyCifDocument;
use crate::cif::{
    BinaryCifError, CifDocument, CifTable, InputError, InputLimits, LoadedBytes, ParseError,
    ProjectionError, ProjectionPlan, SourceBuffer, decode_binary, decode_input,
    is_binary_container, parse_source, project_binary_named, project_source, read_input_file,
};
use crate::python::table::{PredicateSpec, PyCifTable, build_optional_plan};

pub(crate) enum InputSource {
    File(String),
    Bytes { source_name: String, bytes: Vec<u8> },
}

pub(crate) enum ReadOutput {
    Document(CifDocument),
    Table(CifTable),
}

pub(crate) enum ReadFailure {
    Input(InputError),
    Parse(ParseError),
    Projection(ProjectionError),
    Binary {
        source_name: String,
        error: BinaryCifError,
    },
}

enum LoadedSource {
    Text(SourceBuffer),
    Binary { source_name: String, bytes: Vec<u8> },
}

#[pyfunction]
#[pyo3(signature = (file, category=None, columns=None, predicates=Vec::new(), schema=None))]
fn read_file(
    py: Python<'_>,
    file: String,
    category: Option<String>,
    columns: Option<Vec<String>>,
    predicates: Vec<PredicateSpec>,
    schema: Option<String>,
) -> PyResult<Py<PyAny>> {
    read_input(
        py,
        InputSource::File(file),
        category,
        columns,
        predicates,
        schema,
    )
}

#[pyfunction]
#[pyo3(signature = (source_name, bytes, category=None, columns=None, predicates=Vec::new(), schema=None))]
fn read_bytes(
    py: Python<'_>,
    source_name: String,
    bytes: Vec<u8>,
    category: Option<String>,
    columns: Option<Vec<String>>,
    predicates: Vec<PredicateSpec>,
    schema: Option<String>,
) -> PyResult<Py<PyAny>> {
    read_input(
        py,
        InputSource::Bytes { source_name, bytes },
        category,
        columns,
        predicates,
        schema,
    )
}

fn read_input(
    py: Python<'_>,
    input: InputSource,
    category: Option<String>,
    columns: Option<Vec<String>>,
    predicates: Vec<PredicateSpec>,
    schema: Option<String>,
) -> PyResult<Py<PyAny>> {
    let plan = build_optional_plan(category, columns, predicates, schema.clone())
        .map_err(error_to_python)?;
    let output = py
        .detach(move || execute_read(input, plan))
        .map_err(error_to_python)?;
    output_to_python(py, output, false, schema)
}

pub(crate) fn execute_read(
    input: InputSource,
    plan: Option<ProjectionPlan>,
) -> Result<ReadOutput, ReadFailure> {
    let source = load_source(input)?;
    match (source, plan) {
        (LoadedSource::Text(source), None) => parse_source(source)
            .map(ReadOutput::Document)
            .map_err(ReadFailure::Parse),
        (LoadedSource::Text(source), Some(plan)) => project_source(source, plan)
            .map(ReadOutput::Table)
            .map_err(ReadFailure::Projection),
        (LoadedSource::Binary { source_name, bytes }, None) => decode_binary(&bytes)
            .map(ReadOutput::Document)
            .map_err(|error| ReadFailure::Binary { source_name, error }),
        (LoadedSource::Binary { source_name, bytes }, Some(plan)) => {
            project_binary_named(&bytes, plan, &source_name)
                .map(ReadOutput::Table)
                .map_err(|error| ReadFailure::Binary { source_name, error })
        }
    }
}

fn load_source(input: InputSource) -> Result<LoadedSource, ReadFailure> {
    let loaded = match input {
        InputSource::File(file) => read_input_file(file, InputLimits::default()),
        InputSource::Bytes { source_name, bytes } => {
            decode_input(source_name, bytes, InputLimits::default())
        }
    }
    .map_err(ReadFailure::Input)?;
    classify_source(loaded)
}

fn classify_source(loaded: LoadedBytes) -> Result<LoadedSource, ReadFailure> {
    if is_binary_container(&loaded.bytes) {
        return Ok(LoadedSource::Binary {
            source_name: loaded.source_name,
            bytes: loaded.bytes,
        });
    }
    SourceBuffer::from_owned_bytes(loaded.source_name, loaded.bytes)
        .map(LoadedSource::Text)
        .map_err(ReadFailure::Parse)
}

pub(crate) fn output_to_python(
    py: Python<'_>,
    output: ReadOutput,
    include_provenance: bool,
    schema: Option<String>,
) -> PyResult<Py<PyAny>> {
    match output {
        ReadOutput::Document(document) => {
            Ok(Py::new(py, PyCifDocument::new(Arc::new(document), schema))?.into_any())
        }
        ReadOutput::Table(table) => {
            Ok(Py::new(py, PyCifTable::new(Arc::new(table), include_provenance))?.into_any())
        }
    }
}

pub(crate) type ErrorFields = (
    String,
    String,
    Option<String>,
    Option<usize>,
    Option<usize>,
    Option<usize>,
    Option<usize>,
);

pub(crate) fn error_fields(error: &ReadFailure) -> ErrorFields {
    match error {
        ReadFailure::Input(error) => (
            error.code_str().to_owned(),
            error.message().to_owned(),
            Some(error.source_name().to_owned()),
            None,
            None,
            None,
            None,
        ),
        ReadFailure::Parse(error) => parse_error_fields(error),
        ReadFailure::Projection(error) => error.parse_error().map_or_else(
            || {
                (
                    error.code().as_str().to_owned(),
                    error.message().to_owned(),
                    None,
                    None,
                    None,
                    None,
                    None,
                )
            },
            parse_error_fields,
        ),
        ReadFailure::Binary { source_name, error } => (
            error.code().as_str().to_owned(),
            error.message().to_owned(),
            Some(source_name.clone()),
            None,
            None,
            None,
            None,
        ),
    }
}

pub(crate) fn error_to_python(error: ReadFailure) -> PyErr {
    PyValueError::new_err(error_fields(&error))
}

fn parse_error_fields(error: &ParseError) -> ErrorFields {
    let span = error.span();
    (
        error.code_str().to_owned(),
        error.message().to_owned(),
        Some(error.source_name().to_owned()),
        Some(span.line()),
        Some(span.column()),
        Some(span.byte_start()),
        Some(span.byte_end()),
    )
}

pub(super) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(read_file, module)?)?;
    module.add_function(wrap_pyfunction!(read_bytes, module)?)?;
    Ok(())
}
