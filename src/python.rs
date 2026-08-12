use std::sync::Arc;

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::{PyAny, PyCapsule, PyModule};

use crate::cif::{
    ArrowMissingPolicy, BinaryCifError, CifDocument, CifTable, InputError, InputLimits,
    LoadedBytes, MissingKind, ParseError, Predicate, ProjectionError, ProjectionPlan, SchemaName,
    SourceBuffer, decode_binary, decode_input, encode_binary, export_arrow_stream,
    is_binary_container, parse_source, project_binary_named, project_source, read_input_file,
    validate_document as validate_cif_document, write_canonical, write_preserving,
};
use crate::modelcif::{
    MirrorPolicy, ModelCifModel, build_model_with_registry as build_modelcif_with_registry,
    canonical_document as canonical_modelcif_document,
    canonical_document_with_mirror as canonical_modelcif_document_with_mirror,
    validate_document as validate_modelcif,
};
use crate::pdbx::{
    ComponentRegistry, PdbxModel, ProfileValidationReport,
    build_component_registry as load_component_registry, build_model_with_registry,
    canonical_document, validate_document as validate_pdbx,
};

pub(crate) type PredicateSpec = (String, String, Vec<String>);

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

#[pyclass(name = "_CifDocument", frozen)]
struct PyCifDocument {
    document: Arc<CifDocument>,
}

#[pymethods]
impl PyCifDocument {
    #[getter]
    fn block_count(&self) -> usize {
        self.document.blocks().len()
    }

    fn to_canonical(&self, py: Python<'_>) -> PyResult<String> {
        let document = Arc::clone(&self.document);
        py.detach(move || write_canonical(&document))
            .map_err(|error| PyValueError::new_err((error.code_str(), error.message().to_owned())))
    }

    fn to_preserving(&self, py: Python<'_>) -> PyResult<String> {
        let document = Arc::clone(&self.document);
        py.detach(move || write_preserving(&document))
            .map_err(|error| PyValueError::new_err((error.code_str(), error.message().to_owned())))
    }

    fn to_binary(&self, py: Python<'_>) -> PyResult<Vec<u8>> {
        let document = Arc::clone(&self.document);
        py.detach(move || encode_binary(&document))
            .map_err(|error| {
                PyValueError::new_err((
                    error.code().as_str().to_owned(),
                    error.message().to_owned(),
                ))
            })
    }

    fn __repr__(&self) -> String {
        format!("_CifDocument(block_count={})", self.document.blocks().len())
    }
}

#[derive(Clone)]
enum NativeMmcifModel {
    Pdbx(Arc<PdbxModel>),
    ModelCif(Arc<ModelCifModel>),
}

impl NativeMmcifModel {
    fn coordinates(&self) -> &PdbxModel {
        match self {
            Self::Pdbx(model) => model,
            Self::ModelCif(model) => model.coordinates(),
        }
    }
}

#[pyclass(name = "_MmcifModel", frozen)]
struct PyMmcifModel {
    model: NativeMmcifModel,
}

#[pyclass(name = "_ComponentRegistry", frozen)]
struct PyComponentRegistry {
    registry: Arc<ComponentRegistry>,
}

#[pymethods]
impl PyComponentRegistry {
    fn __len__(&self) -> usize {
        self.registry.len()
    }

    fn __repr__(&self) -> String {
        format!("_ComponentRegistry(components={})", self.registry.len())
    }
}

#[pymethods]
impl PyMmcifModel {
    #[getter]
    fn profile(&self) -> &'static str {
        match self.model {
            NativeMmcifModel::Pdbx(_) => "pdbx",
            NativeMmcifModel::ModelCif(_) => "modelcif",
        }
    }

    #[getter]
    fn entry_id(&self) -> &str {
        self.model.coordinates().entry_id()
    }

    #[getter]
    fn entity_count(&self) -> usize {
        self.model.coordinates().entities().len()
    }

    #[getter]
    fn asym_unit_count(&self) -> usize {
        self.model.coordinates().asym_units().len()
    }

    #[getter]
    fn component_count(&self) -> usize {
        self.model.coordinates().components().len()
    }

    #[getter]
    fn atom_site_count(&self) -> usize {
        self.model.coordinates().atom_sites().len()
    }

    #[getter]
    fn connection_count(&self) -> usize {
        self.model.coordinates().connections().len()
    }

    #[getter]
    fn entity_kinds(&self) -> Vec<&str> {
        self.model
            .coordinates()
            .entities()
            .iter()
            .map(|entity| entity.kind().as_str())
            .collect()
    }

    #[getter]
    fn component_ids(&self) -> Vec<&str> {
        self.model
            .coordinates()
            .components()
            .iter()
            .map(|component| component.id())
            .collect()
    }

    #[getter]
    fn prediction_model_count(&self) -> usize {
        match &self.model {
            NativeMmcifModel::Pdbx(_) => 0,
            NativeMmcifModel::ModelCif(model) => model.models().len(),
        }
    }

    #[getter]
    fn target_entity_count(&self) -> usize {
        match &self.model {
            NativeMmcifModel::Pdbx(_) => 0,
            NativeMmcifModel::ModelCif(model) => model.targets().len(),
        }
    }

    #[getter]
    fn template_count(&self) -> usize {
        match &self.model {
            NativeMmcifModel::Pdbx(_) => 0,
            NativeMmcifModel::ModelCif(model) => model.templates().len(),
        }
    }

    #[getter]
    fn qa_metric_count(&self) -> usize {
        match &self.model {
            NativeMmcifModel::Pdbx(_) => 0,
            NativeMmcifModel::ModelCif(model) => model.qa_metrics().len(),
        }
    }

    #[getter]
    fn qa_value_count(&self) -> usize {
        match &self.model {
            NativeMmcifModel::Pdbx(_) => 0,
            NativeMmcifModel::ModelCif(model) => model.qa_values().len(),
        }
    }

    #[getter]
    fn software_names(&self) -> Vec<&str> {
        match &self.model {
            NativeMmcifModel::Pdbx(_) => Vec::new(),
            NativeMmcifModel::ModelCif(model) => model
                .software()
                .iter()
                .map(|software| software.name())
                .collect(),
        }
    }

    #[getter]
    fn qa_metric_names(&self) -> Vec<&str> {
        match &self.model {
            NativeMmcifModel::Pdbx(_) => Vec::new(),
            NativeMmcifModel::ModelCif(model) => model
                .qa_metrics()
                .iter()
                .map(|metric| metric.name())
                .collect(),
        }
    }

    #[getter]
    fn qa_metric_modes(&self) -> Vec<&str> {
        match &self.model {
            NativeMmcifModel::Pdbx(_) => Vec::new(),
            NativeMmcifModel::ModelCif(model) => model
                .qa_metrics()
                .iter()
                .map(|metric| metric.mode())
                .collect(),
        }
    }

    #[pyo3(signature = (mirror_local_qa_metric=None))]
    fn to_document(&self, mirror_local_qa_metric: Option<i64>) -> PyResult<PyCifDocument> {
        let document = match (&self.model, mirror_local_qa_metric) {
            (NativeMmcifModel::Pdbx(model), None) => canonical_document(model),
            (NativeMmcifModel::Pdbx(_), Some(_)) => {
                return Err(PyValueError::new_err(
                    "mirror_local_qa_metric requires a ModelCIF model",
                ));
            }
            (NativeMmcifModel::ModelCif(model), Some(metric_id)) => {
                canonical_modelcif_document_with_mirror(model, MirrorPolicy::LocalMetric(metric_id))
                    .map_err(|error| {
                        PyValueError::new_err((error.code().to_owned(), error.message().to_owned()))
                    })?
            }
            (NativeMmcifModel::ModelCif(model), None) => canonical_modelcif_document(model),
        };
        Ok(PyCifDocument {
            document: Arc::new(document),
        })
    }

    fn source_document(&self) -> PyCifDocument {
        PyCifDocument {
            document: Arc::new(self.model.coordinates().source_document().clone()),
        }
    }

    fn __repr__(&self) -> String {
        format!(
            "_MmcifModel(profile={:?}, entry_id={:?}, entities={}, atoms={})",
            self.profile(),
            self.entry_id(),
            self.entity_count(),
            self.atom_site_count()
        )
    }
}

type DiagnosticFields = (String, String, String, Vec<String>);
type ValidationFields = (String, String, Vec<String>, Vec<DiagnosticFields>);

#[pyfunction]
fn validate_document(
    py: Python<'_>,
    document: &PyCifDocument,
    schema: &str,
) -> PyResult<ValidationFields> {
    let schema = schema
        .parse::<SchemaName>()
        .map_err(|error| PyValueError::new_err((error.code(), error.message().to_owned())))?;
    let document = Arc::clone(&document.document);
    let report = py
        .detach(move || validate_cif_document(&document, schema))
        .map_err(|error| PyValueError::new_err((error.code(), error.message().to_owned())))?;
    let diagnostics = report
        .diagnostics()
        .iter()
        .map(|diagnostic| {
            (
                diagnostic.code().to_owned(),
                diagnostic.severity().as_str().to_owned(),
                diagnostic.message().to_owned(),
                diagnostic.context().to_vec(),
            )
        })
        .collect();
    Ok((
        report.schema_name().to_owned(),
        report.dictionary_version().to_owned(),
        report.coverage().to_vec(),
        diagnostics,
    ))
}

#[pyfunction]
#[pyo3(signature = (document, profile, registry=None))]
fn build_mmcif_model(
    py: Python<'_>,
    document: &PyCifDocument,
    profile: &str,
    registry: Option<&PyComponentRegistry>,
) -> PyResult<PyMmcifModel> {
    let profile = profile
        .parse::<SchemaName>()
        .map_err(|error| PyValueError::new_err((error.code(), error.message().to_owned())))?;
    let document = Arc::clone(&document.document);
    let registry = registry.map(|value| Arc::clone(&value.registry));
    let model = match profile {
        SchemaName::Pdbx => {
            let model = py
                .detach(move || build_model_with_registry(&document, registry.as_deref()))
                .map_err(model_error_to_python)?;
            NativeMmcifModel::Pdbx(Arc::new(model))
        }
        SchemaName::ModelCif => {
            let model = py
                .detach(move || build_modelcif_with_registry(&document, registry.as_deref()))
                .map_err(model_error_to_python)?;
            NativeMmcifModel::ModelCif(Arc::new(model))
        }
    };
    Ok(PyMmcifModel { model })
}

#[pyfunction]
fn build_component_registry(
    py: Python<'_>,
    document: &PyCifDocument,
) -> PyResult<PyComponentRegistry> {
    let document = Arc::clone(&document.document);
    let registry = py
        .detach(move || load_component_registry(&document))
        .map_err(model_error_to_python)?;
    Ok(PyComponentRegistry {
        registry: Arc::new(registry),
    })
}

fn model_error_to_python(error: crate::pdbx::SemanticError) -> PyErr {
    PyValueError::new_err((
        error.code().to_owned(),
        error.message().to_owned(),
        error.context().to_vec(),
    ))
}

#[pyfunction]
fn validate_mmcif_document(
    py: Python<'_>,
    document: &PyCifDocument,
    profile: &str,
) -> PyResult<ValidationFields> {
    let profile = profile
        .parse::<SchemaName>()
        .map_err(|error| PyValueError::new_err((error.code(), error.message().to_owned())))?;
    let document = Arc::clone(&document.document);
    match profile {
        SchemaName::Pdbx => {
            let report = py
                .detach(move || validate_pdbx(&document))
                .map_err(|error| {
                    PyValueError::new_err((error.code(), error.message().to_owned()))
                })?;
            Ok(profile_fields("pdbx", &report))
        }
        SchemaName::ModelCif => {
            let report = py
                .detach(move || validate_modelcif(&document))
                .map_err(|error| {
                    PyValueError::new_err((error.code(), error.message().to_owned()))
                })?;
            Ok(profile_fields("modelcif", &report))
        }
    }
}

#[pyfunction]
fn validate_mmcif_model(py: Python<'_>, model: &PyMmcifModel) -> PyResult<ValidationFields> {
    let model = model.model.clone();
    match model {
        NativeMmcifModel::Pdbx(model) => {
            let report = py
                .detach(move || validate_pdbx(&canonical_document(&model)))
                .map_err(|error| {
                    PyValueError::new_err((error.code(), error.message().to_owned()))
                })?;
            Ok(profile_fields("pdbx", &report))
        }
        NativeMmcifModel::ModelCif(model) => {
            let report = py
                .detach(move || validate_modelcif(&canonical_modelcif_document(&model)))
                .map_err(|error| {
                    PyValueError::new_err((error.code(), error.message().to_owned()))
                })?;
            Ok(profile_fields("modelcif", &report))
        }
    }
}

fn profile_fields(profile: &str, report: &ProfileValidationReport) -> ValidationFields {
    let diagnostics = report
        .diagnostics()
        .iter()
        .map(|diagnostic| {
            (
                diagnostic.code().to_owned(),
                diagnostic.severity().as_str().to_owned(),
                diagnostic.message().to_owned(),
                diagnostic.context().to_vec(),
            )
        })
        .collect();
    (
        profile.to_owned(),
        report.dictionary_version().to_owned(),
        report
            .coverage()
            .iter()
            .map(|value| (*value).to_owned())
            .collect(),
        diagnostics,
    )
}

#[pyclass(name = "_CifTable", frozen)]
struct PyCifTable {
    table: Arc<CifTable>,
    missing: ArrowMissingPolicy,
    include_provenance: bool,
}

#[pymethods]
impl PyCifTable {
    #[getter]
    fn category(&self) -> &str {
        self.table.category()
    }

    #[getter]
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

    fn _with_missing(&self, missing: &str) -> PyResult<Self> {
        Ok(Self {
            table: Arc::clone(&self.table),
            missing: parse_missing_policy(missing)?,
            include_provenance: self.include_provenance,
        })
    }

    #[pyo3(signature = (requested_schema=None))]
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
            "_CifTable(category={:?}, rows={}, columns={})",
            self.table.category(),
            self.table.row_count(),
            self.table.columns().len()
        )
    }
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
    let plan =
        build_optional_plan(category, columns, predicates, schema).map_err(error_to_python)?;
    let output = py
        .detach(move || execute_read(input, plan))
        .map_err(error_to_python)?;
    output_to_python(py, output, false)
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

pub(crate) fn output_to_python(
    py: Python<'_>,
    output: ReadOutput,
    include_provenance: bool,
) -> PyResult<Py<PyAny>> {
    match output {
        ReadOutput::Document(document) => Ok(Py::new(
            py,
            PyCifDocument {
                document: Arc::new(document),
            },
        )?
        .into_any()),
        ReadOutput::Table(table) => Ok(Py::new(
            py,
            PyCifTable {
                table: Arc::new(table),
                missing: ArrowMissingPolicy::Collapse,
                include_provenance,
            },
        )?
        .into_any()),
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

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyCifDocument>()?;
    module.add_class::<PyCifTable>()?;
    module.add_class::<PyMmcifModel>()?;
    module.add_class::<PyComponentRegistry>()?;
    module.add_function(wrap_pyfunction!(read_file, module)?)?;
    module.add_function(wrap_pyfunction!(read_bytes, module)?)?;
    module.add_function(wrap_pyfunction!(validate_document, module)?)?;
    module.add_function(wrap_pyfunction!(build_mmcif_model, module)?)?;
    module.add_function(wrap_pyfunction!(build_component_registry, module)?)?;
    module.add_function(wrap_pyfunction!(validate_mmcif_document, module)?)?;
    module.add_function(wrap_pyfunction!(validate_mmcif_model, module)?)?;
    crate::python_scan::register(module)?;
    Ok(())
}
