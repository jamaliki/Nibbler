//! The `MmcifModel` and component `Registry` classes: semantic PDBx and ModelCIF models,
//! their canonical documents, and their validation.

mod assembly;
mod chemistry;
mod conformer;
mod objects;

use std::sync::Arc;

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyModule;

use super::document::{PyCifDocument, ValidationFields, validation_fields};
use crate::cif::SchemaName;
use crate::modelcif::{
    ModelCifModel, build_model_with_registry as build_modelcif_with_registry,
    canonical_document as canonical_modelcif_document,
    canonical_document_with_local_metric as canonical_modelcif_document_with_local_metric,
    validate_document as validate_modelcif,
};
use crate::pdbx::{
    ComponentRegistry, PdbxModel, build_component_registry as load_component_registry,
    build_model_with_registry, canonical_document, validate_document as validate_pdbx,
};

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

/// An immutable, source-backed macromolecular coordinate model.
#[pyclass(name = "MmcifModel", module = "nibbler", frozen)]
struct PyMmcifModel {
    model: NativeMmcifModel,
}

#[pymethods]
impl PyMmcifModel {
    #[getter]
    /// Return the semantic profile used to construct this model.
    fn profile(&self) -> &'static str {
        match self.model {
            NativeMmcifModel::Pdbx(_) => "pdbx",
            NativeMmcifModel::ModelCif(_) => "modelcif",
        }
    }

    #[getter]
    /// Return the PDBx entry identifier.
    fn entry_id(&self) -> &str {
        self.model.coordinates().entry_id()
    }

    #[getter]
    /// Return the number of semantic entities.
    fn entity_count(&self) -> usize {
        self.model.coordinates().entity_count()
    }

    #[getter]
    /// Return the number of entity instances.
    fn asym_unit_count(&self) -> usize {
        self.model.coordinates().asym_unit_count()
    }

    #[getter]
    /// Return the number of resolved or explicitly unresolved components.
    fn component_count(&self) -> usize {
        self.model.coordinates().component_count()
    }

    #[getter]
    /// Return the number of coordinate atom sites.
    fn atom_site_count(&self) -> usize {
        self.model.coordinates().atom_site_count()
    }

    #[getter]
    /// Return the number of explicit inter-site connections.
    fn connection_count(&self) -> usize {
        self.model.coordinates().connection_count()
    }

    #[getter]
    /// Return entity kinds in source order.
    fn entity_kinds(&self) -> Vec<&str> {
        self.model.coordinates().entity_kinds().collect()
    }

    #[getter]
    /// Return component identifiers in stable order.
    fn component_ids(&self) -> Vec<&str> {
        self.model.coordinates().component_ids().collect()
    }

    #[getter]
    /// Return the number of deposited prediction models.
    fn prediction_model_count(&self) -> usize {
        match &self.model {
            NativeMmcifModel::Pdbx(_) => 0,
            NativeMmcifModel::ModelCif(model) => model.prediction_model_count(),
        }
    }

    #[getter]
    /// Return the number of declared prediction targets.
    fn target_entity_count(&self) -> usize {
        match &self.model {
            NativeMmcifModel::Pdbx(_) => 0,
            NativeMmcifModel::ModelCif(model) => model.target_entity_count(),
        }
    }

    #[getter]
    /// Return the number of structural templates.
    fn template_count(&self) -> usize {
        match &self.model {
            NativeMmcifModel::Pdbx(_) => 0,
            NativeMmcifModel::ModelCif(model) => model.template_count(),
        }
    }

    #[getter]
    /// Return the number of quality-metric definitions.
    fn qa_metric_count(&self) -> usize {
        match &self.model {
            NativeMmcifModel::Pdbx(_) => 0,
            NativeMmcifModel::ModelCif(model) => model.qa_metric_count(),
        }
    }

    #[getter]
    /// Return the number of global, local, and pairwise QA values.
    fn qa_value_count(&self) -> usize {
        match &self.model {
            NativeMmcifModel::Pdbx(_) => 0,
            NativeMmcifModel::ModelCif(model) => model.qa_value_count(),
        }
    }

    #[getter]
    /// Return prediction software names in source order.
    fn software_names(&self) -> Vec<&str> {
        match &self.model {
            NativeMmcifModel::Pdbx(_) => Vec::new(),
            NativeMmcifModel::ModelCif(model) => model.software_names().collect(),
        }
    }

    #[getter]
    /// Return quality-metric names in source order.
    fn qa_metric_names(&self) -> Vec<&str> {
        match &self.model {
            NativeMmcifModel::Pdbx(_) => Vec::new(),
            NativeMmcifModel::ModelCif(model) => model.qa_metric_names().collect(),
        }
    }

    #[getter]
    /// Return quality-metric modes in source order.
    fn qa_metric_modes(&self) -> Vec<&str> {
        match &self.model {
            NativeMmcifModel::Pdbx(_) => Vec::new(),
            NativeMmcifModel::ModelCif(model) => model.qa_metric_modes().collect(),
        }
    }

    /// Return the generic document from which this model was built.
    fn source_document(&self) -> PyCifDocument {
        PyCifDocument::new(
            Arc::new(self.model.coordinates().source_document().clone()),
            Some(self.profile().to_owned()),
        )
    }

    fn __repr__(&self) -> String {
        format!(
            "MmcifModel(profile={:?}, entry_id={:?}, entities={}, atoms={})",
            self.profile(),
            self.entry_id(),
            self.entity_count(),
            self.atom_site_count()
        )
    }
}

/// Immutable chemical component definitions loaded by `nibbler.components.read`.
#[pyclass(name = "Registry", module = "nibbler.components", frozen)]
struct PyComponentRegistry {
    registry: Arc<ComponentRegistry>,
}

#[pymethods]
impl PyComponentRegistry {
    fn __len__(&self) -> usize {
        self.registry.len()
    }

    fn __repr__(&self) -> String {
        format!("Registry(components={})", self.registry.len())
    }
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
    let source = document;
    let document = Arc::clone(source.document());
    let registry = registry.map(|value| Arc::clone(&value.registry));
    let model = match profile {
        SchemaName::Pdbx if registry.is_none() => {
            let model = py
                .detach(|| source.pdbx_model())
                .map_err(model_error_to_python)?;
            NativeMmcifModel::Pdbx(model)
        }
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
    let document = Arc::clone(document.document());
    let registry = py
        .detach(move || load_component_registry(&document))
        .map_err(model_error_to_python)?;
    Ok(PyComponentRegistry {
        registry: Arc::new(registry),
    })
}

#[pyfunction]
#[pyo3(signature = (model, mirror_local_qa_metric=None))]
fn _model_document(
    model: &PyMmcifModel,
    mirror_local_qa_metric: Option<i64>,
) -> PyResult<PyCifDocument> {
    let document = match (&model.model, mirror_local_qa_metric) {
        (NativeMmcifModel::Pdbx(model), None) => canonical_document(model),
        (NativeMmcifModel::Pdbx(_), Some(_)) => {
            return Err(PyValueError::new_err(
                "mirror_local_qa_metric requires a ModelCIF model",
            ));
        }
        (NativeMmcifModel::ModelCif(model), Some(metric_id)) => {
            canonical_modelcif_document_with_local_metric(model, metric_id).map_err(|error| {
                PyValueError::new_err((error.code().to_owned(), error.message().to_owned()))
            })?
        }
        (NativeMmcifModel::ModelCif(model), None) => canonical_modelcif_document(model),
    };
    Ok(PyCifDocument::new(
        Arc::new(document),
        Some(model.profile().to_owned()),
    ))
}

/// Remove a document's H, D and T atom sites and their anisotropic rows.
///
/// A document that loses no row is returned itself, with any model it has built.
#[pyfunction]
fn hydrogen_free_document(
    py: Python<'_>,
    document: Py<PyCifDocument>,
) -> PyResult<Py<PyCifDocument>> {
    let (source, schema) = {
        let borrowed = document.borrow(py);
        (
            Arc::clone(borrowed.document()),
            borrowed.schema_name().map(str::to_owned),
        )
    };
    let result = py
        .detach(|| crate::pdbx::hydrogen_free_document(&source))
        .map_err(model_error_to_python)?;
    if std::ptr::eq(result.blocks(), source.blocks()) {
        return Ok(document);
    }
    Py::new(py, PyCifDocument::new(Arc::new(result), schema))
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
    let document = Arc::clone(document.document());
    match profile {
        SchemaName::Pdbx => {
            let report = py
                .detach(move || validate_pdbx(&document))
                .map_err(|error| {
                    PyValueError::new_err((error.code(), error.message().to_owned()))
                })?;
            Ok(validation_fields(&report))
        }
        SchemaName::ModelCif => {
            let report = py
                .detach(move || validate_modelcif(&document))
                .map_err(|error| {
                    PyValueError::new_err((error.code(), error.message().to_owned()))
                })?;
            Ok(validation_fields(&report))
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
            Ok(validation_fields(&report))
        }
        NativeMmcifModel::ModelCif(model) => {
            let report = py
                .detach(move || validate_modelcif(&canonical_modelcif_document(&model)))
                .map_err(|error| {
                    PyValueError::new_err((error.code(), error.message().to_owned()))
                })?;
            Ok(validation_fields(&report))
        }
    }
}

#[pyfunction]
#[pyo3(signature = (model, atom_rows=None, model_number=1, probe_radius=1.4, point_count=256, radii=None))]
fn ligand_sasa_model(
    py: Python<'_>,
    model: &PyMmcifModel,
    atom_rows: Option<Vec<usize>>,
    model_number: i64,
    probe_radius: f64,
    point_count: usize,
    radii: Option<std::collections::HashMap<String, f64>>,
) -> PyResult<(Vec<usize>, Vec<f64>)> {
    let radii = radii.unwrap_or_default();
    py.detach(|| {
        crate::pdbx::model_ligand_sasa(
            model.model.coordinates(),
            atom_rows.as_deref(),
            model_number,
            probe_radius,
            point_count,
            &radii,
        )
    })
    .map_err(PyValueError::new_err)
}

pub(super) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(ligand_sasa_model, module)?)?;
    module.add_class::<PyMmcifModel>()?;
    module.add_class::<PyComponentRegistry>()?;
    module.add_function(wrap_pyfunction!(build_mmcif_model, module)?)?;
    module.add_function(wrap_pyfunction!(_model_document, module)?)?;
    module.add_function(wrap_pyfunction!(hydrogen_free_document, module)?)?;
    module.add_function(wrap_pyfunction!(build_component_registry, module)?)?;
    module.add_function(wrap_pyfunction!(validate_mmcif_document, module)?)?;
    module.add_function(wrap_pyfunction!(validate_mmcif_model, module)?)?;
    chemistry::register(module)?;
    conformer::register(module)?;
    assembly::register(module)
}
