//! The `CifDocument` class and generic dictionary validation of documents.

pub(crate) mod read;

use std::sync::{Arc, Mutex, MutexGuard, PoisonError};

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::{PyBytes, PyDict, PyModule};

use self::read::{ReadFailure, error_to_python};
use super::table::{PyCifTable, build_optional_plan};
use crate::cif::{
    CifDocument, SchemaName, encode_binary, project_document,
    validate_document as validate_cif_document, write_canonical, write_preserving,
};
use crate::pdbx::{PdbxModel, SiteArray, build_model_with_registry};

/// An immutable, order-preserving generic CIF document.
#[pyclass(name = "CifDocument", module = "nibbler", frozen)]
pub(crate) struct PyCifDocument {
    document: Arc<CifDocument>,
    schema: Option<String>,
    // The PDBx model of `document` built without a registry. Calls that need it share
    // one build, so the file's atom sites are decoded once.
    model: Mutex<ModelSlot>,
}

/// The registry-free PDBx model of a document, or what builds it.
enum ModelSlot {
    Empty,
    Built(Arc<PdbxModel>),
    /// The error of the build that failed. The document does not change, so a new build
    /// would fail the same way.
    Failed(crate::pdbx::SemanticError),
    /// The model of the document Reduce3 read, and the `_atom_site` rows of this document
    /// that copy its rows (see [`crate::pdbx::build_model_with_copies`]).
    #[cfg_attr(not(feature = "reduce3"), allow(dead_code))]
    Copied {
        source: Arc<PdbxModel>,
        rows: Vec<u32>,
        made: Vec<u64>,
    },
}

impl PyCifDocument {
    pub(crate) const fn new(document: Arc<CifDocument>, schema: Option<String>) -> Self {
        Self {
            document,
            schema,
            model: Mutex::new(ModelSlot::Empty),
        }
    }

    /// Wrap a document together with its registry-free PDBx model, if the caller has it.
    pub(crate) fn with_model(
        document: Arc<CifDocument>,
        schema: Option<String>,
        model: Option<Arc<PdbxModel>>,
    ) -> Self {
        Self {
            document,
            schema,
            model: Mutex::new(model.map_or(ModelSlot::Empty, ModelSlot::Built)),
        }
    }

    /// Wrap a document that Reduce3 wrote from the document of `source`. Row `i` of its
    /// `_atom_site` loop copies source row `rows[i]` in the columns that `made[i]` leaves
    /// clear, so its model copies those values from the sites of `source`.
    #[cfg_attr(not(feature = "reduce3"), allow(dead_code))]
    pub(crate) const fn with_copies(
        document: Arc<CifDocument>,
        schema: Option<String>,
        source: Arc<PdbxModel>,
        rows: Vec<u32>,
        made: Vec<u64>,
    ) -> Self {
        Self {
            document,
            schema,
            model: Mutex::new(ModelSlot::Copied { source, rows, made }),
        }
    }

    fn model_slot(&self) -> MutexGuard<'_, ModelSlot> {
        // A slot holds a whole model or none, so a panic elsewhere cannot leave it partial.
        self.model.lock().unwrap_or_else(PoisonError::into_inner)
    }

    /// Return the registry-free PDBx model of this document, building it on first use.
    pub(crate) fn pdbx_model(&self) -> Result<Arc<PdbxModel>, crate::pdbx::SemanticError> {
        let mut slot = self.model_slot();
        if let Some(model) = self.built(&mut slot) {
            return Ok(model);
        }
        if let ModelSlot::Failed(error) = &*slot {
            return Err(error.clone());
        }
        match build_model_with_registry(&self.document, None) {
            Ok(model) => {
                let model = Arc::new(model);
                *slot = ModelSlot::Built(Arc::clone(&model));
                Ok(model)
            }
            Err(error) => {
                *slot = ModelSlot::Failed(error.clone());
                Err(error)
            }
        }
    }

    /// Return the registry-free PDBx model if an earlier call built it, or build it from
    /// the model of the document Reduce3 read. The document keeps the model.
    #[cfg_attr(not(feature = "reduce3"), allow(dead_code))]
    pub(crate) fn built_pdbx_model(&self) -> Option<Arc<PdbxModel>> {
        let mut slot = self.model_slot();
        self.built(&mut slot)
    }

    /// Call `read` with the model, rows and masks of the `Copied` slot of a document that
    /// Reduce3 wrote, while its own model is unbuilt.
    #[cfg(feature = "reduce3")]
    pub(crate) fn source_copies<T>(
        &self,
        read: impl FnOnce(&Arc<PdbxModel>, &[u32], &[u64]) -> T,
    ) -> Option<T> {
        match &*self.model_slot() {
            ModelSlot::Copied { source, rows, made } => Some(read(source, rows, made)),
            _ => None,
        }
    }

    /// Remove the registry-free PDBx model, so the document stops holding it.
    #[cfg_attr(not(feature = "reduce3"), allow(dead_code))]
    pub(crate) fn release_pdbx_model(&self) {
        *self.model_slot() = ModelSlot::Empty;
    }

    fn built(&self, slot: &mut ModelSlot) -> Option<Arc<PdbxModel>> {
        let model = match std::mem::replace(slot, ModelSlot::Empty) {
            ModelSlot::Empty => return None,
            ModelSlot::Failed(error) => {
                *slot = ModelSlot::Failed(error);
                return None;
            }
            ModelSlot::Built(model) => model,
            // On an error the slot stays empty, and a full build reports the error.
            ModelSlot::Copied { source, rows, made } => Arc::new(
                crate::pdbx::build_model_with_copies(&self.document, &source, &rows, &made).ok()?,
            ),
        };
        *slot = ModelSlot::Built(Arc::clone(&model));
        Some(model)
    }

    pub(crate) const fn document(&self) -> &Arc<CifDocument> {
        &self.document
    }

    pub(crate) fn schema_name(&self) -> Option<&str> {
        self.schema.as_deref()
    }
}

#[pymethods]
impl PyCifDocument {
    #[getter]
    /// Return the number of data and global blocks.
    fn block_count(&self) -> usize {
        self.document.blocks().len()
    }

    #[getter]
    /// Return the dictionary selector requested at read time, if any.
    fn schema(&self) -> Option<&str> {
        self.schema.as_deref()
    }

    /// Return `_atom_site` items of this document's PDBx model laid out for numpy, or
    /// None when the model cannot be built or does not hold every item.
    ///
    /// The document builds the model once and keeps it, so later model reads such as
    /// conformer selection reuse it. A text item maps to `(width, codes)`, where `codes`
    /// holds `width` little-endian UCS-4 code points for each row, padded with zeros. Present
    /// text is the text a projection with the PDBx schema reads, and a missing value reads
    /// `?` or `.`. `label_seq_id` and `pdbx_PDB_model_num` map to `(values, kinds)`:
    /// little-endian `int64` values, 0 where missing, and the missing kind of each row, 0
    /// when present, 1 for `?` and 2 for `.`.
    fn atom_site_arrays<'py>(
        &self,
        py: Python<'py>,
        items: Vec<String>,
    ) -> PyResult<Option<Bound<'py, PyDict>>> {
        let arrays = py.detach(|| {
            let model = self.pdbx_model().ok()?;
            crate::pdbx::atom_site_arrays(&self.document, &model, &items)
        });
        let Some(arrays) = arrays else {
            return Ok(None);
        };
        let result = PyDict::new(py);
        for (item, array) in items.iter().zip(arrays) {
            match array {
                SiteArray::Text { width, codes } => {
                    result.set_item(item, (width, PyBytes::new(py, &codes)))?;
                }
                SiteArray::Integer { values, kinds } => {
                    result.set_item(item, (PyBytes::new(py, &values), PyBytes::new(py, &kinds)))?;
                }
            }
        }
        Ok(Some(result))
    }

    /// Project a category from this document without serialization or parsing.
    #[pyo3(signature = (category, *, columns=None, schema=None))]
    fn project(
        &self,
        py: Python<'_>,
        category: String,
        columns: Option<Vec<String>>,
        schema: Option<String>,
    ) -> PyResult<PyCifTable> {
        let plan = build_optional_plan(Some(category), columns, Vec::new(), schema)
            .map_err(error_to_python)?
            .ok_or_else(|| PyValueError::new_err("Projection requires a category"))?;
        let document = Arc::clone(&self.document);
        let table = py
            .detach(move || project_document(&document, plan))
            .map_err(|error| error_to_python(ReadFailure::Projection(error)))?;
        Ok(PyCifTable::new(Arc::new(table), false))
    }

    /// Serialize the document using Nibbler's canonical CIF policy.
    fn to_canonical(&self, py: Python<'_>) -> PyResult<String> {
        let document = Arc::clone(&self.document);
        py.detach(move || write_canonical(&document))
            .map_err(|error| PyValueError::new_err((error.code_str(), error.message().to_owned())))
    }

    /// Serialize with source ordering and valid original value lexemes.
    fn to_preserving(&self, py: Python<'_>) -> PyResult<String> {
        let document = Arc::clone(&self.document);
        py.detach(move || write_preserving(&document))
            .map_err(|error| PyValueError::new_err((error.code_str(), error.message().to_owned())))
    }

    /// Serialize the document as BinaryCIF 0.3 MessagePack.
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
        format!("CifDocument(block_count={})", self.document.blocks().len())
    }
}

pub(crate) type DiagnosticFields = (String, String, String, Vec<String>);

pub(crate) type ValidationFields = (String, String, Vec<String>, Vec<DiagnosticFields>);

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
    Ok(validation_fields(&report))
}

pub(crate) fn validation_fields(report: &crate::cif::ValidationReport) -> ValidationFields {
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
        report.schema_name().to_owned(),
        report.dictionary_version().to_owned(),
        report.coverage().to_vec(),
        diagnostics,
    )
}

pub(super) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyCifDocument>()?;
    module.add_function(wrap_pyfunction!(validate_document, module)?)?;
    read::register(module)
}
