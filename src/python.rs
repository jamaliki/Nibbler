mod objects;

use std::sync::{Arc, Mutex, MutexGuard, PoisonError};

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::{PyAny, PyBytes, PyCapsule, PyDict, PyModule};

use crate::cif::{
    ArrowMissingPolicy, BinaryCifError, CifDocument, CifTable, InputError, InputLimits,
    LoadedBytes, MissingKind, ParseError, Predicate, ProjectionError, ProjectionPlan, SchemaName,
    SourceBuffer, decode_binary, decode_input, encode_binary, export_arrow_stream,
    is_binary_container, parse_source, project_binary_named, project_document, project_source,
    read_input_file, validate_document as validate_cif_document, write_canonical, write_preserving,
};
use crate::modelcif::{
    ModelCifModel, build_model_with_registry as build_modelcif_with_registry,
    canonical_document as canonical_modelcif_document,
    canonical_document_with_local_metric as canonical_modelcif_document_with_local_metric,
    validate_document as validate_modelcif,
};
use crate::pdbx::{
    ComponentRegistry, PdbxModel, SiteArray, assembly_document,
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
    fn with_model(
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
    pub(crate) fn source_copies<T>(&self, read: impl FnOnce(&Arc<PdbxModel>, &[u32], &[u64]) -> T) -> Option<T> {
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
        Ok(PyCifTable {
            table: Arc::new(table),
            missing: ArrowMissingPolicy::Collapse,
            include_provenance: false,
        })
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

/// Immutable chemical component definitions loaded by `nibbler.components.read`.
#[pyclass(name = "Registry", module = "nibbler.components", frozen)]
struct PyComponentRegistry {
    registry: Arc<ComponentRegistry>,
}

#[pyclass(name = "PreparedChemistry", module = "nibbler.mmcif", frozen)]
struct PyPreparedChemistry {
    graph: crate::pdbx::PreparedChemistry,
}

#[pyclass(name = "ConformerSelection", module = "nibbler.mmcif", frozen)]
struct PyConformerSelection {
    selected: crate::pdbx::ConformerSelection,
    source: Option<Arc<CifDocument>>,
    // The source's registry-free model, passed on when the selection keeps every row.
    model: Option<Arc<PdbxModel>>,
}

#[pymethods]
impl PyConformerSelection {
    #[getter]
    fn atom_rows(&self) -> Vec<usize> {
        self.selected.atom_rows.clone()
    }
    #[getter]
    fn search_nodes(&self) -> usize {
        self.selected.search_nodes
    }
    #[getter]
    fn model_number(&self) -> Option<i64> {
        self.selected.model_number
    }
    #[getter]
    fn schema_version(&self) -> &'static str {
        self.selected.schema_version
    }
    #[getter]
    fn document(&self, py: Python<'_>) -> PyResult<PyCifDocument> {
        let source = self
            .source
            .as_ref()
            .ok_or_else(|| PyValueError::new_err("Array selection has no source document"))?;
        let document = py
            .detach(|| crate::pdbx::selected_atom_document(source, &self.selected.atom_rows))
            .map_err(model_error_to_python)?;
        // Selecting every row returns the source blocks unchanged, so the source model holds.
        let model = self
            .model
            .as_ref()
            .filter(|_| std::ptr::eq(document.blocks(), source.blocks()))
            .map(Arc::clone);
        Ok(PyCifDocument::with_model(
            Arc::new(document),
            Some("pdbx".into()),
            model,
        ))
    }
}

#[pyfunction]
#[pyo3(signature = (model, model_number=None, seed=0, policy="occupancy"))]
fn select_mmcif_conformer(
    py: Python<'_>,
    model: &PyMmcifModel,
    model_number: Option<i64>,
    seed: u64,
    policy: &str,
) -> PyResult<PyConformerSelection> {
    let policy = conformer_policy(policy)?;
    let selected = py
        .detach(|| {
            crate::pdbx::select_conformer_with_policy(
                model.model.coordinates(),
                model_number,
                seed,
                policy,
            )
        })
        .map_err(model_error_to_python)?;
    Ok(PyConformerSelection {
        selected,
        source: Some(Arc::new(
            model.model.coordinates().source_document().clone(),
        )),
        model: None,
    })
}

#[pyfunction]
#[pyo3(signature = (document, model_number=None, seed=0, policy="occupancy"))]
fn select_document_conformer(
    py: Python<'_>,
    document: &PyCifDocument,
    model_number: Option<i64>,
    seed: u64,
    policy: &str,
) -> PyResult<PyConformerSelection> {
    let policy = conformer_policy(policy)?;
    let (selected, model) = py
        .detach(|| {
            // A built model decodes the same atom sites and links as the document path,
            // and later calls on this document reuse it. A document that does not form
            // a model still selects from its atom_site table alone.
            match document.pdbx_model() {
                Ok(model) => {
                    crate::pdbx::select_conformer_with_policy(&model, model_number, seed, policy)
                        .map(|selected| (selected, Some(model)))
                }
                Err(_) => crate::pdbx::select_conformer_document_with_policy(
                    &document.document,
                    model_number,
                    seed,
                    policy,
                )
                .map(|selected| (selected, None)),
            }
        })
        .map_err(model_error_to_python)?;
    Ok(PyConformerSelection {
        selected,
        source: Some(Arc::clone(&document.document)),
        model,
    })
}

fn conformer_policy(value: &str) -> PyResult<crate::pdbx::ConformerPolicy> {
    match value {
        "occupancy" => Ok(crate::pdbx::ConformerPolicy::Occupancy),
        "random_clash_aware" => Ok(crate::pdbx::ConformerPolicy::RandomClashAware),
        _ => Err(PyValueError::new_err("Unknown conformer policy")),
    }
}

#[pyfunction]
#[pyo3(signature = (coordinates, occupancy, residue_ids, chain_ids, sequence_ids, atom_names, alternate_ids, bonds, clash_mask, seed=0))]
#[allow(clippy::too_many_arguments)]
fn select_conformer_arrays(
    py: Python<'_>,
    coordinates: Vec<[f64; 3]>,
    occupancy: Vec<f64>,
    residue_ids: Vec<usize>,
    chain_ids: Vec<String>,
    sequence_ids: Vec<Option<i64>>,
    atom_names: Vec<String>,
    alternate_ids: Vec<Option<String>>,
    bonds: Vec<[usize; 2]>,
    clash_mask: Vec<bool>,
    seed: u64,
) -> PyResult<PyConformerSelection> {
    let count = coordinates.len();
    if [
        occupancy.len(),
        residue_ids.len(),
        chain_ids.len(),
        sequence_ids.len(),
        atom_names.len(),
        alternate_ids.len(),
        clash_mask.len(),
    ]
    .iter()
    .any(|&n| n != count)
    {
        return Err(PyValueError::new_err(
            "Conformer arrays must have equal lengths",
        ));
    }
    let atoms: Vec<_> = coordinates
        .into_iter()
        .zip(occupancy)
        .zip(residue_ids)
        .zip(chain_ids)
        .zip(sequence_ids)
        .zip(atom_names)
        .zip(alternate_ids)
        .zip(clash_mask)
        .map(
            |(
                ((((((coordinates, occupancy), residue), chain), sequence), atom_name), alternate),
                clash_eligible,
            )| {
                crate::pdbx::ConformerAtom {
                    coordinates,
                    occupancy,
                    residue,
                    chain,
                    sequence,
                    atom_name,
                    alternate,
                    clash_eligible,
                }
            },
        )
        .collect();
    let selected = py
        .detach(|| crate::pdbx::select_atom_conformers(&atoms, &bonds, seed))
        .map_err(model_error_to_python)?;
    Ok(PyConformerSelection {
        selected,
        source: None,
        model: None,
    })
}

#[pyfunction]
fn deduplicate_assembly_arrays(
    py: Python<'_>,
    positions: Vec<[f64; 3]>,
    source_atoms: Vec<usize>,
    units: Vec<usize>,
    copies: Vec<usize>,
    eligible: Vec<bool>,
    occupied: Vec<bool>,
) -> PyResult<(Vec<usize>, Vec<[usize; 2]>)> {
    let result = py
        .detach(|| {
            crate::pdbx::deduplicate_assembly_copies(
                &positions,
                &source_atoms,
                &units,
                &copies,
                &eligible,
                &occupied,
            )
        })
        .map_err(model_error_to_python)?;
    Ok((result.retained_rows, result.duplicate_pairs))
}

#[pymethods]
impl PyPreparedChemistry {
    #[getter]
    fn schema_version(&self) -> &'static str {
        self.graph.schema_version()
    }

    #[getter]
    fn atom_count(&self) -> usize {
        self.graph.atoms().len()
    }

    #[getter]
    fn bond_count(&self) -> usize {
        self.graph.bonds().len()
    }

    fn to_dict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, pyo3::types::PyAny>> {
        objects::to_python(py, &self.graph)
    }
}

#[pyfunction]
#[pyo3(signature = (model, registry=None, model_number=1, alternate_id=None, chem_data=None, use_observed_hydrogens=true))]
fn prepare_mmcif_chemistry(
    py: Python<'_>,
    model: &PyMmcifModel,
    registry: Option<&PyComponentRegistry>,
    model_number: i64,
    alternate_id: Option<String>,
    chem_data: Option<std::path::PathBuf>,
    use_observed_hydrogens: bool,
) -> PyResult<PyPreparedChemistry> {
    let graph = py
        .detach(|| {
            if registry.is_some() && chem_data.is_some() {
                return Err(crate::pdbx::SemanticError::new("PDBX_LINK_PRODUCT",
                    "Supply registry or chem_data, not both", Vec::new()));
            }
            let library = chem_data.as_deref().map(|root| {
                crate::pdbx::chemical_library_registry(model.model.coordinates(), root)
                    .map_err(|message| crate::pdbx::SemanticError::new("PDBX_LINK_PRODUCT", message, Vec::new()))
            }).transpose()?;
            crate::pdbx::prepare_product(
                model.model.coordinates(),
                registry.map(|r| r.registry.as_ref()).or(library.as_ref()),
                model_number,
                alternate_id.as_deref(),
                use_observed_hydrogens,
            )
        })
        .map_err(model_error_to_python)?;
    Ok(PyPreparedChemistry { graph })
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

/// Write out every copy of one biological assembly of a model's coordinates.
#[pyfunction]
#[pyo3(signature = (model, assembly_id=None))]
fn build_assembly_document(
    py: Python<'_>,
    model: &PyMmcifModel,
    assembly_id: Option<String>,
) -> PyResult<PyCifDocument> {
    let document = py
        .detach(|| assembly_document(model.model.coordinates(), assembly_id.as_deref()))
        .map_err(model_error_to_python)?;
    Ok(PyCifDocument::new(
        Arc::new(document),
        Some("pdbx".to_owned()),
    ))
}

type DiagnosticFields = (String, String, String, Vec<String>);
type AssemblyMappedDocument = (PyCifDocument, Vec<(usize, String, String)>);

#[pyfunction]
#[pyo3(signature = (model, assembly_id=None))]
fn build_assembly_with_mapping(
    py: Python<'_>,
    model: &PyMmcifModel,
    assembly_id: Option<String>,
) -> PyResult<AssemblyMappedDocument> {
    let source = model.model.coordinates();
    // The expansion also builds its own model from the source sites, so a later
    // molecular view of the result need not decode the written atom sites again.
    let (document, merges, expanded) = py
        .detach(|| crate::pdbx::assembly_with_model(source, assembly_id.as_deref()))
        .map_err(model_error_to_python)?;
    let mapping = crate::pdbx::assembly_mapping_names(source, merges);
    Ok((
        PyCifDocument::with_model(
            Arc::new(document),
            Some("pdbx".to_owned()),
            expanded.map(Arc::new),
        ),
        mapping,
    ))
}
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
    Ok(validation_fields(&report))
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
    let document = Arc::clone(&source.document);
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

fn validation_fields(report: &crate::cif::ValidationReport) -> ValidationFields {
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

/// A projected CIF category with lossless internal missing states.
#[pyclass(name = "CifTable", module = "nibbler", frozen)]
struct PyCifTable {
    table: Arc<CifTable>,
    missing: ArrowMissingPolicy,
    include_provenance: bool,
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
    schema: Option<String>,
) -> PyResult<Py<PyAny>> {
    match output {
        ReadOutput::Document(document) => Ok(Py::new(
            py,
            PyCifDocument::new(Arc::new(document), schema),
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
    module.add_function(wrap_pyfunction!(ligand_sasa_model, module)?)?;
    module.add_class::<PyCifDocument>()?;
    module.add_class::<PyCifTable>()?;
    module.add_class::<PyMmcifModel>()?;
    module.add_class::<PyComponentRegistry>()?;
    module.add_class::<PyPreparedChemistry>()?;
    module.add_class::<PyConformerSelection>()?;
    module.add_function(wrap_pyfunction!(select_mmcif_conformer, module)?)?;
    module.add_function(wrap_pyfunction!(select_document_conformer, module)?)?;
    module.add_function(wrap_pyfunction!(select_conformer_arrays, module)?)?;
    module.add_function(wrap_pyfunction!(deduplicate_assembly_arrays, module)?)?;
    module.add_function(wrap_pyfunction!(prepare_mmcif_chemistry, module)?)?;
    module.add_function(wrap_pyfunction!(read_file, module)?)?;
    module.add_function(wrap_pyfunction!(read_bytes, module)?)?;
    module.add_function(wrap_pyfunction!(validate_document, module)?)?;
    module.add_function(wrap_pyfunction!(build_mmcif_model, module)?)?;
    module.add_function(wrap_pyfunction!(_model_document, module)?)?;
    module.add_function(wrap_pyfunction!(build_assembly_document, module)?)?;
    module.add_function(wrap_pyfunction!(build_assembly_with_mapping, module)?)?;
    module.add_function(wrap_pyfunction!(build_component_registry, module)?)?;
    module.add_function(wrap_pyfunction!(validate_mmcif_document, module)?)?;
    module.add_function(wrap_pyfunction!(validate_mmcif_model, module)?)?;
    crate::python_scan::register(module)?;
    #[cfg(feature = "reduce3")]
    crate::python_reduce::register(module)?;
    #[cfg(feature = "reduce3")]
    crate::python_reduce_parents::register(module)?;
    Ok(())
}
