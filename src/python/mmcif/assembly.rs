//! Biological assembly expansion and duplicate-copy removal.

use std::sync::Arc;

use pyo3::prelude::*;
use pyo3::types::PyModule;

use super::{PyMmcifModel, model_error_to_python};
use crate::pdbx::assembly_document;
use crate::python::document::PyCifDocument;

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

pub(super) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(build_assembly_document, module)?)?;
    module.add_function(wrap_pyfunction!(build_assembly_with_mapping, module)?)?;
    module.add_function(wrap_pyfunction!(deduplicate_assembly_arrays, module)?)?;
    Ok(())
}
