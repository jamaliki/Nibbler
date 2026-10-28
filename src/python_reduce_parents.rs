//! Hydrogen parents from the same fixed-mode dictionaries as Reduce3.

use crate::observations;
use crate::observations::parents::{Parents, dictionary as dictionary_parents};
use crate::python::PyCifDocument;
use crate::reduce::load_monomer_library;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyModule;
use std::collections::HashMap;
use std::path::PathBuf;
use std::sync::Arc;

#[pyfunction]
#[pyo3(signature = (document, *, chem_data=None, model=1))]
fn reduce_hydrogen_parents(
    py: Python<'_>,
    document: &PyCifDocument,
    chem_data: Option<PathBuf>,
    model: i64,
) -> PyResult<HashMap<String, Parents>> {
    let source = Arc::clone(document.document());
    py.detach(move || {
        let library = load_monomer_library(chem_data.as_deref()).map_err(|error| error.to_string())?;
        let mut result: HashMap<String, Parents> = HashMap::new();
        for (component, heavy_names) in observations::ligand_groups(&source, model)? {
            let Some(dictionary) = reduce3::hplace::hydrogen_dictionary(&library, &component, &heavy_names, false) else {
                continue;
            };
            let parents = dictionary_parents(&library, &dictionary)?;
            if let Some(previous) = result.get(&component) {
                if previous != parents.as_ref() {
                    return Err(format!("Reduce3 selected conflicting hydrogen dictionaries for component {component}"));
                }
            } else {
                result.insert(component, parents.as_ref().clone());
            }
        }
        Ok(result)
    }).map_err(PyValueError::new_err)
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(reduce_hydrogen_parents, module)?)?;
    module.add_function(wrap_pyfunction!(reduce_required_parent_sources, module)?)?;
    module.add_function(wrap_pyfunction!(reduce_chemical_sources, module)?)?;
    Ok(())
}

#[pyfunction]
#[pyo3(signature = (document, *, chem_data, model=1))]
fn reduce_chemical_sources(
    py: Python<'_>,
    document: &PyCifDocument,
    chem_data: PathBuf,
    model: i64,
) -> PyResult<Vec<PathBuf>> {
    py.detach(|| observations::chemical_sources(document.document(), model, &chem_data))
        .map_err(PyValueError::new_err)
}

#[pyfunction]
#[pyo3(signature = (document, expected_keys, *, model=1, chem_data=None))]
fn reduce_required_parent_sources(
    py: Python<'_>,
    document: &PyCifDocument,
    expected_keys: &str,
    model: i64,
    chem_data: Option<PathBuf>,
) -> PyResult<Vec<u32>> {
    let source = Arc::clone(document.document());
    py.detach(|| {
        observations::required_parent_sources(&source, expected_keys, model, chem_data.as_deref())
    })
    .map_err(PyValueError::new_err)
}
