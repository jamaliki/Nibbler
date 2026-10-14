//! Hydrogen parents from the same fixed-mode dictionaries as Reduce3.

use crate::observations;
use crate::python::PyCifDocument;
use crate::reduce::load_monomer_library;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyModule;
use reduce3::monlib::{Comp, MonLib};
use std::collections::HashMap;
use std::path::PathBuf;
use std::sync::{Arc, Mutex, OnceLock, PoisonError};

type Parents = HashMap<String, String>;
type CacheKey = (PathBuf, PathBuf, String);
static PARENTS: OnceLock<Mutex<HashMap<CacheKey, Arc<Parents>>>> = OnceLock::new();

fn hydrogen(element: &str) -> bool {
    element.eq_ignore_ascii_case("H")
        || element.eq_ignore_ascii_case("D")
        || element.eq_ignore_ascii_case("T")
}

fn dictionary_parents(library: &MonLib, dictionary: &Comp) -> Result<Arc<Parents>, String> {
    let key = (
        library.root.clone(),
        dictionary.source.clone(),
        dictionary.id.clone(),
    );
    let cache = PARENTS.get_or_init(Mutex::default);
    if let Some(parents) = cache
        .lock()
        .unwrap_or_else(PoisonError::into_inner)
        .get(&key)
    {
        return Ok(Arc::clone(parents));
    }
    let mut parents = HashMap::new();
    for bond in &dictionary.bonds {
        let (Some(first), Some(second)) = (dictionary.atom(&bond.a1), dictionary.atom(&bond.a2))
        else {
            return Err(format!(
                "Reduce3 dictionary {} has a bond with an absent atom",
                dictionary.id
            ));
        };
        let pair = match (hydrogen(&first.type_symbol), hydrogen(&second.type_symbol)) {
            (true, false) => Some((&bond.a1, &bond.a2)),
            (false, true) => Some((&bond.a2, &bond.a1)),
            _ => None,
        };
        if let Some((child, parent)) = pair
            && let Some(previous) = parents.insert(child.clone(), parent.clone())
            && previous != *parent
        {
            return Err(format!(
                "Reduce3 dictionary {} hydrogen {child} has conflicting parents",
                dictionary.id
            ));
        }
    }
    let parents = Arc::new(parents);
    cache
        .lock()
        .unwrap_or_else(PoisonError::into_inner)
        .insert(key, Arc::clone(&parents));
    Ok(parents)
}

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
    Ok(())
}

#[pyfunction]
#[pyo3(signature = (document, expected_keys, *, model=1))]
fn reduce_required_parent_sources(
    py: Python<'_>,
    document: &PyCifDocument,
    expected_keys: &str,
    model: i64,
) -> PyResult<Vec<u32>> {
    let source = Arc::clone(document.document());
    py.detach(|| observations::required_parent_sources(&source, expected_keys, model))
        .map_err(PyValueError::new_err)
}
