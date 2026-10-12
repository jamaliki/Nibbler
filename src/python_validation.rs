//! Python entry point for native model validation.

use std::path::PathBuf;
use std::sync::Arc;

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyModule;

use crate::python::PyCifDocument;
use crate::validation::{Options, load_library, run};

#[pyfunction]
#[pyo3(signature = (document, chem_data, nuclear=false, use_cdl=true, sigma_cutoff=5.0, include_graph=false))]
fn validate_geometry_document(
    py: Python<'_>,
    document: &PyCifDocument,
    chem_data: PathBuf,
    nuclear: bool,
    use_cdl: bool,
    sigma_cutoff: f64,
    include_graph: bool,
) -> PyResult<String> {
    let source = Arc::clone(document.document());
    py.detach(move || {
        let library = load_library(&chem_data)?;
        let result = run(
            &source,
            &library,
            &Options {
                nuclear,
                use_cdl,
                sigma_cutoff,
                include_graph,
            },
        )?;
        serde_json::to_string(&result).map_err(|error| error.to_string())
    })
    .map_err(PyValueError::new_err)
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(validate_geometry_document, module)?)?;
    Ok(())
}
