//! Python holds one frozen source while the packer loads component features.

use crate::observations;
#[cfg(feature = "reduce3")]
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
#[cfg(feature = "reduce3")]
use std::path::PathBuf;

#[pyclass(module = "nibbler.observations", frozen)]
pub(super) struct PreparedObservationSource {
    pub(super) source: observations::PreparedSource,
}

#[pymethods]
impl PreparedObservationSource {
    fn ligand_components(&self, py: Python<'_>) -> Vec<String> {
        py.detach(|| self.source.components())
    }

    #[cfg(feature = "reduce3")]
    #[pyo3(signature = (*, chem_data))]
    fn missing_hydrogen_dictionaries(&self, py: Python<'_>, chem_data: PathBuf) -> PyResult<Vec<String>> {
        py.detach(|| self.source.missing_hydrogen_dictionaries(&chem_data)).map_err(PyValueError::new_err)
    }
}
