//! The `PreparedChemistry` class: the linked product graph of one model.

use pyo3::prelude::*;
use pyo3::types::PyModule;

use super::{PyComponentRegistry, PyMmcifModel, model_error_to_python, objects};

#[pyclass(name = "PreparedChemistry", module = "nibbler.mmcif", frozen)]
struct PyPreparedChemistry {
    graph: crate::pdbx::PreparedChemistry,
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
                return Err(crate::pdbx::SemanticError::new(
                    "PDBX_LINK_PRODUCT",
                    "Supply registry or chem_data, not both",
                    Vec::new(),
                ));
            }
            let library = chem_data
                .as_deref()
                .map(|root| {
                    crate::pdbx::chemical_library_registry(model.model.coordinates(), root).map_err(
                        |message| {
                            crate::pdbx::SemanticError::new(
                                "PDBX_LINK_PRODUCT",
                                message,
                                Vec::new(),
                            )
                        },
                    )
                })
                .transpose()?;
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

pub(super) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyPreparedChemistry>()?;
    module.add_function(wrap_pyfunction!(prepare_mmcif_chemistry, module)?)?;
    Ok(())
}
