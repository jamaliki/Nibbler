//! The `ConformerSelection` class and conformer selection from models, documents and arrays.

use std::sync::Arc;

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyModule;

use super::{PyMmcifModel, model_error_to_python};
use crate::cif::CifDocument;
use crate::pdbx::PdbxModel;
use crate::python::document::PyCifDocument;

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
                    document.document(),
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
        source: Some(Arc::clone(document.document())),
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

pub(super) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyConformerSelection>()?;
    module.add_function(wrap_pyfunction!(select_mmcif_conformer, module)?)?;
    module.add_function(wrap_pyfunction!(select_document_conformer, module)?)?;
    module.add_function(wrap_pyfunction!(select_conformer_arrays, module)?)?;
    Ok(())
}
