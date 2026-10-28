//! Python perceives ligand domains between one native read and one native finish.

use crate::observations::molecular::{self, Failure, View};
use pyo3::IntoPyObjectExt;
use pyo3::exceptions::{PyIndexError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyByteArray, PyDict, PyList, PyString};
use std::path::PathBuf;

/// Raise the Python exception class of the step that failed.
pub(super) fn raise(error: Failure) -> PyErr {
    match error {
        Failure::Value(message) => PyValueError::new_err(message),
        Failure::Chemistry(error) => PyValueError::new_err((
            error.code().to_owned(),
            error.message().to_owned(),
            error.context().to_vec(),
        )),
        Failure::Index(message) => PyIndexError::new_err(message),
    }
}

#[pyclass(module = "nibbler.observations", frozen)]
pub(crate) struct MolecularView {
    pub(crate) view: View,
}

#[pymethods]
impl MolecularView {
    /// Atom and bond dictionaries of each bounded ligand domain, caps last.
    fn perception_graphs<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyList>> {
        let result = PyList::empty(py);
        for (atoms, bonds) in self.view.perception_graphs() {
            let records = PyList::empty(py);
            for (index, item) in atoms.iter().enumerate() {
                let atom = item.atom;
                let record = PyDict::new(py);
                record.set_item("asym_id", atom.asym_id.as_str())?;
                record.set_item("auth_sequence_id", atom.auth_sequence_id.as_str())?;
                record.set_item("insertion_code", atom.insertion_code.as_deref())?;
                record.set_item("component_id", atom.component_id.as_str())?;
                if item.cap {
                    let name = format!("boundary-H-{index}");
                    record.set_item("element", "H")?;
                    record.set_item("atom_id", &name)?;
                    record.set_item("observed_atom_id", &name)?;
                    record.set_item("formal_charge", 0)?;
                    record.set_item("hydrogen_count", 0)?;
                    record.set_item("active", true)?;
                } else {
                    record.set_item("element", atom.element.as_str())?;
                    record.set_item("atom_id", atom.atom_id.as_str())?;
                    record.set_item("observed_atom_id", item.observed)?;
                    record.set_item("formal_charge", atom.formal_charge)?;
                    record.set_item("hydrogen_count", atom.hydrogen_count)?;
                    record.set_item("active", atom.active)?;
                }
                record.set_item("source_row", item.source_row)?;
                records.append(record)?;
            }
            let edges = PyList::empty(py);
            for &(first, second, order) in bonds {
                let record = PyDict::new(py);
                record.set_item("active", true)?;
                record.set_item("order", order)?;
                record.set_item("atom_indices", [first, second])?;
                edges.append(record)?;
            }
            let graph = PyDict::new(py);
            graph.set_item("atoms", records)?;
            graph.set_item("bonds", edges)?;
            result.append(graph)?;
        }
        Ok(result)
    }

    fn ligand_components(&self, py: Python<'_>) -> Vec<String> {
        py.detach(|| self.view.ligand_components())
    }

    #[pyo3(signature = (*, chem_data))]
    fn missing_hydrogen_dictionaries(
        &self,
        py: Python<'_>,
        chem_data: PathBuf,
    ) -> PyResult<Vec<String>> {
        py.detach(|| self.view.missing_hydrogen_dictionaries(&chem_data))
            .map_err(PyValueError::new_err)
    }
}

/// Typed buffers, token keys, retained legacy keys and audit counts.
pub(super) type FinishedValue<'py> = (
    Bound<'py, PyDict>,
    Bound<'py, PyAny>,
    Bound<'py, PyAny>,
    Bound<'py, PyDict>,
);

/// One key list as JSON text or as nested Python lists.
fn keys<'py>(py: Python<'py>, keys: molecular::Keys) -> PyResult<Bound<'py, PyAny>> {
    match keys {
        molecular::Keys::Text(text) => Ok(PyString::new(py, &text).into_any()),
        molecular::Keys::Lists(lists) => lists.into_bound_py_any(py),
    }
}

/// Convert one finished pack; key lists come back as the finish options chose.
pub(super) fn finished<'py>(
    py: Python<'py>,
    value: molecular::Finished,
) -> PyResult<FinishedValue<'py>> {
    let (keys, required) = (keys(py, value.keys)?, keys(py, value.required)?);
    let buffers = PyDict::new(py);
    for (name, array) in value.arrays {
        buffers.set_item(
            name,
            (array.dtype, array.shape, PyByteArray::new(py, &array.bytes)),
        )?;
    }
    let stats = value.stats;
    let audit = PyDict::new(py);
    audit.set_item("tokens", stats.tokens)?;
    audit.set_item("source_atoms", stats.source_atoms)?;
    audit.set_item(
        "unsupported_known_charge_targets",
        stats.unsupported_known_charge_targets,
    )?;
    audit.set_item("unknown_sequence_tokens", stats.unknown_sequence_tokens)?;
    audit.set_item("represented_atoms", stats.represented_atoms)?;
    let exclusions = PyDict::new(py);
    for (name, count) in stats.exclusions {
        exclusions.set_item(name, count)?;
    }
    audit.set_item("exclusions", exclusions)?;
    audit.set_item("nucleotide_residues", stats.nucleotide_residues)?;
    audit.set_item(
        "unknown_source_atom_tokens",
        stats.unknown_source_atom_tokens,
    )?;
    audit.set_item("ncaa_sidechain_tokens", stats.ncaa_sidechain_tokens)?;
    audit.set_item(
        "modified_nucleotide_atom_tokens",
        stats.modified_nucleotide_atom_tokens,
    )?;
    audit.set_item(
        "fully_atomized_modified_nucleotides",
        stats.fully_atomized_modified_nucleotides,
    )?;
    Ok((buffers, keys, required, audit))
}
