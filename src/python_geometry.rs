//! Python bindings for explicit native restraint evaluation.

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList, PyModule};

use crate::geometry::{Angle, Bond, Chirality, Plane, RestraintGraph};
use crate::python::PyCifDocument;

#[pyfunction]
#[pyo3(signature = (coordinates, radii, atom_indices, probe_radius=1.4, point_count=256))]
/// Return selected atom areas in A^2 from coordinates and van der Waals radii in A.
///
/// All input atoms can occlude the selected atoms. The caller must select models,
/// conformers, waters and hydrogens before supplying arrays. Results follow
/// atom_indices order. More sphere points improve the Shrake-Rupley estimate at
/// greater cost. Native work releases the GIL.
fn atom_sasa(
    py: Python<'_>,
    coordinates: Vec<[f64; 3]>,
    radii: Vec<f64>,
    atom_indices: Vec<usize>,
    probe_radius: f64,
    point_count: usize,
) -> PyResult<Vec<f64>> {
    py.detach(|| {
        crate::geometry::sasa(
            &coordinates,
            &radii,
            &atom_indices,
            probe_radius,
            point_count,
        )
    })
    .map_err(PyValueError::new_err)
}

#[pyfunction]
#[pyo3(signature = (document, atom_rows=None, model_number=1, probe_radius=1.4, point_count=256, radii=None))]
fn ligand_sasa_document(
    py: Python<'_>,
    document: &PyCifDocument,
    atom_rows: Option<Vec<usize>>,
    model_number: i64,
    probe_radius: f64,
    point_count: usize,
    radii: Option<std::collections::HashMap<String, f64>>,
) -> PyResult<(Vec<usize>, Vec<f64>)> {
    let radii = radii.unwrap_or_default();
    py.detach(|| {
        crate::pdbx::ligand_sasa(
            document.document(),
            atom_rows.as_deref(),
            model_number,
            probe_radius,
            point_count,
            &radii,
        )
    })
    .map_err(PyValueError::new_err)
}

#[pyclass(name = "RestraintGraph", module = "nibbler.geometry", frozen)]
struct PyRestraintGraph {
    graph: RestraintGraph,
}

#[pymethods]
impl PyRestraintGraph {
    #[new]
    fn new(
        atom_ids: Vec<String>,
        hydrogen_flags: Vec<bool>,
        bonds: Vec<Bond>,
        angles: Vec<Angle>,
        chiralities: Vec<Chirality>,
        planes: Vec<Plane>,
    ) -> PyResult<Self> {
        Ok(Self {
            graph: RestraintGraph::new(
                atom_ids,
                hydrogen_flags,
                bonds,
                angles,
                chiralities,
                planes,
            )
            .map_err(PyValueError::new_err)?,
        })
    }

    #[pyo3(signature = (coordinates, sigma_cutoff=5.0))]
    fn evaluate<'py>(
        &self,
        py: Python<'py>,
        coordinates: Vec<[f64; 3]>,
        sigma_cutoff: f64,
    ) -> PyResult<Bound<'py, PyDict>> {
        let result = py
            .detach(|| self.graph.evaluate(&coordinates, sigma_cutoff))
            .map_err(PyValueError::new_err)?;
        let output = PyDict::new(py);
        output.set_item("near_cutoff", result.near_cutoff)?;
        let checks = PyDict::new(py);
        for (name, count) in result.counts {
            let check = PyDict::new(py);
            check.set_item("status", "checked")?;
            check.set_item("count", count)?;
            checks.set_item(name, check)?;
        }
        let failures = PyList::empty(py);
        for finding in result.findings {
            let row = PyDict::new(py);
            row.set_item("check", finding.check)?;
            row.set_item("proxy_index", finding.proxy_index)?;
            row.set_item(
                "atom_ids",
                finding
                    .indices
                    .iter()
                    .map(|&i| &self.graph.atom_ids[i])
                    .collect::<Vec<_>>(),
            )?;
            row.set_item("indices", finding.indices)?;
            row.set_item("model", finding.model)?;
            row.set_item("target", finding.target)?;
            row.set_item("delta", finding.delta)?;
            row.set_item("sigma", finding.sigma)?;
            row.set_item("score", finding.score)?;
            row.set_item("residual", finding.residual)?;
            failures.append(row)?;
        }
        for (i, parents) in result.hydrogen_failures {
            let row = PyDict::new(py);
            row.set_item("check", "hydrogens")?;
            row.set_item("atom", &self.graph.atom_ids[i])?;
            row.set_item("indices", vec![i])?;
            row.set_item("parents", parents)?;
            row.set_item("reason", "hydrogen_parent_count")?;
            failures.append(row)?;
        }
        output.set_item("checks", checks)?;
        output.set_item("failures", failures)?;
        Ok(output)
    }
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(atom_sasa, module)?)?;
    module.add_function(wrap_pyfunction!(ligand_sasa_document, module)?)?;
    module.add_class::<PyRestraintGraph>()?;
    Ok(())
}
