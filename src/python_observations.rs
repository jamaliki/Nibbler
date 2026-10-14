//! Python receives typed buffers after native observation packing.

use crate::observations::{self, HydrogenParents, LigandFeatures, Schema, Tokens};
use crate::python::PyCifDocument;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::{PyByteArray, PyDict, PyModule};
use std::collections::HashMap;
use std::sync::Arc;

#[pyclass(module = "nibbler.observations")]
struct ObservationPacker {
    schema: Schema,
    elements: Vec<String>,
    excluded: Vec<String>,
    features: LigandFeatures,
    hydrogen_parents: HydrogenParents,
}

#[pymethods]
impl ObservationPacker {
    #[new]
    fn new(schema: &str, elements: Vec<String>, excluded: Vec<String>) -> PyResult<Self> {
        Ok(Self {
            schema: Schema::read(schema).map_err(PyValueError::new_err)?,
            elements,
            excluded,
            features: HashMap::new(),
            hydrogen_parents: HashMap::new(),
        })
    }

    fn add_component(&mut self, component: String, features: HashMap<String, [i16; 4]>) {
        self.features.insert(component, features);
    }

    fn add_hydrogen_parents(&mut self, component: String, parents: HashMap<String, String>) {
        self.hydrogen_parents.insert(component, parents);
    }

    #[staticmethod]
    fn crop_permutation(
        py: Python<'_>, before: &[u8], after: &[u8], chains: &[u8], ligand: &[u8],
    ) -> PyResult<Vec<i64>> {
        py.detach(|| observations::crops::permutation(before, after, chains, ligand))
            .map_err(PyValueError::new_err)
    }

    #[staticmethod]
    #[pyo3(signature = (positions, chains, ligand, usable, *, crop_size=384))]
    fn crop_plans<'py>(
        py: Python<'py>, positions: &[u8], chains: &[u8], ligand: &[u8], usable: &[u8], crop_size: usize,
    ) -> PyResult<(Bound<'py, PyDict>, f64)> {
        let value = py.detach(|| observations::crops::plan(positions, chains, ligand, usable, crop_size))
            .map_err(PyValueError::new_err)?;
        let result = PyDict::new(py);
        let count = value.centers.len();
        let put = |name, dtype, shape: Vec<usize>, bytes: Vec<u8>| -> PyResult<()> {
            result.set_item(name, (dtype, shape, PyByteArray::new(py, &bytes)))
        };
        for (name, values) in [("crop_center", value.centers), ("crop_tie", value.ties),
            ("crop_chain_words", vec![value.words as i64]), ("crop_offsets", vec![0, count as i64])] {
            put(name, "=i8", vec![values.len()], values.iter().flat_map(|value| value.to_ne_bytes()).collect())?;
        }
        put("crop_radius", "=f4", vec![count], value.radii.iter().flat_map(|value| value.to_ne_bytes()).collect())?;
        put("crop_chain_bits", "=u8", vec![value.bits.len()], value.bits.iter().flat_map(|value| value.to_ne_bytes()).collect())?;
        put("crop_probability", "=f8", vec![count], value.probabilities.iter().flat_map(|value| value.to_ne_bytes()).collect())?;
        Ok((result, value.mass))
    }

    #[pyo3(signature = (document, *, model=1))]
    fn ligand_components(
        &self,
        py: Python<'_>,
        document: &PyCifDocument,
        model: i64,
    ) -> PyResult<Vec<String>> {
        py.detach(|| observations::components(document.document(), model))
            .map_err(PyValueError::new_err)
    }

    #[pyo3(signature = (document, *, model=1, complete=false, topology_policy="observed_connectivity", protein_only=false, complete_tokens=false, expected_keys=None, validate=false, with_keys=true))]
    #[allow(clippy::too_many_arguments)]
    fn pack<'py>(
        &self,
        py: Python<'py>,
        document: &PyCifDocument,
        model: i64,
        complete: bool,
        topology_policy: &str,
        protein_only: bool,
        complete_tokens: bool,
        expected_keys: Option<&str>,
        validate: bool,
        with_keys: bool,
    ) -> PyResult<(Bound<'py, PyDict>, Vec<[String; 5]>)> {
        let result = py
            .detach(|| {
                let mut result = observations::pack(
                    document.document(),
                    &self.schema,
                    &self.features,
                    &self.hydrogen_parents,
                    &self.elements,
                    &self.excluded,
                    &observations::PackOptions {
                        model,
                        complete,
                        topology_policy,
                        protein_only,
                    },
                )?;
                if complete_tokens {
                    result.complete(&self.schema);
                }
                if let Some(keys) = expected_keys {
                    result.align(keys)?;
                }
                if validate {
                    result.validate(&self.schema)?;
                }
                Ok::<_, String>(result)
            })
            .map_err(PyValueError::new_err)?;
        let arrays = buffers(py, &result, protein_only)?;
        Ok((arrays, if with_keys { result.keys } else { Vec::new() }))
    }

    fn protein_document(
        &self,
        py: Python<'_>,
        positions: &[u8],
        masks: &[u8],
        identities: &[u8],
        identifier: &str,
    ) -> PyResult<(PyCifDocument, usize, String)> {
        let (document, chains, sequence) = py
            .detach(|| {
                observations::afdb::document(&self.schema, positions, masks, identities, identifier)
            })
            .map_err(PyValueError::new_err)?;
        Ok((
            PyCifDocument::new(Arc::new(document), None),
            chains,
            sequence,
        ))
    }
}

fn buffers<'py>(
    py: Python<'py>,
    value: &Tokens,
    protein_only: bool,
) -> PyResult<Bound<'py, PyDict>> {
    let result = PyDict::new(py);
    let count = value.identity.len();
    let put = |name: &str, dtype: &str, shape: Vec<usize>, bytes: Vec<u8>| -> PyResult<()> {
        result.set_item(name, (dtype, shape, PyByteArray::new(py, &bytes)))
    };
    let integer = |values: &[i64]| {
        values
            .iter()
            .flat_map(|value| value.to_ne_bytes())
            .collect()
    };
    let float = |values: &[f32]| {
        values
            .iter()
            .flat_map(|value| value.to_ne_bytes())
            .collect()
    };
    put("pos", "=f4", vec![count, 32, 3], float(&value.positions))?;
    put("mask", "?", vec![count, 32], value.masks.clone())?;
    put(
        "atom_presence_known",
        "?",
        vec![count, 32],
        value.presence.clone(),
    )?;
    put("topology", "=i8", vec![count, 3], integer(&value.topology))?;
    for (name, values) in [
        ("aatype", &value.identity),
        ("res_idx", &value.sequence),
        ("chain_idx", &value.chains),
        ("chemical_state", &value.state),
    ] {
        put(name, "=i8", vec![count], integer(values))?;
    }
    put("state_known", "?", vec![count], value.state_known.clone())?;
    if !protein_only {
        for (name, values) in [
            ("is_ligand", &value.ligand),
            ("ligand_charge_known", &value.charge_known),
            ("ligand_hydrogen_known", &value.hydrogen_known),
            ("covalently_bound_to_protein", &value.covalent),
        ] {
            put(name, "?", vec![count], values.clone())?;
        }
        put("ligand_charge", "=f4", vec![count], float(&value.charges))?;
        put(
            "ligand_hybridization",
            "=i8",
            vec![count],
            integer(&value.hybridization),
        )?;
        let edges: Vec<i64> = value
            .attachments
            .iter()
            .map(|&(a, _)| a as i64)
            .chain(value.attachments.iter().map(|&(_, b)| b as i64))
            .collect();
        put(
            "protein_attachment_index",
            "=i8",
            vec![2, value.attachments.len()],
            integer(&edges),
        )?;
    }
    Ok(result)
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<ObservationPacker>()
}
