//! Python receives typed buffers after native observation packing.

use crate::observations::{self, HydrogenParents, LigandFeatures, Schema, Tokens};
use crate::python::PyCifDocument;
use pyo3::buffer::{PyBuffer, PyUntypedBuffer};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::{PyByteArray, PyDict, PyMemoryView, PyModule};
use std::collections::HashMap;
use std::path::PathBuf;
use std::sync::Arc;

mod prepared;
use prepared::PreparedObservationSource;

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

    #[pyo3(signature = (payload, *, min_chain_residues=50, crop_size=384))]
    fn afdb_eligibility<'py>(
        &self,
        py: Python<'py>,
        payload: &Bound<'py, PyAny>,
        min_chain_residues: usize,
        crop_size: usize,
    ) -> PyResult<Bound<'py, PyDict>> {
        let input = PyUntypedBuffer::get(payload)?;
        if !input.readonly()
            || input.dimensions() != 1
            || !input.is_c_contiguous()
            || input.item_size() != 1
        {
            return Err(PyValueError::new_err(
                "AFDB eligibility requires a read-only contiguous byte buffer",
            ));
        }
        let view = PyMemoryView::from(payload)?.call_method1("cast", ("B",))?;
        let buffer = PyBuffer::<u8>::get(&view)?;
        let bytes = buffer
            .as_slice(py)
            .ok_or_else(|| PyValueError::new_err("Invalid AFDB byte buffer"))?;
        let value = observations::afdb_eligibility::scan(
            &self.schema,
            bytes.len(),
            |i| bytes[i].get(),
            min_chain_residues,
            crop_size,
        )
        .map_err(PyValueError::new_err)?;
        let count = value.eligible.len();
        let result = PyDict::new(py);
        result.set_item(
            "eligible",
            ("?", vec![count], PyByteArray::new(py, &value.eligible)),
        )?;
        for (name, values) in [
            ("complete_residues", value.complete_residues),
            ("max_chain_residues", value.max_chain_residues),
            ("reconstructed_chains", value.reconstructed_chains),
        ] {
            let bytes: Vec<_> = values
                .iter()
                .flat_map(|value| value.to_ne_bytes())
                .collect();
            result.set_item(name, ("=i8", vec![count], PyByteArray::new(py, &bytes)))?;
        }
        Ok(result)
    }

    #[staticmethod]
    fn crop_permutation(
        py: Python<'_>,
        before: &[u8],
        after: &[u8],
        chains: &[u8],
        ligand: &[u8],
    ) -> PyResult<Vec<i64>> {
        py.detach(|| observations::crops::permutation(before, after, chains, ligand))
            .map_err(PyValueError::new_err)
    }

    #[staticmethod]
    #[pyo3(signature = (positions, chains, ligand, usable, *, crop_size=384, centre_eligible=None, chain_minimum=None, chain_protein=None))]
    fn crop_plans<'py>(
        py: Python<'py>,
        positions: &[u8],
        chains: &[u8],
        ligand: &[u8],
        usable: &[u8],
        crop_size: usize,
        centre_eligible: Option<&[u8]>,
        chain_minimum: Option<&[u8]>,
        chain_protein: Option<&[u8]>,
    ) -> PyResult<(Bound<'py, PyDict>, f64)> {
        let value = py
            .detach(|| observations::crops::plan_with_chain_policy(positions, chains, ligand, usable, crop_size,
                centre_eligible, chain_minimum, chain_protein))
            .map_err(PyValueError::new_err)?;
        let result = PyDict::new(py);
        let count = value.centers.len();
        let put = |name, dtype, shape: Vec<usize>, bytes: Vec<u8>| -> PyResult<()> {
            result.set_item(name, (dtype, shape, PyByteArray::new(py, &bytes)))
        };
        for (name, values) in [
            ("crop_center", value.centers),
            ("crop_tie", value.ties),
            ("crop_chain_words", vec![value.words as i64]),
            ("crop_offsets", vec![0, count as i64]),
        ] {
            put(
                name,
                "=i8",
                vec![values.len()],
                values
                    .iter()
                    .flat_map(|value| value.to_ne_bytes())
                    .collect(),
            )?;
        }
        put(
            "crop_radius",
            "=f4",
            vec![count],
            value
                .radii
                .iter()
                .flat_map(|value| value.to_ne_bytes())
                .collect(),
        )?;
        put(
            "crop_chain_bits",
            "=u8",
            vec![value.bits.len()],
            value
                .bits
                .iter()
                .flat_map(|value| value.to_ne_bytes())
                .collect(),
        )?;
        put(
            "crop_probability",
            "=f8",
            vec![count],
            value
                .probabilities
                .iter()
                .flat_map(|value| value.to_ne_bytes())
                .collect(),
        )?;
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

    #[cfg(feature = "reduce3")]
    #[pyo3(signature = (document, *, chem_data, model=1))]
    fn missing_hydrogen_dictionaries(&self, py: Python<'_>, document: &PyCifDocument, chem_data: PathBuf, model: i64) -> PyResult<Vec<String>> {
        py.detach(|| observations::missing_hydrogen_dictionaries(document.document(), model, &chem_data))
            .map_err(PyValueError::new_err)
    }

    #[pyo3(signature = (document, *, model=1))]
    fn prepare_source(&self, py: Python<'_>, document: &PyCifDocument, model: i64) -> PyResult<PreparedObservationSource> {
        let source = py.detach(|| observations::PreparedSource::read(document.document(), model))
            .map_err(PyValueError::new_err)?;
        Ok(PreparedObservationSource { source })
    }

    #[pyo3(signature = (document, *, model=1, complete=false, topology_policy="observed_connectivity", protein_only=false, complete_tokens=false, expected_keys=None, validate=false, with_keys=true, chem_data=None, target_uncertainties=None, product_targets=None))]
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
        chem_data: Option<PathBuf>,
        target_uncertainties: Option<&str>,
        product_targets: Option<&str>,
    ) -> PyResult<(Bound<'py, PyDict>, Vec<[String; 5]>)> {
        self.finish_pack(py, || observations::pack(document.document(), &self.schema, &self.features,
            &self.hydrogen_parents, &self.elements, &self.excluded, &observations::PackOptions {
                model, complete, topology_policy, protein_only, chem_data: chem_data.as_deref(),
            }), protein_only, complete_tokens, expected_keys, validate, with_keys, target_uncertainties, product_targets)
    }

    #[pyo3(signature = (source, *, complete=false, topology_policy="observed_connectivity", protein_only=false, complete_tokens=false, expected_keys=None, validate=false, with_keys=true, chem_data=None, target_uncertainties=None, product_targets=None))]
    #[allow(clippy::too_many_arguments)]
    fn pack_prepared<'py>(
        &self, py: Python<'py>, source: &PreparedObservationSource, complete: bool,
        topology_policy: &str, protein_only: bool, complete_tokens: bool,
        expected_keys: Option<&str>, validate: bool, with_keys: bool,
        chem_data: Option<PathBuf>, target_uncertainties: Option<&str>,
        product_targets: Option<&str>,
    ) -> PyResult<(Bound<'py, PyDict>, Vec<[String; 5]>)> {
        self.finish_pack(py, || source.source.pack(&self.schema, &self.features,
            &self.hydrogen_parents, &self.elements, &self.excluded, &observations::PackOptions {
                model: source.source.model, complete, topology_policy, protein_only,
                chem_data: chem_data.as_deref(),
            }), protein_only, complete_tokens, expected_keys, validate, with_keys, target_uncertainties, product_targets)
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

impl ObservationPacker {
    #[allow(clippy::too_many_arguments)]
    fn finish_pack<'py>(&self, py: Python<'py>, build: impl FnOnce() -> Result<Tokens, String> + Send,
        protein_only: bool, complete_tokens: bool, expected_keys: Option<&str>, validate: bool,
        with_keys: bool, target_uncertainties: Option<&str>, product_targets: Option<&str>) -> PyResult<(Bound<'py, PyDict>, Vec<[String; 5]>)> {
        let result = py.detach(|| {
            let mut result = build()?;
            if let Some(json) = product_targets {
                let values: Vec<observations::product::ProductTarget> =
                    serde_json::from_str(json).map_err(|error| error.to_string())?;
                observations::product::apply(&mut result, &self.schema, &values)?;
            }
            if let Some(json) = target_uncertainties {
                let values: Vec<observations::uncertainty::TargetUncertainty> =
                    serde_json::from_str(json).map_err(|error| error.to_string())?;
                observations::uncertainty::apply(&mut result, &self.schema, &values)?;
            }
            if complete_tokens { result.complete(&self.schema)?; }
            if let Some(keys) = expected_keys { result.align(keys)?; }
            result.protein_gaps();
            if validate { result.validate(&self.schema)?; }
            Ok::<_, String>(result)
        }).map_err(PyValueError::new_err)?;
        let arrays = buffers(py, &result, protein_only)?;
        Ok((arrays, if with_keys { result.keys } else { Vec::new() }))
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
        if !value.bridge_edges.is_empty() {
            let edges: Vec<i64> = value.bridge_edges.iter().map(|&(a, _, _)| a as i64)
                .chain(value.bridge_edges.iter().map(|&(_, b, _)| b as i64)).collect();
            put("ligand_bridge_index", "=i8", vec![2, value.bridge_edges.len()], integer(&edges))?;
            put("ligand_bridge_type", "u1", vec![value.bridge_edges.len()],
                value.bridge_edges.iter().map(|&(_, _, kind)| kind).collect())?;
            put("crop_positions", "=f4", vec![count, 3], float(&value.crop_positions))?;
            put("crop_center_eligible", "?", vec![count], value.crop_center_eligible.clone())?;
        }
    }
    Ok(result)
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PreparedObservationSource>()?;
    module.add_class::<ObservationPacker>()
}
