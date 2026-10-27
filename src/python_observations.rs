//! Python receives typed buffers after native observation packing.

use crate::observations::{self, HydrogenParents, LigandFeatures, Schema, TokenTable, Tokens};
use crate::python::PyCifDocument;
use pyo3::buffer::{PyBuffer, PyUntypedBuffer};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::{PyByteArray, PyDict, PyMemoryView, PyModule};
use std::collections::HashMap;
use std::path::PathBuf;
use std::sync::Arc;

#[cfg(feature = "reduce3")]
pub(crate) mod molecular;
mod prepared;
use prepared::PreparedObservationSource;

#[pyclass(module = "nibbler.observations")]
struct ObservationPacker {
    schema: Schema,
    table: TokenTable,
    excluded: Vec<String>,
    features: LigandFeatures,
    hydrogen_parents: HydrogenParents,
    #[cfg(feature = "reduce3")]
    molecular: Option<observations::molecular::MolecularSchema>,
}

#[pymethods]
impl ObservationPacker {
    #[new]
    fn new(schema: &str, token_table: &Bound<'_, PyDict>, excluded: Vec<String>) -> PyResult<Self> {
        Ok(Self {
            schema: Schema::read(schema).map_err(PyValueError::new_err)?,
            table: token_table_from(token_table)?,
            excluded,
            features: HashMap::new(),
            hydrogen_parents: HashMap::new(),
            #[cfg(feature = "reduce3")]
            molecular: None,
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
    #[allow(clippy::too_many_arguments)]
    #[pyo3(signature = (positions, chains, ligand, usable, *, crop_size=384, centre_eligible=None, chain_minimum=None, chain_protein=None, ligand_limits=Some(observations::crops::LIGAND_LIMITS)))]
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
        ligand_limits: Option<(u32, u32)>,
    ) -> PyResult<(Bound<'py, PyDict>, f64)> {
        let value = py
            .detach(|| observations::crops::plan_with_chain_policy(positions, chains, ligand, usable, crop_size,
                centre_eligible, chain_minimum, chain_protein, ligand_limits))
            .map_err(PyValueError::new_err)?;
        let result = PyDict::new(py);
        let count = value.centers.len();
        let put = |name, dtype, shape: Vec<usize>, bytes: Vec<u8>| -> PyResult<()> {
            result.set_item(name, (dtype, shape, PyByteArray::new(py, &bytes)))
        };
        let counts = |values: Vec<u32>| values.into_iter().map(i64::from).collect::<Vec<_>>();
        let uncapped = ligand_limits.is_none();
        for (name, values) in [
            ("crop_center", value.centers),
            ("crop_tie", value.ties),
            ("crop_chain_words", vec![value.words as i64]),
            ("crop_offsets", vec![0, count as i64]),
            ("crop_ligand_tokens", counts(value.ligand_tokens)),
            ("crop_ligand_chains", counts(value.ligand_chains)),
        ]
        .into_iter()
        .chain(uncapped.then(|| [
            ("centre_ligand_tokens", counts(value.centre_ligand_tokens)),
            ("centre_ligand_chains", counts(value.centre_ligand_chains)),
        ]).into_iter().flatten())
        {
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

    /// Store the Atom84 contract used by the single-read molecular preparation.
    #[cfg(feature = "reduce3")]
    fn set_molecular_schema(&mut self, schema: &str) -> PyResult<()> {
        self.molecular = Some(observations::molecular::MolecularSchema::read(schema).map_err(PyValueError::new_err)?);
        Ok(())
    }

    /// Read the source table and the linked product of one document once.
    ///
    /// `base` is a view of the document Reduce3 read to write this one. The view then shares
    /// its product when this document gives the same one.
    #[cfg(feature = "reduce3")]
    #[pyo3(signature = (document, *, model=1, chem_data, base=None))]
    fn molecular_view(&self, py: Python<'_>, document: &PyCifDocument, model: i64,
        chem_data: PathBuf, base: Option<PyRef<'_, molecular::MolecularView>>) -> PyResult<molecular::MolecularView> {
        let schema = self.molecular.as_ref()
            .ok_or_else(|| PyValueError::new_err("The packer has no molecular schema"))?;
        // The document keeps its model, so Reduce3 and the model of its output can reuse
        // the decoded atom sites. A model that fails to build is built again by the read,
        // which reports source errors first.
        let source = document;
        let document = source.document().clone();
        let base = base.as_ref().map(|base| &base.view);
        let view = py.detach(|| {
            // A document that Reduce3 wrote from the document of `base` reads its view from
            // the model of that document, without a model of its own.
            #[cfg(feature = "reduce3")]
            let copied = base.and_then(|base| source.source_copies(|copied, rows, made| {
                observations::molecular::View::read_copied(document.clone(), copied, rows, made, model, &chem_data,
                    schema, base)
            })?);
            #[cfg(not(feature = "reduce3"))]
            let copied = None;
            match copied {
                Some(view) => Ok(view),
                None => observations::molecular::View::read_with_model(document, source.pdbx_model().ok(), model,
                    &chem_data, schema, base),
            }
        })
        .map_err(molecular::raise)?;
        Ok(molecular::MolecularView { view })
    }

    /// Build legacy tokens, the source snapshot and typed Atom84 rows from one view.
    ///
    /// With `index`, only the token layout that the heavy index stores and hashes comes back.
    /// With `reference_labels`, the source arrays also hold `source_fixed_component`,
    /// `source_fixed_element` and `source_fixed_name`: UTF-8 labels zero-padded to eight
    /// bytes, as `S8` arrays. They are left out when a label is longer or ends in a NUL.
    #[cfg(feature = "reduce3")]
    #[pyo3(signature = (view, perceived, *, complete=false, sasa=true, target_uncertainties=None, expected_keys=None,
        keys_json=false, index=false, reference_labels=false))]
    #[allow(clippy::too_many_arguments)]
    fn finish_molecular<'py>(&self, py: Python<'py>, view: &molecular::MolecularView, perceived: Vec<Vec<u8>>,
        complete: bool, sasa: bool, target_uncertainties: Option<&str>, expected_keys: Option<&str>, keys_json: bool,
        index: bool, reference_labels: bool) -> PyResult<molecular::FinishedValue<'py>> {
        let schema = self.molecular.as_ref()
            .ok_or_else(|| PyValueError::new_err("The packer has no molecular schema"))?;
        let packer = observations::molecular::Packer {
            schema: &self.schema, features: &self.features, hydrogen_parents: &self.hydrogen_parents,
            table: &self.table, excluded: &self.excluded, molecular: schema,
        };
        let options = observations::molecular::Options {
            complete, sasa, uncertainties: target_uncertainties, expected: expected_keys, keys_json, index,
            reference_labels,
        };
        let value = py.detach(|| observations::molecular::finish(&view.view, &perceived, &packer, &options))
            .map_err(molecular::raise)?;
        molecular::finished(py, value)
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
            &self.hydrogen_parents, &self.table, &self.excluded, &observations::PackOptions {
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
            &self.hydrogen_parents, &self.table, &self.excluded, &observations::PackOptions {
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

/// Read `{"elements": {symbol: token}, "other": token, "hydrogen": token}`.
fn token_table_from(value: &Bound<'_, PyDict>) -> PyResult<TokenTable> {
    let (mut elements, mut other, mut hydrogen) = (None, None, None);
    for (key, item) in value.iter() {
        match key.extract::<String>()?.as_str() {
            "elements" => elements = Some(item.extract::<HashMap<String, i64>>()?),
            "other" => other = Some(item.extract::<i64>()?),
            "hydrogen" => hydrogen = Some(item.extract::<i64>()?),
            name => return Err(PyValueError::new_err(format!("Unknown ligand token table field {name:?}"))),
        }
    }
    let missing = |name: &str| PyValueError::new_err(format!("The ligand token table has no {name}"));
    TokenTable::new(elements.ok_or_else(|| missing("elements"))?, other.ok_or_else(|| missing("other token"))?,
        hydrogen.ok_or_else(|| missing("hydrogen token"))?).map_err(PyValueError::new_err)
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
                observations::product::apply(&mut result, &self.schema, &self.table, &values)?;
            }
            if let Some(json) = target_uncertainties {
                let values: Vec<observations::uncertainty::TargetUncertainty> =
                    serde_json::from_str(json).map_err(|error| error.to_string())?;
                observations::uncertainty::apply(&mut result, &self.schema, &values)?;
            }
            if complete_tokens { result.complete(&self.schema)?; }
            if let Some(keys) = expected_keys { result.align(keys)?; }
            result.protein_gaps();
            if validate { result.validate(&self.schema, &self.table)?; }
            Ok::<_, String>(result)
        }).map_err(PyValueError::new_err)?;
        let arrays = buffers(py, &result, protein_only)?;
        let keys = if with_keys { result.keys.iter().map(|key| key.each_ref().map(String::from)).collect() } else { Vec::new() };
        Ok((arrays, keys))
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
    #[cfg(feature = "reduce3")]
    module.add_class::<molecular::MolecularView>()?;
    module.add_class::<ObservationPacker>()
}
