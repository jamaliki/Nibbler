//! Python holds one frozen source while the packer loads component features.

use crate::observations;
#[cfg(feature = "reduce3")]
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::{PyByteArray, PyDict};
use std::collections::HashMap;
#[cfg(feature = "reduce3")]
use std::path::PathBuf;

#[pyclass(module = "nibbler.observations", frozen)]
pub(super) struct PreparedObservationSource {
    pub(super) source: observations::PreparedSource,
}

#[pymethods]
impl PreparedObservationSource {
    fn atom_table<'py>(&self, py: Python<'py>) -> PyResult<(Bound<'py, PyDict>, HashMap<String, Vec<String>>)> {
        let source = &self.source.source;
        let count = source.atoms.len();
        let buffers = PyDict::new(py);
        let mut parents = vec![-1_i64; count];
        for &(child, parent) in &source.hydrogen_bonds {
            if let Some(parent) = parent {
                parents[child] = if parents[child] == -1 || parents[child] == parent as i64 {
                    parent as i64
                } else { -2 };
            }
        }
        for value in &mut parents { if *value < 0 { *value = -1; } }
        let put = |name, dtype, shape: Vec<usize>, data: Vec<u8>| {
            buffers.set_item(name, (dtype, shape, PyByteArray::new(py, &data)))
        };
        let integers: Vec<(&str, Vec<i64>)> = vec![
            ("source_atom_row", source.atoms.iter().map(|atom| i64::from(atom.source_row)).collect()),
            ("source_atom_group", source.atoms.iter().map(|atom| atom.group as i64).collect()),
            ("source_atom_sequence", source.atoms.iter().map(|atom| atom.sequence).collect()),
            ("source_hydrogen_parent", parents),
        ];
        for (name, values) in integers {
            put(name, "=i8", vec![count], values.iter().flat_map(|v| v.to_ne_bytes()).collect())?;
        }
        put("source_atom_pos", "=f4", vec![count, 3], source.atoms.iter()
            .flat_map(|atom| atom.position).flat_map(|v| v.to_ne_bytes()).collect())?;
        put("source_atom_occupancy", "=f4", vec![count], source.atoms.iter()
            .flat_map(|atom| atom.occupancy.to_ne_bytes()).collect())?;
        for (name, values) in [
            ("source_atom_mask", source.atoms.iter().map(|atom| u8::from(atom.position.iter().all(|v| v.is_finite()))).collect()),
            ("source_atom_is_protein", source.atoms.iter().map(|atom| u8::from(atom.protein)).collect()),
        ] {
            put(name, "?", vec![count], values)?;
        }
        let edges = source.connections.len();
        put("source_connection_index", "=i8", vec![2, edges], (0..2).flat_map(|side|
            source.connections.iter().map(move |edge| edge.atoms[side].map_or(-1, |row| row as i64)))
            .flat_map(|v| v.to_ne_bytes()).collect())?;
        put("source_connection_row", "=i8", vec![edges], source.connections.iter()
            .flat_map(|edge| (edge.row as i64).to_ne_bytes()).collect())?;
        let mut labels = HashMap::new();
        for name in ["name", "element", "component", "chain", "auth_chain", "auth_sequence",
            "insertion", "entity", "polymer_type", "alternate"] {
            labels.insert(format!("source_atom_{name}"), source.atoms.iter().map(|atom| match name {
                "name" => atom.name.clone(), "element" => atom.element.clone(),
                "component" => atom.component.clone(), "chain" => atom.chain.clone(),
                "auth_chain" => atom.auth_chain.clone(), "auth_sequence" => atom.auth_sequence.clone(),
                "insertion" => atom.insertion.clone(), "entity" => atom.entity.clone(),
                "alternate" => atom.alternate.clone(),
                _ => source.entity_types.get(&atom.entity).cloned().unwrap_or_default(),
            }).collect());
        }
        for name in ["type", "order", "symmetry1", "symmetry2"] {
            labels.insert(format!("source_connection_{name}"), source.connections.iter().map(|edge| match name {
                "type" => edge.kind.clone(), "order" => edge.order.clone(),
                "symmetry1" => edge.symmetries[0].clone(), _ => edge.symmetries[1].clone(),
            }).collect());
        }
        Ok((buffers, labels))
    }

    fn ligand_components(&self, py: Python<'_>) -> Vec<String> {
        py.detach(|| self.source.components())
    }

    #[cfg(feature = "reduce3")]
    #[pyo3(signature = (*, chem_data))]
    fn missing_hydrogen_dictionaries(&self, py: Python<'_>, chem_data: PathBuf) -> PyResult<Vec<String>> {
        py.detach(|| self.source.missing_hydrogen_dictionaries(&chem_data)).map_err(PyValueError::new_err)
    }
}
