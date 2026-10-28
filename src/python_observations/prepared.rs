//! Python holds one frozen source while the packer loads component features.

use crate::cif::FastMap;
use crate::observations;
use crate::pdbx::ShortText;
#[cfg(feature = "reduce3")]
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::{PyByteArray, PyDict, PyList, PyString};
#[cfg(feature = "reduce3")]
use std::path::PathBuf;

#[pyclass(module = "nibbler.observations", frozen)]
pub(super) struct PreparedObservationSource {
    pub(super) source: observations::PreparedSource,
}

#[pymethods]
impl PreparedObservationSource {
    fn atom_table<'py>(
        &self,
        py: Python<'py>,
    ) -> PyResult<(Bound<'py, PyDict>, Bound<'py, PyDict>)> {
        let source = &self.source.source;
        let count = source.atoms.len();
        let buffers = PyDict::new(py);
        let mut parents = vec![-1_i64; count];
        for &(child, parent) in &source.hydrogen_bonds {
            if let Some(parent) = parent {
                parents[child] = if parents[child] == -1 || parents[child] == parent as i64 {
                    parent as i64
                } else {
                    -2
                };
            }
        }
        for value in &mut parents {
            if *value < 0 {
                *value = -1;
            }
        }
        let put = |name, dtype, shape: Vec<usize>, data: Vec<u8>| {
            buffers.set_item(name, (dtype, shape, PyByteArray::new(py, &data)))
        };
        let integers: [(&str, Vec<i64>); 4] = [
            (
                "source_atom_row",
                source
                    .atoms
                    .iter()
                    .map(|atom| i64::from(atom.source_row))
                    .collect(),
            ),
            (
                "source_atom_group",
                source.atoms.iter().map(|atom| atom.group as i64).collect(),
            ),
            (
                "source_atom_sequence",
                source.atoms.iter().map(|atom| atom.sequence).collect(),
            ),
            ("source_hydrogen_parent", parents),
        ];
        for (name, values) in integers {
            put(
                name,
                "=i8",
                vec![count],
                bytes(count, values.iter().map(|v| v.to_ne_bytes())),
            )?;
        }
        put(
            "source_atom_pos",
            "=f4",
            vec![count, 3],
            bytes(
                3 * count,
                source
                    .atoms
                    .iter()
                    .flat_map(|atom| atom.position)
                    .map(f32::to_ne_bytes),
            ),
        )?;
        put(
            "source_atom_occupancy",
            "=f4",
            vec![count],
            bytes(
                count,
                source.atoms.iter().map(|atom| atom.occupancy.to_ne_bytes()),
            ),
        )?;
        for (name, values) in [
            (
                "source_atom_mask",
                source
                    .atoms
                    .iter()
                    .map(|atom| u8::from(atom.position.iter().all(|v| v.is_finite())))
                    .collect(),
            ),
            (
                "source_atom_is_protein",
                source
                    .atoms
                    .iter()
                    .map(|atom| u8::from(atom.protein))
                    .collect(),
            ),
        ] {
            put(name, "?", vec![count], values)?;
        }
        let edges = source.connections.len();
        put(
            "source_connection_index",
            "=i8",
            vec![2, edges],
            bytes(
                2 * edges,
                (0..2)
                    .flat_map(|side| {
                        source
                            .connections
                            .iter()
                            .map(move |edge| edge.atoms[side].map_or(-1, |row| row as i64))
                    })
                    .map(i64::to_ne_bytes),
            ),
        )?;
        put(
            "source_connection_row",
            "=i8",
            vec![edges],
            bytes(
                edges,
                source
                    .connections
                    .iter()
                    .map(|edge| (edge.row as i64).to_ne_bytes()),
            ),
        )?;
        // Labels repeat heavily, so each distinct value becomes one Python string that
        // the lists share.
        let labels = PyDict::new(py);
        let mut strings: FastMap<&str, Bound<'py, PyString>> = FastMap::default();
        let atoms = &source.atoms;
        let atom_columns: [(&str, AtomLabel); 9] = [
            ("source_atom_name", |atom| &atom.name),
            ("source_atom_element", |atom| &atom.element),
            ("source_atom_component", |atom| &atom.component),
            ("source_atom_chain", |atom| &atom.chain),
            ("source_atom_auth_chain", |atom| &atom.auth_chain),
            ("source_atom_auth_sequence", |atom| &atom.auth_sequence),
            ("source_atom_insertion", |atom| &atom.insertion),
            ("source_atom_entity", |atom| &atom.entity),
            ("source_atom_alternate", |atom| &atom.alternate),
        ];
        for (name, field) in atom_columns {
            let items = shared_strings(py, &mut strings, atoms.iter().map(field));
            labels.set_item(name, PyList::new(py, items)?)?;
        }
        // Atoms of one entity are adjacent, so one lookup serves each run.
        let mut entity: Option<(&ShortText, &str)> = None;
        let polymer_types = atoms.iter().map(|atom| match entity {
            Some((previous, kind)) if *previous == atom.entity => kind,
            _ => {
                let kind = source
                    .entity_types
                    .get(atom.entity.as_str())
                    .map_or("", String::as_str);
                entity = Some((&atom.entity, kind));
                kind
            }
        });
        let items = shared_strings(py, &mut strings, polymer_types);
        labels.set_item("source_atom_polymer_type", PyList::new(py, items)?)?;
        let connections = &source.connections;
        let connection_columns: [(&str, ConnectionLabel); 4] = [
            ("source_connection_type", |edge| &edge.kind),
            ("source_connection_order", |edge| &edge.order),
            ("source_connection_symmetry1", |edge| &edge.symmetries[0]),
            ("source_connection_symmetry2", |edge| &edge.symmetries[1]),
        ];
        for (name, field) in connection_columns {
            let items = shared_strings(py, &mut strings, connections.iter().map(field));
            labels.set_item(name, PyList::new(py, items)?)?;
        }
        Ok((buffers, labels))
    }

    fn ligand_components(&self, py: Python<'_>) -> Vec<String> {
        py.detach(|| self.source.components())
    }

    #[cfg(feature = "reduce3")]
    #[pyo3(signature = (*, chem_data))]
    fn missing_hydrogen_dictionaries(
        &self,
        py: Python<'_>,
        chem_data: PathBuf,
    ) -> PyResult<Vec<String>> {
        py.detach(|| self.source.missing_hydrogen_dictionaries(&chem_data))
            .map_err(PyValueError::new_err)
    }
}

/// Read one label of an atom.
type AtomLabel = fn(&observations::source::Atom) -> &ShortText;
/// Read one label of a connection.
type ConnectionLabel = fn(&observations::source::Connection) -> &str;

/// Return one Python string per value, sharing one object among equal values.
fn shared_strings<'py, 'a, T>(
    py: Python<'py>,
    strings: &mut FastMap<&'a str, Bound<'py, PyString>>,
    values: impl Iterator<Item = &'a T>,
) -> Vec<Bound<'py, PyString>>
where
    T: ?Sized + PartialEq + AsRef<str> + 'a,
{
    // Most columns repeat the previous row's value, which skips the map.
    let mut last: Option<(&T, Bound<'py, PyString>)> = None;
    values
        .map(|value| {
            if let Some((previous, string)) = &last
                && **previous == *value
            {
                return string.clone();
            }
            let text = value.as_ref();
            let string = strings
                .entry(text)
                .or_insert_with(|| PyString::new(py, text))
                .clone();
            last = Some((value, string.clone()));
            string
        })
        .collect()
}

/// Concatenate fixed-width native-endian values into one byte buffer.
fn bytes<const N: usize>(count: usize, values: impl Iterator<Item = [u8; N]>) -> Vec<u8> {
    let mut output = Vec::with_capacity(count * N);
    for value in values {
        output.extend_from_slice(&value);
    }
    output
}
