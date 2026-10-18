//! Selected source atoms with product chemistry, SASA targets and their known masks.

use super::arrays::{self, Array};
use super::json;
use super::product::{Columns, Residues};
use crate::cif::CifDocument;
use crate::observations::source::{Connection, Source};
use crate::pdbx::ShortText;
use crate::observations::uncertainty::TargetUncertainty;
use crate::pdbx::{ProductAtom, ProductBond};
use std::collections::{HashMap, HashSet};

/// Numeric source fields in the order Python's source snapshot stores them.
pub(super) struct Snapshot {
    pub row: Vec<i64>,
    pub group: Vec<i64>,
    pub sequence: Vec<i64>,
    pub hydrogen_parent: Vec<i64>,
    pub pos: Vec<[f32; 3]>,
    pub occupancy: Vec<f32>,
    pub mask: Vec<bool>,
    pub is_protein: Vec<bool>,
    pub connections: Vec<[i64; 2]>,
    pub connection_row: Vec<i64>,
    pub charge: Vec<f32>,
    pub charge_known: Vec<bool>,
    pub hybridization: Vec<i64>,
    pub hybridization_known: Vec<bool>,
    pub features: Vec<[i16; 16]>,
    pub connectivity_known: Vec<bool>,
    pub hydrogen_count: Vec<i64>,
    pub hydrogen_known: Vec<bool>,
    pub provenance: String,
    pub bonds: Vec<[i64; 2]>,
    pub bond_order: Vec<i8>,
    pub bond_order_known: Vec<bool>,
    pub sasa: Vec<f32>,
    pub sasa_known: Vec<bool>,
    pub sasa_provenance: String,
    pub target_known: Vec<bool>,
    /// Upper-cased element is H, D or T.
    pub hydrogen: Vec<bool>,
    /// The component has no dictionary chemistry.
    unknown: Vec<bool>,
}

/// Whether the upper-cased element is H, D or T. Only these six characters upper-case to one of
/// them, and no longer text upper-cases to one character. A test against `&str` reads the
/// bytes and skips the UTF-8 check of `as_str`.
fn is_hydrogen(element: &ShortText) -> bool {
    ["H", "D", "T", "h", "d", "t"].iter().any(|value| element == value)
}

/// Reads one text field of a source connection.
type EdgeField = fn(&Connection) -> &ShortText;

/// Longest escaped label, with its leading comma, that a repeat copies as a fixed-size block.
const PIECE: usize = 32;

/// One label list as JSON strings joined by commas, without its brackets.
struct Column<'a, T: ?Sized> {
    out: Vec<u8>,
    previous: Option<&'a T>,
    piece: [u8; PIECE],
    /// Bytes of `piece` in use, or zero when the escaped value does not fit.
    length: usize,
}

impl<'a, T: PartialEq + ?Sized> Column<'a, T> {
    fn new(capacity: usize) -> Self {
        Self { out: Vec::with_capacity(capacity), previous: None, piece: [0; PIECE], length: 0 }
    }

    /// Append one value. Labels often repeat from one atom to the next. A value equal to
    /// the one before copies its escaped text from a fixed-size block, which needs no
    /// `memcpy` call.
    fn push(&mut self, value: &'a T, scratch: &mut String, text: impl FnOnce(&'a T) -> &'a str) {
        if self.length > 0 && self.previous.is_some_and(|last| last == value) {
            self.out.extend_from_slice(&self.piece);
            let end = self.out.len() - PIECE + self.length;
            self.out.truncate(end);
            return;
        }
        let text = text(value);
        let bytes = text.as_bytes();
        // Short text that needs no escapes goes into the block with fixed-size copies.
        if self.previous.is_some() && bytes.len() <= PIECE - 3
            && bytes.iter().all(|&byte| matches!(byte, b' '..=b'~') && byte != b'"' && byte != b'\\') {
            self.piece[0] = b',';
            self.piece[1] = b'"';
            copy_short(&mut self.piece, bytes);
            self.piece[bytes.len() + 2] = b'"';
            self.length = bytes.len() + 3;
            self.out.extend_from_slice(&self.piece);
            let end = self.out.len() - PIECE + self.length;
            self.out.truncate(end);
            self.previous = Some(value);
            return;
        }
        scratch.clear();
        scratch.push(',');
        json::string(scratch, text);
        let bytes = scratch.as_bytes();
        self.out.extend_from_slice(if self.previous.is_some() { bytes } else { &bytes[1..] });
        self.length = if bytes.len() <= PIECE {
            self.piece[..bytes.len()].copy_from_slice(bytes);
            bytes.len()
        } else {
            0
        };
        self.previous = Some(value);
    }
}

/// Copy `text`, of at most `PIECE - 3` bytes, to `piece` from index 2.
///
/// Two fixed-size copies that may overlap compile to plain loads and stores, where a copy
/// of `text.len()` bytes would call `memcpy`.
fn copy_short(piece: &mut [u8; PIECE], text: &[u8]) {
    fn fixed<const N: usize>(piece: &mut [u8; PIECE], text: &[u8], start: usize) {
        let mut block = [0; N];
        block.copy_from_slice(&text[start..start + N]);
        piece[2 + start..2 + start + N].copy_from_slice(&block);
    }
    let length = text.len();
    match length {
        0 => {}
        1..=3 => {
            piece[2] = text[0];
            piece[2 + length / 2] = text[length / 2];
            piece[1 + length] = text[length - 1];
        }
        4..=7 => {
            fixed::<4>(piece, text, 0);
            fixed::<4>(piece, text, length - 4);
        }
        8..=16 => {
            fixed::<8>(piece, text, 0);
            fixed::<8>(piece, text, length - 8);
        }
        _ => {
            fixed::<16>(piece, text, 0);
            fixed::<16>(piece, text, length - 16);
        }
    }
}

/// Append `"name":[values]` to a JSON object, with a comma unless it is the first member.
/// One label list as UTF-8 values zero-padded to eight bytes, which NumPy reads as `S8`.
struct Fixed {
    values: Vec<[u8; 8]>,
    /// Whether every value fits and keeps its last byte under NumPy, which drops
    /// trailing NUL bytes.
    exact: bool,
}

impl Fixed {
    fn new(count: usize) -> Self {
        Self { values: Vec::with_capacity(count), exact: true }
    }

    fn push(&mut self, text: &ShortText) {
        let bytes = text.as_bytes();
        self.exact &= bytes.len() <= 8 && bytes.last() != Some(&0);
        // Folding from the last byte puts the first byte lowest, so the little-endian
        // bytes of the value are the text followed by zeros.
        let value = bytes.iter().rev().fold(0_u64, |value, &byte| value << 8 | u64::from(byte));
        self.values.push(value.to_le_bytes());
    }

    /// The values as an `S8` array, or `None` when one does not survive the padding.
    fn array(self) -> Option<Array> {
        let count = self.values.len();
        self.exact.then(|| Array { dtype: "S8", shape: vec![count], bytes: self.values.into_flattened() })
    }
}

fn member<T: ?Sized>(out: &mut Vec<u8>, scratch: &mut String, name: &str, column: &Column<'_, T>) {
    scratch.clear();
    if out.len() > 1 {
        scratch.push(',');
    }
    json::string(scratch, name);
    scratch.push_str(":[");
    out.extend_from_slice(scratch.as_bytes());
    out.extend_from_slice(&column.out);
    out.push(b']');
}

fn normal(value: &str) -> &str {
    if matches!(value, "." | "?") { "" } else { value }
}

impl Snapshot {
    /// `unknown` names the components without dictionary chemistry.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn read(source: &Source, atoms: &[ProductAtom], bonds: &[ProductBond], rows: &[Option<usize>],
        features: &[[i16; 16]], residues: &Residues, complete: &[bool], columns: Columns,
        provenance: String, unknown: &HashSet<String>) -> Self {
        let count = source.atoms.len();
        let mut hydrogen_parent = vec![-1_i64; count];
        for &(child, parent) in &source.hydrogen_bonds {
            if let Some(parent) = parent {
                hydrogen_parent[child] = if hydrogen_parent[child] == -1 || hydrogen_parent[child] == parent as i64 {
                    parent as i64
                } else {
                    -2
                };
            }
        }
        for value in &mut hydrogen_parent {
            if *value < 0 {
                *value = -1;
            }
        }
        // The atom array is much larger than the cache, so read it once for all fields.
        let mut hydrogen = Vec::with_capacity(count);
        let mut row = Vec::with_capacity(count);
        let mut group = Vec::with_capacity(count);
        let mut sequence = Vec::with_capacity(count);
        let mut pos = Vec::with_capacity(count);
        let mut occupancy = Vec::with_capacity(count);
        let mut mask = Vec::with_capacity(count);
        let mut is_protein = Vec::with_capacity(count);
        let mut unknown_flags = Vec::with_capacity(count);
        // Atoms of one residue share a component, so test each run of one component once.
        let mut previous: Option<(&ShortText, bool)> = None;
        for atom in &source.atoms {
            let flag = match previous {
                Some((component, flag)) if *component == atom.component => flag,
                _ => {
                    let flag = unknown.contains(atom.component.as_str());
                    previous = Some((&atom.component, flag));
                    flag
                }
            };
            unknown_flags.push(flag);
            hydrogen.push(is_hydrogen(&atom.element));
            row.push(i64::from(atom.source_row));
            group.push(atom.group as i64);
            sequence.push(atom.sequence);
            pos.push(atom.position);
            occupancy.push(atom.occupancy);
            mask.push(atom.position.iter().all(|value| value.is_finite()));
            is_protein.push(atom.protein);
        }
        let limit = row.iter().map(|&value| value as usize + 1).max().unwrap_or(0);
        let mut by_row = vec![usize::MAX; limit];
        for (index, &value) in row.iter().enumerate() {
            by_row[value as usize] = index;
        }
        let mut result = Self {
            row,
            group,
            sequence,
            hydrogen_parent,
            pos,
            occupancy,
            mask,
            is_protein,
            connections: source.connections.iter()
                .map(|edge| edge.atoms.map(|atom| atom.map_or(-1, |row| row as i64))).collect(),
            connection_row: source.connections.iter().map(|edge| edge.row as i64).collect(),
            charge: vec![0.; count],
            charge_known: vec![false; count],
            hybridization: vec![0; count],
            hybridization_known: vec![false; count],
            features: vec![[0; 16]; count],
            connectivity_known: vec![false; count],
            hydrogen_count: vec![0; count],
            hydrogen_known: vec![false; count],
            provenance,
            bonds: Vec::new(),
            bond_order: Vec::new(),
            bond_order_known: Vec::new(),
            sasa: vec![0.; count],
            sasa_known: vec![false; count],
            sasa_provenance: String::new(),
            target_known: Vec::new(),
            hydrogen,
            unknown: unknown_flags,
        };
        let mut mapping = vec![usize::MAX; atoms.len()];
        for (index, atom) in atoms.iter().enumerate() {
            let Some(row) = rows[index].and_then(|row| by_row.get(row).copied()).filter(|&row| row != usize::MAX)
            else {
                continue;
            };
            if !atom.active {
                continue;
            }
            mapping[index] = row;
            let values = features[index];
            result.features[row] = values;
            let known = complete[residues.of_atom[index] as usize];
            result.connectivity_known[row] = known;
            let hydrogens = i64::from(values[columns.total_hydrogen]);
            result.hydrogen_count[row] = hydrogens;
            result.hydrogen_known[row] = known && (0..=3).contains(&hydrogens);
            result.charge[row] = f32::from(values[columns.formal_charge]);
            result.charge_known[row] = known && values[columns.formal_charge_known] != 0;
            result.hybridization[row] = i64::from(values[columns.hybridization]);
            result.hybridization_known[row] = known && values[columns.hybridization] != 0;
        }
        let mut parents: Vec<(usize, bool)> = vec![(usize::MAX, false); count];
        for bond in bonds {
            let [first, second] = bond.atom_indices.map(|atom| mapping[atom]);
            if !bond.active || first == usize::MAX || second == usize::MAX {
                continue;
            }
            result.bonds.push([first as i64, second as i64]);
            result.bond_order.push(bond.order as i8);
            if bond.order == 0 {
                result.connectivity_known[first] = false;
                result.connectivity_known[second] = false;
            }
            let (a, b) = (result.hydrogen[first], result.hydrogen[second]);
            if a != b && bond.order != 0 {
                let (child, parent) = if a { (first, second) } else { (second, first) };
                let entry = &mut parents[child];
                if entry.0 == usize::MAX {
                    entry.0 = parent;
                } else if entry.0 != parent {
                    entry.1 = true;
                }
            }
        }
        for (child, &(parent, several)) in parents.iter().enumerate() {
            if parent != usize::MAX {
                result.hydrogen_parent[child] = if several { -1 } else { parent as i64 };
            }
        }
        result.bond_order_known = result.bond_order.iter().map(|&order| order > 0).collect();
        result
    }

    /// Shrake-Rupley SASA of noncanonical protein and nonprotein heavy atoms.
    pub(super) fn sasa(&mut self, source: &Source, document: &CifDocument, model: i64, enabled: bool,
        canonical: &HashSet<String>) {
        // Atoms of one residue share a component, so test each run of one component once.
        let mut previous: Option<(&ShortText, bool, bool)> = None;
        let selected: Vec<usize> = source.atoms.iter().enumerate().filter(|(index, atom)| {
            let (water, standard) = match previous {
                Some((component, water, standard)) if *component == atom.component => (water, standard),
                _ => {
                    let water = matches!(atom.component.to_uppercase().as_str(), "HOH" | "DOD" | "WAT");
                    let standard = canonical.contains(atom.component.as_str());
                    previous = Some((&atom.component, water, standard));
                    (water, standard)
                }
            };
            !self.hydrogen[*index]
                && !water
                && (!atom.protein || !standard)
                && self.mask[*index]
                && self.occupancy[*index] > 0.
        }).map(|(index, _)| index).collect();
        let mut error = None;
        if enabled && !selected.is_empty() {
            let rows: Vec<usize> = selected.iter().map(|&index| self.row[index] as usize).collect();
            match crate::pdbx::ligand_sasa(document, Some(&rows), model, 1.4, 256, &HashMap::new()) {
                Ok((found, areas)) => {
                    let areas: HashMap<usize, f64> = found.into_iter().zip(areas).collect();
                    for (&index, row) in selected.iter().zip(&rows) {
                        self.sasa[index] = areas.get(row).copied().unwrap_or(f64::NAN) as f32;
                        self.sasa_known[index] = true;
                    }
                }
                Err(message) => error = Some(message),
            }
        }
        let known = self.sasa_known.iter().filter(|&&value| value).count();
        let mut text = String::from("{\"context\":\"all-occupied-nonwater-heavy-in-selected-model\",\"enabled\":");
        text.push_str(if enabled { "true" } else { "false" });
        text.push_str(",\"error\":");
        match &error {
            Some(message) => json::string(&mut text, message),
            None => text.push_str("null"),
        }
        text.push_str(&format!(",\"hydrogens_occlude\":false,\"known_count\":{known},\
            \"method\":\"Shrake-Rupley-Fibonacci\",\"model\":{model},\"point_count\":256,\"probe_radius\":1.4,\
            \"radii\":\"Nibbler-Bondi\",\"radii_overrides\":{{}},\"schema\":\"nibbler-source-atoms.v1\",\
            \"selection\":\"noncanonical-protein-and-nonprotein-heavy\",\"target_count\":{},\
            \"units\":\"angstrom^2\"}}", selected.len()));
        self.sasa_provenance = text;
    }

    /// Keep raw presence while withholding uncertain added-hydrogen targets.
    pub(super) fn mask_uncertain_hydrogens(&mut self, source: &Source, uncertainties: &[TargetUncertainty]) {
        self.target_known = self.mask.iter().zip(&self.occupancy).map(|(&mask, &occupancy)| mask && occupancy > 0.)
            .collect();
        let missing: HashSet<[&str; 5]> = uncertainties.iter().map(|row| row.parent_key.each_ref().map(|value| normal(value)))
            .collect();
        if missing.is_empty() {
            return;
        }
        let names: HashSet<&str> = missing.iter().map(|key| key[4]).collect();
        let parent: Vec<bool> = source.atoms.iter().map(|atom| {
            names.contains(normal(&atom.name)) && missing.contains(&[normal(&atom.chain), normal(&atom.auth_sequence),
                normal(&atom.insertion), normal(&atom.component), normal(&atom.name)])
        }).collect();
        let groups: HashSet<i64> = self.group.iter().zip(&parent).filter(|(_, parent)| **parent)
            .map(|(&group, _)| group).collect();
        for atom in 0..parent.len() {
            let owner = self.hydrogen_parent[atom];
            let child = if owner >= 0 {
                parent[owner as usize]
            } else {
                self.hydrogen[atom] && groups.contains(&self.group[atom])
            };
            if child {
                self.target_known[atom] = false;
            }
            if parent[atom] {
                self.hydrogen_known[atom] = false;
                self.charge_known[atom] = false;
            }
        }
    }

    /// Keep unknown-component geometry without attesting dictionary chemistry.
    pub(super) fn mask_unknown_chemistry(&mut self) {
        let flags = std::mem::take(&mut self.unknown);
        if !flags.contains(&true) {
            return;
        }
        for (atom, _) in flags.iter().enumerate().filter(|(_, flag)| **flag) {
            self.charge_known[atom] = false;
            self.hydrogen_known[atom] = false;
            self.hybridization_known[atom] = false;
            self.connectivity_known[atom] = false;
            self.hybridization[atom] = 0;
            if self.hydrogen[atom] {
                self.target_known[atom] = false;
            }
        }
        for (edge, pair) in self.bonds.iter().enumerate() {
            if pair.iter().any(|&atom| flags[atom.max(0) as usize]) {
                self.bond_order_known[edge] = false;
            }
        }
    }

    /// Numeric fields, then labels and exclusion codes as UTF-8 JSON. With `fixed`, the
    /// component, element and name labels also come back as `S8` arrays when they fit.
    pub(super) fn arrays(&self, source: &Source, exclusions: &[String], fixed: bool) -> Vec<(&'static str, Array)> {
        let count = self.row.len();
        // The atom array is much larger than the cache, so read it once for every label column.
        let mut scratch = String::new();
        let mut columns: [Column<'_, ShortText>; 10] = std::array::from_fn(|_| Column::new(count * 6));
        let mut padded: [Fixed; 3] = std::array::from_fn(|_| Fixed::new(if fixed { source.atoms.len() } else { 0 }));
        for atom in &source.atoms {
            let values = [&atom.alternate, &atom.auth_chain, &atom.auth_sequence, &atom.chain, &atom.component,
                &atom.element, &atom.entity, &atom.insertion, &atom.name];
            for (column, value) in columns.iter_mut().zip(values) {
                column.push(value, &mut scratch, ShortText::as_str);
            }
            columns[9].push(&atom.entity, &mut scratch,
                |entity| source.entity_types.get(entity.as_str()).map_or("", String::as_str));
            if fixed {
                for (label, value) in padded.iter_mut().zip([&atom.component, &atom.element, &atom.name]) {
                    label.push(value);
                }
            }
        }
        let mut labels = Vec::with_capacity(columns.iter().map(|column| column.out.len() + 40).sum::<usize>()
            + source.connections.len() * 64);
        labels.push(b'{');
        let names = ["source_atom_alternate", "source_atom_auth_chain", "source_atom_auth_sequence",
            "source_atom_chain", "source_atom_component", "source_atom_element", "source_atom_entity",
            "source_atom_insertion", "source_atom_name", "source_atom_polymer_type"];
        for (name, column) in names.into_iter().zip(&columns) {
            member(&mut labels, &mut scratch, name, column);
        }
        let edges = &source.connections;
        let edge_fields: [(&str, EdgeField); 4] = [
            ("source_connection_order", |edge| &edge.order),
            ("source_connection_symmetry1", |edge| &edge.symmetries[0]),
            ("source_connection_symmetry2", |edge| &edge.symmetries[1]),
            ("source_connection_type", |edge| &edge.kind),
        ];
        for (name, field) in edge_fields {
            let mut column = Column::new(edges.len() * 8);
            for edge in edges {
                column.push(field(edge), &mut scratch, ShortText::as_str);
            }
            member(&mut labels, &mut scratch, name, &column);
        }
        labels.push(b'}');
        let mut codes = String::from("{");
        for (index, name) in exclusions.iter().enumerate() {
            if index > 0 {
                codes.push(',');
            }
            codes.push_str(&format!("\"{index}\":"));
            json::string(&mut codes, name);
        }
        codes.push('}');
        let features: Vec<i16> = self.features.iter().flatten().copied().collect();
        let fixed = if fixed {
            let [component, element, name] = padded.map(Fixed::array);
            component.zip(element).zip(name)
        } else {
            None
        };
        let mut out = vec![
            ("source_atom_row", arrays::int64(&self.row, vec![count])),
            ("source_atom_group", arrays::int64(&self.group, vec![count])),
            ("source_atom_sequence", arrays::int64(&self.sequence, vec![count])),
            ("source_hydrogen_parent", arrays::int64(&self.hydrogen_parent, vec![count])),
            ("source_atom_pos", arrays::points(&self.pos, vec![count, 3])),
            ("source_atom_occupancy", arrays::float32(&self.occupancy, vec![count])),
            ("source_atom_mask", arrays::boolean(&self.mask, vec![count])),
            ("source_atom_is_protein", arrays::boolean(&self.is_protein, vec![count])),
            ("source_connections", arrays::pairs(&self.connections)),
            ("source_connection_row", arrays::int64(&self.connection_row, vec![self.connection_row.len()])),
            ("source_atom_charge", arrays::float32(&self.charge, vec![count])),
            ("source_atom_charge_known", arrays::boolean(&self.charge_known, vec![count])),
            ("source_atom_hybridization", arrays::int64(&self.hybridization, vec![count])),
            ("source_atom_hybridization_known", arrays::boolean(&self.hybridization_known, vec![count])),
            ("source_atom_features", arrays::int16(&features, vec![count, 16])),
            ("source_atom_connectivity_known", arrays::boolean(&self.connectivity_known, vec![count])),
            ("source_atom_hydrogen_count", arrays::int64(&self.hydrogen_count, vec![count])),
            ("source_atom_hydrogen_known", arrays::boolean(&self.hydrogen_known, vec![count])),
            ("source_chemistry_provenance", arrays::text(self.provenance.clone())),
            ("source_bonds", arrays::pairs(&self.bonds)),
            ("source_bond_order", arrays::int8(&self.bond_order, vec![self.bond_order.len()])),
            ("source_bond_order_known", arrays::boolean(&self.bond_order_known, vec![self.bond_order.len()])),
            ("source_sasa", arrays::float32(&self.sasa, vec![count])),
            ("source_sasa_mask", arrays::boolean(&self.sasa_known, vec![count])),
            ("source_sasa_provenance", arrays::text(self.sasa_provenance.clone())),
            ("source_atom_target_known", arrays::boolean(&self.target_known, vec![count])),
            ("source_labels", { let count = labels.len(); arrays::uint8(labels, vec![count]) }),
            ("source_exclusion_codes", arrays::text(codes)),
        ];
        if let Some(((component, element), name)) = fixed {
            out.extend([("source_fixed_component", component), ("source_fixed_element", element),
                ("source_fixed_name", name)]);
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fixed_labels_pad_to_eight_bytes_or_refuse() {
        let mut labels = Fixed::new(3);
        for text in ["CA", "", "A\0B"] {
            labels.push(&ShortText::new(text));
        }
        assert!(labels.array().is_some_and(|array| array.dtype == "S8" && array.shape == [3]
            && array.bytes == b"CA\0\0\0\0\0\0\0\0\0\0\0\0\0\0A\0B\0\0\0\0\0"));
        for text in ["ABCDEFGHI", "A\0"] {
            let mut labels = Fixed::new(1);
            labels.push(&ShortText::new(text));
            assert!(labels.array().is_none());
        }
    }
}
