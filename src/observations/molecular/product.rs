//! Native polymer chemistry and bounded ligand perception domains on one linked product.

use super::json;
use super::{Failure, Perception};
use crate::cif::{BlockCategory, CifDocument, CifValueRef, FastMap};
use crate::pdbx::{BondSource, ProductAtom, ProductBond};
use std::collections::HashMap;

/// AtomFeature columns written by native polymer chemistry and RDKit perception.
#[derive(Clone, Copy, serde::Deserialize)]
pub(crate) struct Columns {
    pub count: usize,
    pub atomic_number: usize,
    pub formal_charge: usize,
    pub total_hydrogen: usize,
    pub hybridization: usize,
    pub formal_charge_known: usize,
}

/// One residue occurrence keyed by chain instance, author number, insertion and component.
pub(super) struct Residues {
    pub of_atom: Vec<u32>,
    pub first: Vec<usize>,
    pub canonical: Vec<bool>,
}

/// Ligand residues, the canonical residues bonded to them and their atom union.
pub(super) struct Domain {
    pub atoms: Vec<usize>,
    pub hosts: Vec<u32>,
}

/// A perception copy whose cut bonds end in hydrogen caps.
pub(super) struct LocalGraph {
    pub ordered: Vec<usize>,
    pub caps: Vec<usize>,
    pub bonds: Vec<(usize, usize, u8)>,
    pub boundary: Vec<usize>,
    pub uncertain: Vec<usize>,
}

/// Map each struct_asym instance to its entity_poly type, as gemmi reads raw values.
pub(super) fn polymer_chains(document: &CifDocument) -> Result<HashMap<String, String>, Failure> {
    let Some(block) = document.blocks().first() else {
        return Ok(HashMap::new());
    };
    let column = |category: &str, item: &str| -> Vec<String> {
        let Some(table) = BlockCategory::find(block.entries(), category) else {
            return Vec::new();
        };
        let Some(column) = table.column(item) else {
            return Vec::new();
        };
        (0..table.row_count())
            .map(|row| match table.value(row, column) {
                Some(CifValueRef::Unknown) => "?".to_owned(),
                Some(CifValueRef::NotApplicable) | None => ".".to_owned(),
                Some(_) => table
                    .text(row, column)
                    .map_or_else(String::new, |text| text.into_owned()),
            })
            .collect()
    };
    let paired =
        |first: Vec<String>, second: Vec<String>| -> Result<Vec<(String, String)>, Failure> {
            if first.len() < second.len() {
                return Err(Failure::Value(
                    "zip() argument 2 is longer than argument 1".into(),
                ));
            }
            if first.len() > second.len() {
                return Err(Failure::Value(
                    "zip() argument 2 is shorter than argument 1".into(),
                ));
            }
            Ok(first.into_iter().zip(second).collect())
        };
    let polymers: HashMap<String, String> = paired(
        column("entity_poly", "entity_id"),
        column("entity_poly", "type"),
    )?
    .into_iter()
    .collect();
    Ok(paired(
        column("struct_asym", "id"),
        column("struct_asym", "entity_id"),
    )?
    .into_iter()
    .map(|(chain, entity)| (chain, polymers.get(&entity).cloned().unwrap_or_default()))
    .collect())
}

fn kind<'a>(chains: &'a HashMap<String, String>, asym: &str) -> &'a str {
    chains.get(asym).map_or("", String::as_str)
}

fn same_residue(first: &ProductAtom, second: &ProductAtom) -> bool {
    first.asym_id == second.asym_id
        && first.auth_sequence_id == second.auth_sequence_id
        && first.insertion_code == second.insertion_code
        && first.component_id == second.component_id
}

pub(super) fn residues(
    atoms: &[ProductAtom],
    chains: &HashMap<String, String>,
    perception: &Perception,
) -> Residues {
    let mut lookup: HashMap<(&str, &str, Option<&str>, &str), u32> = HashMap::new();
    let mut of_atom: Vec<u32> = Vec::with_capacity(atoms.len());
    let mut first = Vec::new();
    for (index, atom) in atoms.iter().enumerate() {
        if index > 0 && same_residue(&atoms[index - 1], atom) {
            of_atom.push(of_atom[index - 1]);
            continue;
        }
        let key = (
            atom.asym_id.as_str(),
            atom.auth_sequence_id.as_str(),
            atom.insertion_code.as_deref(),
            atom.component_id.as_str(),
        );
        let id = *lookup.entry(key).or_insert_with(|| {
            first.push(index);
            (first.len() - 1) as u32
        });
        of_atom.push(id);
    }
    let canonical = first
        .iter()
        .map(|&index| {
            let atom = &atoms[index];
            let kind = kind(chains, &atom.asym_id);
            (kind.contains("polypeptide")
                && perception.protein.contains(atom.component_id.as_str()))
                || (kind.contains("ribonucleotide")
                    && perception.nucleotides.contains(atom.component_id.as_str()))
        })
        .collect();
    Residues {
        of_atom,
        first,
        canonical,
    }
}

/// Python's repr of a residue key tuple, which orders seeds and hosts.
pub(super) fn residue_repr(atom: &ProductAtom) -> String {
    let mut out = String::from("(");
    json::repr(&mut out, &atom.asym_id);
    out.push_str(", ");
    json::repr(&mut out, &atom.auth_sequence_id);
    out.push_str(", ");
    match &atom.insertion_code {
        Some(code) => json::repr(&mut out, code),
        None => out.push_str("None"),
    }
    out.push_str(", ");
    json::repr(&mut out, &atom.component_id);
    out.push(')');
    out
}

/// The cached template hybridization of one component, by atom name.
type Template = HashMap<String, [i64; 4]>;

/// Copy native H/charge targets and unchanged cached template hybridization.
pub(super) fn native_features(
    atoms: &[ProductAtom],
    residues: &Residues,
    perception: &Perception,
    columns: Columns,
) -> Vec<[i16; 16]> {
    let mut features = vec![[0_i16; 16]; atoms.len()];
    // Few elements recur and a residue's atoms share one template, so both lookups are
    // kept: the numbers by element, the template until the residue changes.
    let mut numbers: FastMap<&str, i16> = FastMap::default();
    let mut template: Option<(u32, Option<&Template>)> = None;
    for (index, atom) in atoms.iter().enumerate() {
        let (Some(charge), Some(hydrogens)) = (atom.formal_charge, atom.hydrogen_count) else {
            continue;
        };
        let residue = residues.of_atom[index];
        if !residues.canonical[residue as usize] || !atom.active {
            continue;
        }
        let row = &mut features[index];
        row[columns.atomic_number] = *numbers.entry(atom.element.as_str()).or_insert_with(|| {
            perception
                .numbers
                .get(&atom.element.to_uppercase())
                .copied()
                .unwrap_or(0)
        });
        row[columns.formal_charge] = charge as i16;
        row[columns.formal_charge_known] = 1;
        row[columns.total_hydrogen] = hydrogens as i16;
        let name = atom
            .template_atom_id
            .as_deref()
            .filter(|name| !name.is_empty())
            .unwrap_or(&atom.atom_id);
        let names = match template {
            Some((last, names)) if last == residue => names,
            _ => {
                let names = perception.templates.get(atom.component_id.as_str());
                template = Some((residue, names));
                names
            }
        };
        if let Some(&[known, total, formal, hybridization]) =
            names.and_then(|names| names.get(name))
            && known != 0
            && total == hydrogens as i64
            && formal == i64::from(charge)
        {
            row[columns.hybridization] = hybridization as i16;
        }
    }
    features
}

fn polymer_link(
    first: &ProductAtom,
    second: &ProductAtom,
    bond: &ProductBond,
    chains: &HashMap<String, String>,
) -> bool {
    if bond.source == BondSource::PolymerSequence {
        return true;
    }
    let (a, b) = (kind(chains, &first.asym_id), kind(chains, &second.asym_id));
    let names = (
        first.atom_id.replace('*', "'"),
        second.atom_id.replace('*', "'"),
    );
    let pair = |x: &str, y: &str| (names.0 == x && names.1 == y) || (names.0 == y && names.1 == x);
    (a.contains("polypeptide") && b.contains("polypeptide") && pair("C", "N"))
        || (a.contains("ribonucleotide") && b.contains("ribonucleotide") && pair("O3'", "P"))
}

/// Group each ligand with its covalently linked ligands and direct canonical hosts.
pub(super) fn domains(
    atoms: &[ProductAtom],
    bonds: &[ProductBond],
    residues: &Residues,
    chains: &HashMap<String, String>,
    perception: &Perception,
) -> Vec<Domain> {
    let count = residues.first.len();
    let mut members: Vec<Vec<usize>> = vec![Vec::new(); count];
    for (index, &residue) in residues.of_atom.iter().enumerate() {
        members[residue as usize].push(index);
    }
    let seed: Vec<bool> = (0..count)
        .map(|residue| {
            let component = &atoms[residues.first[residue]].component_id;
            !residues.canonical[residue]
                && !perception.unknown.contains(component.as_str())
                && !perception.water.contains(component.as_str())
        })
        .collect();
    let mut neighbors: Vec<Vec<u32>> = vec![Vec::new(); count];
    for bond in bonds {
        let [a, b] = bond.atom_indices;
        if !bond.active || bond.order == 0 || !atoms[a].active || !atoms[b].active {
            continue;
        }
        let (first, second) = (residues.of_atom[a], residues.of_atom[b]);
        if first != second && !polymer_link(&atoms[a], &atoms[b], bond, chains) {
            neighbors[first as usize].push(second);
            neighbors[second as usize].push(first);
        }
    }
    for list in &mut neighbors {
        list.sort_unstable();
        list.dedup();
    }
    let mut seeds: Vec<(String, u32)> = (0..count as u32)
        .filter(|&residue| seed[residue as usize])
        .map(|residue| {
            (
                residue_repr(&atoms[residues.first[residue as usize]]),
                residue,
            )
        })
        .collect();
    seeds.sort();
    let mut visited = vec![false; count];
    let mut selected_mark = vec![u32::MAX; count];
    let mut host_mark = vec![u32::MAX; count];
    let mut result = Vec::new();
    for (stamp, &(_, start)) in seeds.iter().enumerate() {
        let stamp = stamp as u32;
        if visited[start as usize] {
            continue;
        }
        let (mut pending, mut selected, mut hosts) = (vec![start], Vec::new(), Vec::new());
        while let Some(current) = pending.pop() {
            if selected_mark[current as usize] == stamp {
                continue;
            }
            selected_mark[current as usize] = stamp;
            selected.push(current);
            for &neighbor in &neighbors[current as usize] {
                if seed[neighbor as usize] {
                    pending.push(neighbor);
                } else if residues.canonical[neighbor as usize] {
                    if host_mark[neighbor as usize] != stamp {
                        host_mark[neighbor as usize] = stamp;
                        hosts.push(neighbor);
                    }
                    pending.extend(
                        neighbors[neighbor as usize]
                            .iter()
                            .copied()
                            .filter(|&other| seed[other as usize]),
                    );
                }
            }
        }
        for &residue in &selected {
            visited[residue as usize] = true;
        }
        let mut union: Vec<usize> = selected
            .iter()
            .chain(&hosts)
            .flat_map(|&residue| members[residue as usize].iter().copied())
            .collect();
        union.sort_unstable();
        let mut hosts: Vec<(String, u32)> = hosts
            .into_iter()
            .map(|residue| {
                (
                    residue_repr(&atoms[residues.first[residue as usize]]),
                    residue,
                )
            })
            .collect();
        hosts.sort();
        result.push(Domain {
            atoms: union,
            hosts: hosts.into_iter().map(|(_, residue)| residue).collect(),
        });
    }
    result
}

/// Active, ordered bonds incident to each atom, in bond order.
pub(super) fn incident(atoms: &[ProductAtom], bonds: &[ProductBond]) -> (Vec<usize>, Vec<usize>) {
    let mut offsets = vec![0_usize; atoms.len() + 1];
    for bond in bonds.iter().filter(|bond| bond.active && bond.order != 0) {
        for &atom in &bond.atom_indices {
            offsets[atom + 1] += 1;
        }
    }
    for index in 0..atoms.len() {
        offsets[index + 1] += offsets[index];
    }
    let mut cursor = offsets.clone();
    let mut edges = vec![0_usize; offsets[atoms.len()]];
    for (edge, bond) in bonds
        .iter()
        .enumerate()
        .filter(|(_, bond)| bond.active && bond.order != 0)
    {
        for &atom in &bond.atom_indices {
            edges[cursor[atom]] = edge;
            cursor[atom] += 1;
        }
    }
    (offsets, edges)
}

/// Saturate cut bonds only in a perception copy and mark boundary labels unknown.
pub(super) fn local_graph(
    atoms: &[ProductAtom],
    bonds: &[ProductBond],
    selected: &[usize],
    incident: &(Vec<usize>, Vec<usize>),
    local: &mut [u32],
) -> LocalGraph {
    const OUTSIDE: u32 = u32::MAX;
    for (new, &old) in selected.iter().enumerate() {
        local[old] = new as u32;
    }
    let (offsets, incident) = incident;
    let mut edges: Vec<usize> = selected
        .iter()
        .flat_map(|&atom| incident[offsets[atom]..offsets[atom + 1]].iter().copied())
        .collect();
    edges.sort_unstable();
    edges.dedup();
    let mut graph = LocalGraph {
        ordered: selected.to_vec(),
        caps: Vec::new(),
        bonds: Vec::new(),
        boundary: Vec::new(),
        uncertain: Vec::new(),
    };
    for &edge in &edges {
        let bond = &bonds[edge];
        let [first, second] = bond.atom_indices;
        let (inner_first, inner_second) = (local[first] != OUTSIDE, local[second] != OUTSIDE);
        if inner_first && inner_second {
            graph
                .bonds
                .push((local[first] as usize, local[second] as usize, bond.order));
        } else if inner_first != inner_second {
            let (old, other) = if inner_first {
                (first, second)
            } else {
                (second, first)
            };
            if !atoms[old].active || !atoms[other].active {
                continue;
            }
            graph.boundary.push(old);
            for _ in 0..bond.order {
                let cap = selected.len() + graph.caps.len();
                graph.caps.push(old);
                graph.bonds.push((local[old] as usize, cap, 1));
            }
        }
    }
    graph.boundary.sort_unstable();
    graph.boundary.dedup();
    let boundary = |atom: usize| graph.boundary.binary_search(&atom).is_ok();
    let mut uncertain = graph.boundary.clone();
    for &edge in &edges {
        let [first, second] = bonds[edge].atom_indices;
        if boundary(first) && local[second] != OUTSIDE {
            uncertain.push(second);
        }
        if boundary(second) && local[first] != OUTSIDE {
            uncertain.push(first);
        }
    }
    uncertain.sort_unstable();
    uncertain.dedup();
    graph.uncertain = uncertain;
    for &old in selected {
        local[old] = OUTSIDE;
    }
    graph
}

/// Residues with complete active heavy chemistry, and residues joined by a non-peptide link.
pub(super) fn completeness(
    atoms: &[ProductAtom],
    bonds: &[ProductBond],
    residues: &Residues,
    features: &[[i16; 16]],
    columns: Columns,
) -> (Vec<bool>, Vec<bool>) {
    let count = residues.first.len();
    let (mut present, mut incomplete, mut modified) =
        (vec![false; count], vec![false; count], vec![false; count]);
    for (index, atom) in atoms.iter().enumerate() {
        if !atom.active || matches!(atom.element.as_str(), "H" | "D" | "T") {
            continue;
        }
        let residue = residues.of_atom[index] as usize;
        present[residue] = true;
        if atom.coordinates.is_none()
            || atom.hydrogen_count.is_none()
            || atom.formal_charge.is_none()
            || features[index][columns.atomic_number] == 0
        {
            incomplete[residue] = true;
        }
    }
    for bond in bonds {
        if !bond.active || bond.order == 0 || bond.source == BondSource::PolymerSequence {
            continue;
        }
        let [a, b] = bond.atom_indices;
        let (first, second) = (residues.of_atom[a] as usize, residues.of_atom[b] as usize);
        if first != second {
            let (x, y) = (atoms[a].atom_id.as_str(), atoms[b].atom_id.as_str());
            if !((x == "C" && y == "N") || (x == "N" && y == "C")) {
                modified[first] = true;
                modified[second] = true;
            }
        }
    }
    let complete = present
        .iter()
        .zip(&incomplete)
        .map(|(&present, &incomplete)| present && !incomplete)
        .collect();
    (complete, modified)
}
