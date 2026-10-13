//! Whole-unit overlap checks for copies of the same non-polymer source.

use super::SemanticError;
use serde::Serialize;
use std::collections::{BTreeMap, HashMap, HashSet};

/// Source-row selection and the map for removed assembly atoms.
#[derive(Clone, Debug, Serialize)]
pub struct AssemblyDeduplication {
    /// Input row indices retained in source order.
    pub retained_rows: Vec<usize>,
    /// Removed input row and its retained same-source-atom row.
    pub duplicate_pairs: Vec<[usize; 2]>,
}

/// Remove duplicate non-polymer copies with at least 50% overlap within 1 A.
///
/// Unit identifiers must name the same source chemical unit across copies.
/// Source atom identifiers map removed atoms to the retained copy. Lower copy
/// identifiers take precedence; callers must rank identity operators first.
/// Only eligible units participate. Mark covalently attached or unresolved
/// multi-atom coordination units ineligible unless their atom mapping is exact.
///
/// # Errors
///
/// Rejects mismatched arrays, invalid coordinates and unequal source atom sets.
pub fn deduplicate_assembly_copies(
    positions: &[[f64; 3]],
    source_atoms: &[usize],
    units: &[usize],
    copies: &[usize],
    eligible: &[bool],
    occupied: &[bool],
) -> Result<AssemblyDeduplication, SemanticError> {
    let n = positions.len();
    if [
        source_atoms.len(),
        units.len(),
        copies.len(),
        eligible.len(),
        occupied.len(),
    ]
    .iter()
    .any(|&m| m != n)
        || positions
            .iter()
            .flatten()
            .any(|x| !x.is_finite() || x.abs() > 1e12)
    {
        return Err(error("Assembly overlap arrays or coordinates are invalid"));
    }
    let mut groups: BTreeMap<usize, BTreeMap<usize, Vec<usize>>> = BTreeMap::new();
    for i in 0..n {
        if eligible[i] {
            groups
                .entry(units[i])
                .or_default()
                .entry(copies[i])
                .or_default()
                .push(i);
        }
    }
    let mut removed = HashSet::new();
    let mut pairs = Vec::new();
    for variants in groups.values() {
        if variants.len() < 2 {
            continue;
        }
        let mut cells: HashMap<[i64; 3], Vec<(usize, usize)>> = HashMap::new();
        let mut retained: BTreeMap<usize, BTreeMap<usize, usize>> = BTreeMap::new();
        for (&copy, rows) in variants {
            let source: BTreeMap<_, _> = rows.iter().map(|&i| (source_atoms[i], i)).collect();
            if source.len() != rows.len() {
                return Err(error("Duplicate source atom within one assembly copy"));
            }
            if retained
                .values()
                .any(|other| other.keys().ne(source.keys()))
            {
                return Err(error("Assembly copies have different source atom sets"));
            }
            let mut counts: BTreeMap<usize, usize> = BTreeMap::new();
            let observed: Vec<_> = rows.iter().copied().filter(|&i| occupied[i]).collect();
            for &i in &observed {
                let cell = positions[i].map(|x| x.floor() as i64);
                let mut seen = HashSet::new();
                for dx in -1..=1 {
                    for dy in -1..=1 {
                        for dz in -1..=1 {
                            if let Some(nearby) =
                                cells.get(&[cell[0] + dx, cell[1] + dy, cell[2] + dz])
                            {
                                for &(earlier, j) in nearby {
                                    let d: f64 = positions[i]
                                        .iter()
                                        .zip(positions[j])
                                        .map(|(x, y)| (x - y).powi(2))
                                        .sum();
                                    if d < 1. {
                                        seen.insert(earlier);
                                    }
                                }
                            }
                        }
                    }
                }
                for earlier in seen {
                    *counts.entry(earlier).or_default() += 1;
                }
            }
            let duplicate = counts.iter().find_map(|(&earlier, &count)| {
                (count.saturating_mul(2) >= observed.len() && !observed.is_empty())
                    .then_some(earlier)
            });
            if let Some(target) = duplicate.and_then(|key| retained.get(&key)) {
                for (&atom, &i) in &source {
                    if let Some(&j) = target.get(&atom) {
                        removed.insert(i);
                        pairs.push([i, j]);
                    }
                }
            } else {
                for &i in &observed {
                    cells
                        .entry(positions[i].map(|x| x.floor() as i64))
                        .or_default()
                        .push((copy, i));
                }
                retained.insert(copy, source);
            }
        }
    }
    pairs.sort_unstable();
    Ok(AssemblyDeduplication {
        retained_rows: (0..n).filter(|i| !removed.contains(i)).collect(),
        duplicate_pairs: pairs,
    })
}

fn error(message: &str) -> SemanticError {
    SemanticError::new("PDBX_ASSEMBLY_OVERLAP", message, Vec::new())
}
