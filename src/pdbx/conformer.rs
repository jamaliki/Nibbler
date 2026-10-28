//! Clash-aware conformer selection with bounded search and sparse contacts.

use super::SemanticError;

mod contacts;
mod source;

use contacts::contacts;
use serde::Serialize;
pub use source::{
    select_conformer, select_conformer_document, select_conformer_document_with_policy,
    select_conformer_with_policy,
};
use std::collections::{BTreeMap, HashMap, HashSet};

/// One input site. Group identifiers name residue occurrences.
#[derive(Clone, Debug)]
pub struct ConformerAtom {
    /// Residue occurrence index.
    pub residue: usize,
    /// Source chain instance.
    pub chain: String,
    /// Deposited sequence index, when known.
    pub sequence: Option<i64>,
    /// Atom name used for backbone bond exclusions.
    pub atom_name: String,
    /// Alternate identifier, or None for shared atoms.
    pub alternate: Option<String>,
    /// Supplied Cartesian coordinates in angstroms.
    pub coordinates: [f64; 3],
    /// Deposited occupancy.
    pub occupancy: f64,
    /// Whether this site contributes to heavy-atom contact constraints.
    pub clash_eligible: bool,
}

/// Selected source rows and reproducible policy metadata.
#[derive(Clone, Debug, Serialize)]
pub struct ConformerSelection {
    /// Version of the selection rules.
    pub schema_version: &'static str,
    /// Indices into the supplied atom list, in source order.
    pub atom_rows: Vec<usize>,
    /// Chosen alternate per variable residue occurrence.
    pub assignments: BTreeMap<usize, String>,
    /// Search choices examined, including backtracking.
    pub search_nodes: usize,
    /// Seed used only to order equal-scoring choices.
    pub seed: u64,
    /// Selected deposited model number; None for the array API.
    pub model_number: Option<i64>,
}

#[derive(Default)]
struct Group {
    residue: usize,
    labels: Vec<String>,
    occupancy: Vec<f64>,
    blocked: HashSet<usize>,
    neighbors: HashSet<usize>,
}

type Pair = (usize, usize);
type Constraints = HashMap<Pair, HashSet<Pair>>;

/// Search ordering does not change source occupancies or clash constraints.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ConformerPolicy {
    /// Prefer higher summed residue occupancy; use the seed for ties.
    Occupancy,
    /// Seeded random alternative order, excluding zero-occupancy alternatives.
    RandomClashAware,
}

impl ConformerPolicy {
    pub(super) fn schema(self) -> &'static str {
        match self {
            Self::Occupancy => "clash-aware-altloc.v1",
            Self::RandomClashAware => "random-clash-aware-altloc.v1",
        }
    }
}

fn fail(message: impl Into<String>) -> SemanticError {
    SemanticError::new("PDBX_CONFORMER_CONFLICT", message, vec!["atom_site".into()])
}

/// Select residue-level alternatives while avoiding nonbonded clashes below 1.5 A.
///
/// Bond exclusions name source-row pairs. Proximal residues on different chains
/// prefer a shared alternate letter when their occupancy scores tie. Shared atoms
/// remain selected. Search tries higher occupancy first and fails after 100,000
/// choices; it never renames a component or returns a known clashing assignment.
///
/// # Errors
///
/// Rejects nonfinite coordinates, invalid exclusions, and unsatisfiable or overly
/// large correlated components. No file or network access occurs.
pub fn select_atom_conformers(
    atoms: &[ConformerAtom],
    excluded_bonds: &[[usize; 2]],
    seed: u64,
) -> Result<ConformerSelection, SemanticError> {
    select_atom_conformers_with_policy(atoms, excluded_bonds, seed, ConformerPolicy::Occupancy)
}

/// Select alternatives with an explicit search-order policy.
///
/// # Errors
///
/// Returns the validation and search errors from [`select_atom_conformers`].
pub fn select_atom_conformers_with_policy(
    atoms: &[ConformerAtom],
    excluded_bonds: &[[usize; 2]],
    seed: u64,
    policy: ConformerPolicy,
) -> Result<ConformerSelection, SemanticError> {
    select_with_links(atoms, excluded_bonds, &[], seed, policy)
}

pub(super) fn select_with_links(
    atoms: &[ConformerAtom],
    excluded_bonds: &[[usize; 2]],
    linked_alternatives: &[[usize; 2]],
    seed: u64,
    policy: ConformerPolicy,
) -> Result<ConformerSelection, SemanticError> {
    if atoms.iter().any(|a| {
        !a.occupancy.is_finite()
            || a.coordinates
                .iter()
                .any(|x| !x.is_finite() || x.abs() > 1e12)
    }) {
        return Err(fail("Conformer coordinates or occupancies are invalid"));
    }
    if excluded_bonds
        .iter()
        .any(|b| b.iter().any(|&i| i >= atoms.len()))
    {
        return Err(fail("Bond exclusion index is out of range"));
    }
    let exclusions: HashSet<_> = excluded_bonds
        .iter()
        .map(|&[a, b]| canonical(a, b))
        .collect();
    let (mut groups, memberships) = groups(atoms);
    if groups.is_empty() {
        return Ok(ConformerSelection {
            schema_version: policy.schema(),
            atom_rows: (0..atoms.len()).collect(),
            assignments: BTreeMap::new(),
            search_nodes: 0,
            seed,
            model_number: None,
        });
    }
    let (mut constraints, preferences) = contacts(atoms, &memberships, &exclusions, &mut groups);
    for &[first, second] in linked_alternatives {
        if let (Some((g, v)), Some((h, w))) = (memberships[first], memberships[second]) {
            if g == h {
                continue;
            }
            groups[g].neighbors.insert(h);
            groups[h].neighbors.insert(g);
            for a in 0..groups[g].labels.len() {
                for b in 0..groups[h].labels.len() {
                    if (a == v) != (b == w) {
                        let (pair, variants) = if g < h {
                            ((g, h), (a, b))
                        } else {
                            ((h, g), (b, a))
                        };
                        constraints.entry(pair).or_default().insert(variants);
                    }
                }
            }
        }
    }
    let mut assigned = vec![None; groups.len()];
    let mut nodes = 0;
    let mut visited = HashSet::new();
    for start in 0..groups.len() {
        if visited.contains(&start) {
            continue;
        }
        let mut component = Vec::new();
        let mut queue = vec![start];
        while let Some(g) = queue.pop() {
            if visited.insert(g) {
                component.push(g);
                queue.extend(groups[g].neighbors.iter().copied());
            }
        }
        if component.len() > 512 {
            return Err(SemanticError::new(
                "PDBX_CONFORMER_LIMIT",
                "Correlated conformer component exceeds 512 residues",
                vec!["atom_site".into()],
            ));
        }
        component.sort_by_key(|&g| {
            (
                std::cmp::Reverse(groups[g].neighbors.len()),
                groups[g].residue,
            )
        });
        if !solve(
            &component,
            0,
            &groups,
            &constraints,
            &preferences,
            &mut assigned,
            &mut nodes,
            seed,
            policy,
        ) {
            if nodes > 100_000 {
                return Err(SemanticError::new(
                    "PDBX_CONFORMER_LIMIT",
                    "Conformer search exceeded 100000 choices",
                    vec!["atom_site".into()],
                ));
            }
            return Err(fail(format!(
                "No clash-free conformer assignment for residue {}",
                groups[start].residue
            )));
        }
    }
    let atom_rows = memberships
        .iter()
        .enumerate()
        .filter_map(|(i, member)| {
            member
                .is_none_or(|(g, v)| assigned[g] == Some(v))
                .then_some(i)
        })
        .collect();
    let assignments = groups
        .iter()
        .enumerate()
        .filter_map(|(g, group)| {
            assigned[g]
                .and_then(|v| group.labels.get(v))
                .map(|label| (group.residue, label.clone()))
        })
        .collect();
    Ok(ConformerSelection {
        schema_version: policy.schema(),
        atom_rows,
        assignments,
        search_nodes: nodes,
        seed,
        model_number: None,
    })
}

fn canonical(a: usize, b: usize) -> Pair {
    if a < b { (a, b) } else { (b, a) }
}

fn groups(atoms: &[ConformerAtom]) -> (Vec<Group>, Vec<Option<Pair>>) {
    let mut by_residue: BTreeMap<usize, BTreeMap<&str, f64>> = BTreeMap::new();
    for atom in atoms {
        if let Some(alt) = &atom.alternate {
            *by_residue
                .entry(atom.residue)
                .or_default()
                .entry(alt)
                .or_default() += atom.occupancy;
        }
    }
    let mut groups = Vec::new();
    let mut indices = HashMap::new();
    for (residue, values) in by_residue {
        if values.len() < 2 {
            continue;
        }
        let (labels, occupancy): (Vec<_>, Vec<_>) =
            values.into_iter().map(|(s, v)| (s.to_owned(), v)).unzip();
        indices.insert(residue, groups.len());
        groups.push(Group {
            residue,
            labels,
            occupancy,
            ..Group::default()
        });
    }
    let memberships = atoms
        .iter()
        .map(|a| {
            a.alternate.as_ref().and_then(|alt| {
                indices.get(&a.residue).and_then(|&g| {
                    groups[g]
                        .labels
                        .iter()
                        .position(|s| s == alt)
                        .map(|v| (g, v))
                })
            })
        })
        .collect();
    (groups, memberships)
}

#[allow(clippy::too_many_arguments)]
fn solve(
    component: &[usize],
    depth: usize,
    groups: &[Group],
    constraints: &Constraints,
    preferences: &HashSet<Pair>,
    assigned: &mut [Option<usize>],
    nodes: &mut usize,
    seed: u64,
    policy: ConformerPolicy,
) -> bool {
    let Some(&g) = component.get(depth) else {
        return true;
    };
    let mut options: Vec<_> = (0..groups[g].labels.len())
        .filter(|v| {
            !groups[g].blocked.contains(v)
                && (policy == ConformerPolicy::Occupancy || groups[g].occupancy[*v] > 0.)
        })
        .collect();
    let preference = |v: usize| {
        groups[g]
            .neighbors
            .iter()
            .filter(|&&h| {
                preferences.contains(&canonical(g, h))
                    && assigned[h].is_some_and(|w| groups[h].labels[w] == groups[g].labels[v])
            })
            .count()
    };
    options.sort_by(|&a, &b| {
        if policy == ConformerPolicy::RandomClashAware {
            return tie(seed ^ 0xa0761d6478bd642f, groups[g].residue, a).cmp(&tie(
                seed ^ 0xa0761d6478bd642f,
                groups[g].residue,
                b,
            ));
        }
        groups[g].occupancy[b]
            .total_cmp(&groups[g].occupancy[a])
            .then_with(|| preference(b).cmp(&preference(a)))
            .then_with(|| tie(seed, groups[g].residue, a).cmp(&tie(seed, groups[g].residue, b)))
    });
    for v in options {
        *nodes += 1;
        if *nodes > 100_000 {
            return false;
        }
        let blocked = groups[g].neighbors.iter().any(|&h| {
            assigned[h].is_some_and(|w| {
                let (pair, variants) = if g < h {
                    ((g, h), (v, w))
                } else {
                    ((h, g), (w, v))
                };
                constraints
                    .get(&pair)
                    .is_some_and(|pairs| pairs.contains(&variants))
            })
        });
        if blocked {
            continue;
        }
        assigned[g] = Some(v);
        if solve(
            component,
            depth + 1,
            groups,
            constraints,
            preferences,
            assigned,
            nodes,
            seed,
            policy,
        ) {
            return true;
        }
        assigned[g] = None;
    }
    false
}

fn tie(seed: u64, group: usize, variant: usize) -> u64 {
    if seed == 0 {
        return variant as u64;
    }
    let mut x = seed ^ (group as u64).wrapping_mul(0x9e3779b97f4a7c15) ^ variant as u64;
    x = (x ^ (x >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
    (x ^ (x >> 27)).wrapping_mul(0x94d049bb133111eb)
}
