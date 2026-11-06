//! Alternate letters carried by different residues that occupy one physical site.
//!
//! A residue whose labelled atoms share one letter is not a variable group, so the
//! residue-level search never weighs it against anything. Depositions that store
//! the competing conformers of a site as separate residues, ligands or chains
//! (one letter each) therefore kept every conformer, and the view held the site
//! twice. This pass finds those sites and keeps one letter per site.
//!
//! The rule mirrors the residue-level one. Two single-letter residues with
//! different letters are alternatives of one site when a heavy atom of one lies
//! closer than 1.5 A to a heavy atom of the other, ignoring authored bonds,
//! backbone C-N links and unoccupied atoms. Alternatives that touch, directly or
//! through a chain of such contacts, form one site. The site keeps the letter with
//! the highest summed occupancy over its residues; the seed breaks ties and orders
//! the random policy. Every labelled atom of a losing letter in the site is
//! dropped, unlabelled atoms stay, and no coordinate changes. Residues with
//! several letters remain with the residue-level search.

use super::contacts::backbone_bond;
use super::{ConformerAtom, ConformerPolicy, Group, Pair, canonical, solve};
use std::collections::{BTreeMap, HashMap, HashSet};

/// Atoms removed because their site keeps a different letter.
pub(super) struct Duplicates {
    /// One flag per input atom: true when the atom is dropped.
    pub dropped: Vec<bool>,
    /// Choices examined while picking a letter for each site.
    pub nodes: usize,
}

/// Squared distance below which different letters cannot be coexisting conformers.
const CLASH: f64 = 2.25;

pub(super) fn resolve(
    atoms: &[ConformerAtom],
    exclusions: &HashSet<Pair>,
    seed: u64,
    policy: ConformerPolicy,
) -> Duplicates {
    let mut result = Duplicates {
        dropped: vec![false; atoms.len()],
        nodes: 0,
    };
    let single = single_letter_residues(atoms);
    if single.is_empty() {
        return result;
    }
    let sites = sites(atoms, &single, exclusions);
    if sites.is_empty() {
        return result;
    }
    let mut totals: BTreeMap<usize, BTreeMap<&str, f64>> = BTreeMap::new();
    for atom in atoms {
        if let (Some(letter), Some(&site)) = (&atom.alternate, sites.get(&atom.residue)) {
            *totals.entry(site).or_default().entry(letter).or_default() += atom.occupancy;
        }
    }
    let mut groups = Vec::new();
    let mut index = HashMap::new();
    for (site, letters) in totals {
        let (labels, occupancy): (Vec<_>, Vec<_>) = letters
            .into_iter()
            .map(|(label, value)| (label.to_owned(), value))
            .unzip();
        index.insert(site, groups.len());
        groups.push(Group {
            residue: site,
            labels,
            occupancy,
            ..Group::default()
        });
    }
    let mut assigned = vec![None; groups.len()];
    for g in 0..groups.len() {
        let solved = solve(
            &[g],
            0,
            &groups,
            &HashMap::new(),
            &HashSet::new(),
            &mut assigned,
            &mut result.nodes,
            seed,
            policy,
        );
        if !solved {
            assigned[g] = None;
        }
    }
    for (i, atom) in atoms.iter().enumerate() {
        let (Some(letter), Some(&site)) = (&atom.alternate, sites.get(&atom.residue)) else {
            continue;
        };
        let g = index[&site];
        if let Some(winner) = assigned[g] {
            result.dropped[i] = groups[g].labels[winner] != *letter;
        }
    }
    result
}

/// Residues whose labelled atoms all carry one letter, with that letter.
fn single_letter_residues(atoms: &[ConformerAtom]) -> HashMap<usize, &str> {
    let mut letters: HashMap<usize, Option<&str>> = HashMap::new();
    for atom in atoms {
        if let Some(letter) = atom.alternate.as_deref() {
            letters
                .entry(atom.residue)
                .and_modify(|seen| {
                    if *seen != Some(letter) {
                        *seen = None;
                    }
                })
                .or_insert(Some(letter));
        }
    }
    letters
        .into_iter()
        .filter_map(|(residue, letter)| letter.map(|letter| (residue, letter)))
        .collect()
}

/// Map each residue that touches a different letter to a site identifier.
///
/// The identifier is the smallest residue of the connected set, so it is stable.
fn sites(
    atoms: &[ConformerAtom],
    single: &HashMap<usize, &str>,
    exclusions: &HashSet<Pair>,
) -> HashMap<usize, usize> {
    let eligible = |a: &ConformerAtom| {
        a.alternate.is_some()
            && a.clash_eligible
            && a.occupancy > 0.
            && single.contains_key(&a.residue)
    };
    let cell_of = |atom: &ConformerAtom| atom.coordinates.map(|x| (x / 5.).floor() as i64);
    let mut cells: HashMap<[i64; 3], Vec<usize>> = HashMap::new();
    for (i, atom) in atoms.iter().enumerate().filter(|(_, a)| eligible(a)) {
        cells.entry(cell_of(atom)).or_default().push(i);
    }
    let mut parent: HashMap<usize, usize> = HashMap::new();
    for (i, atom) in atoms.iter().enumerate().filter(|(_, a)| eligible(a)) {
        let cell = cell_of(atom);
        for dx in -1..=1 {
            for dy in -1..=1 {
                for dz in -1..=1 {
                    let Some(nearby) = cells.get(&[cell[0] + dx, cell[1] + dy, cell[2] + dz])
                    else {
                        continue;
                    };
                    for &j in nearby.iter().filter(|&&j| j > i) {
                        let other = &atoms[j];
                        if atom.residue == other.residue
                            || atom.alternate == other.alternate
                            || exclusions.contains(&canonical(i, j))
                            || backbone_bond(atom, other)
                        {
                            continue;
                        }
                        let distance: f64 = atom
                            .coordinates
                            .iter()
                            .zip(other.coordinates)
                            .map(|(x, y)| (x - y).powi(2))
                            .sum();
                        if distance < CLASH {
                            union(&mut parent, atom.residue, other.residue);
                        }
                    }
                }
            }
        }
    }
    let members: Vec<usize> = parent.keys().copied().collect();
    members
        .into_iter()
        .map(|residue| (residue, find(&mut parent, residue)))
        .collect()
}

/// Root of a residue's set. Roots are always the smallest residue of the set.
fn find(parent: &mut HashMap<usize, usize>, residue: usize) -> usize {
    let mut root = residue;
    while let Some(&next) = parent.get(&root) {
        if next == root {
            break;
        }
        root = next;
    }
    let mut node = residue;
    while node != root {
        let next = parent[&node];
        parent.insert(node, root);
        node = next;
    }
    root
}

fn union(parent: &mut HashMap<usize, usize>, a: usize, b: usize) {
    parent.entry(a).or_insert(a);
    parent.entry(b).or_insert(b);
    let (a, b) = (find(parent, a), find(parent, b));
    if a != b {
        parent.insert(a.max(b), a.min(b));
    }
}
