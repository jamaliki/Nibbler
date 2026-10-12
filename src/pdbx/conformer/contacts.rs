//! Sparse contacts define conformer conflicts and cross-chain preferences.

use super::{ConformerAtom, Constraints, Group, Pair, canonical};
use std::collections::{HashMap, HashSet};

pub(super) fn contacts(
    atoms: &[ConformerAtom],
    membership: &[Option<Pair>],
    exclusions: &HashSet<Pair>,
    groups: &mut [Group],
) -> (Constraints, HashSet<Pair>) {
    let mut constraints: Constraints = HashMap::new();
    let mut preferences = HashSet::new();
    let mut cells: HashMap<[i64; 3], Vec<usize>> = HashMap::new();
    for (i, atom) in atoms
        .iter()
        .enumerate()
        .filter(|(_, a)| a.clash_eligible && a.occupancy > 0.)
    {
        cells
            .entry(atom.coordinates.map(|x| (x / 5.).floor() as i64))
            .or_default()
            .push(i);
    }
    for (i, atom) in atoms
        .iter()
        .enumerate()
        .filter(|(i, a)| membership[*i].is_some() && a.clash_eligible && a.occupancy > 0.)
    {
        let cell = atom.coordinates.map(|x| (x / 5.).floor() as i64);
        for dx in -1..=1 {
            for dy in -1..=1 {
                for dz in -1..=1 {
                    let Some(nearby) = cells.get(&[cell[0] + dx, cell[1] + dy, cell[2] + dz])
                    else {
                        continue;
                    };
                    for &j in nearby {
                        if i == j
                            || (membership[j].is_some() && j < i)
                            || atom.residue == atoms[j].residue
                        {
                            continue;
                        }
                        let distance: f64 = atom
                            .coordinates
                            .iter()
                            .zip(atoms[j].coordinates)
                            .map(|(x, y)| (x - y).powi(2))
                            .sum();
                        if distance > 25. {
                            continue;
                        }
                        let Some((g, v)) = membership[i] else {
                            continue;
                        };
                        if let Some((h, w)) = membership[j] {
                            if atom.chain != atoms[j].chain
                                && !groups[g]
                                    .labels
                                    .iter()
                                    .all(|label| !groups[h].labels.contains(label))
                            {
                                preferences.insert(canonical(g, h));
                                groups[g].neighbors.insert(h);
                                groups[h].neighbors.insert(g);
                            }
                            if distance < 2.25
                                && !exclusions.contains(&canonical(i, j))
                                && !backbone_bond(atom, &atoms[j])
                            {
                                let (pair, variants) = if g < h {
                                    ((g, h), (v, w))
                                } else {
                                    ((h, g), (w, v))
                                };
                                constraints.entry(pair).or_default().insert(variants);
                                groups[g].neighbors.insert(h);
                                groups[h].neighbors.insert(g);
                            }
                        } else if distance < 2.25
                            && !exclusions.contains(&canonical(i, j))
                            && !backbone_bond(atom, &atoms[j])
                        {
                            groups[g].blocked.insert(v);
                        }
                    }
                }
            }
        }
    }
    (constraints, preferences)
}

fn backbone_bond(a: &ConformerAtom, b: &ConformerAtom) -> bool {
    if a.chain != b.chain {
        return false;
    }
    match (
        a.atom_name.as_str(),
        b.atom_name.as_str(),
        a.sequence,
        b.sequence,
    ) {
        ("C", "N", Some(x), Some(y)) => x.checked_add(1) == Some(y),
        ("N", "C", Some(x), Some(y)) => y.checked_add(1) == Some(x),
        _ => false,
    }
}
