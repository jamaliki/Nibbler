use super::{PreparedChemistry, SemanticError, rearrange, valence};
use std::collections::BTreeSet;

pub(super) fn resolve(graph: &mut PreparedChemistry) -> Result<(), SemanticError> {
    let mut linked: Vec<_> = graph
        .atoms
        .iter()
        .enumerate()
        .filter_map(|(i, a)| a.impacted.then_some(i))
        .collect();
    // Resolve carbon additions before donor H displacement; the complete
    // product may transfer that H to the former pi-bond partner.
    linked.sort_by_key(|&a| graph.atoms[a].element != "C");
    for &a in &linked {
        if graph.atoms[a].hydrogen_count.is_none() || graph.atoms[a].formal_charge.is_none() {
            return Err(graph.error(a, "Missing linked component chemistry"));
        }
    }
    remove_leaving_groups(graph, &linked)?;
    restore_lost_order(graph, &linked)?;
    for &a in &linked {
        if graph.atoms[a].element == "C" && overvalent(graph, a) {
            rearrange::lower_pi(graph, a, true)?;
        }
    }
    rearrange::rebalance_terminals(graph, &linked)?;
    for &a in &linked {
        if !graph.atoms[a].active {
            continue;
        }
        if graph.atoms[a].element == "S" && graph.degree(a) == 3 && displace_h(graph, a, 1) > 0 {
            continue;
        }
        if !overvalent(graph, a) {
            continue;
        }
        if graph.atoms[a].element == "C" {
            rearrange::carbon_product(graph, a)?;
        } else if matches!(graph.atoms[a].element.as_str(), "N" | "O") {
            let target = valence::neutral_degree(&graph.atoms[a].element).unwrap_or(0) as i32
                + graph.atoms[a].formal_charge.unwrap_or(0);
            let count = graph
                .degree(a)
                .saturating_sub(usize::try_from(target).unwrap_or(0));
            displace_h(graph, a, count);
        }
    }
    normalize_charges(graph);
    normalize_amides(graph)?;
    for a in 0..graph.atoms.len() {
        if graph.atoms[a].active
            && graph.atoms[a].impacted
            && !valid(graph, a, graph.degree(a), graph.atoms[a].formal_charge)
        {
            return Err(graph.error(
                a,
                &format!(
                    "Unsupported product valence {} and charge {:?}",
                    graph.degree(a),
                    graph.atoms[a].formal_charge
                ),
            ));
        }
        if graph.atoms[a].hydrogen_count.is_some() {
            graph.atoms[a].hydrogen_count = Some(if graph.atoms[a].active {
                graph
                    .neighbors(a)
                    .filter(|&(n, _)| graph.atoms[n].is_hydrogen())
                    .count()
            } else {
                0
            });
        }
    }
    Ok(())
}

pub(super) fn valid(
    graph: &PreparedChemistry,
    a: usize,
    degree: usize,
    charge: Option<i32>,
) -> bool {
    charge.is_some_and(|q| valence::accepted(&graph.atoms[a].element, q, degree))
}

pub(super) fn overvalent(graph: &PreparedChemistry, a: usize) -> bool {
    valence::neutral_degree(&graph.atoms[a].element).is_some_and(|base| {
        let charge = if graph.atoms[a].element == "C" {
            0
        } else {
            graph.atoms[a].formal_charge.unwrap_or(0)
        };
        let limit = usize::try_from(base as i32 + charge).unwrap_or(0);
        graph.degree(a) > limit
    })
}

pub(super) fn displace_h(graph: &mut PreparedChemistry, a: usize, count: usize) -> usize {
    let atoms: Vec<_> = graph.hydrogens(a).into_iter().take(count).collect();
    let removed = atoms.len();
    if removed > 0 {
        graph.remove(&atoms, "hydrogen_displacement");
    }
    removed
}

fn leaving_groups(graph: &PreparedChemistry, parent: usize) -> Vec<Vec<usize>> {
    let mut groups = Vec::new();
    let mut visited = BTreeSet::from([parent]);
    for (start, _) in graph.neighbors(parent) {
        if graph.atoms[start].residue != graph.atoms[parent].residue || visited.contains(&start) {
            continue;
        }
        let mut group = Vec::new();
        let mut queue = vec![start];
        while let Some(a) = queue.pop() {
            if !visited.insert(a) {
                continue;
            }
            group.push(a);
            for (b, _) in graph.neighbors(a) {
                if b != parent && graph.atoms[b].residue == graph.atoms[parent].residue {
                    queue.push(b);
                }
            }
        }
        let heavy: Vec<_> = group
            .iter()
            .copied()
            .filter(|&a| !graph.atoms[a].is_hydrogen())
            .collect();
        let flagged = if heavy.is_empty() { &group } else { &heavy };
        if !flagged.is_empty() && flagged.iter().all(|&a| graph.atoms[a].leaving) {
            groups.push(group);
        }
    }
    groups.sort_by_key(|g| g.iter().min().copied());
    groups
}

fn remove_leaving_groups(
    graph: &mut PreparedChemistry,
    linked: &[usize],
) -> Result<(), SemanticError> {
    for &a in linked {
        let inter: Vec<_> = graph.adjacency[a]
            .iter()
            .filter_map(|&k| {
                let b = &graph.bonds[k];
                let other = if b.atom_indices[0] == a {
                    b.atom_indices[1]
                } else {
                    b.atom_indices[0]
                };
                (b.active && b.order > 0 && graph.atoms[other].residue != graph.atoms[a].residue)
                    .then_some(k)
            })
            .collect();
        let mut budget: usize = inter
            .iter()
            .map(|&k| usize::from(graph.bonds[k].order))
            .sum();
        let polymer = inter.iter().any(|&k| {
            graph.bonds[k].source == super::BondSource::PolymerSequence || {
                let b = &graph.bonds[k];
                let names = (
                    &graph.atoms[b.atom_indices[0]].atom_id,
                    &graph.atoms[b.atom_indices[1]].atom_id,
                );
                matches!(
                    (names.0.as_str(), names.1.as_str()),
                    ("C", "N") | ("N", "C") | ("O3'", "P") | ("P", "O3'")
                )
            }
        });
        let mut groups = leaving_groups(graph, a);
        groups.sort_by_key(|g| g.iter().any(|&i| graph.atoms[i].coordinates.is_some()));
        for group in groups {
            if budget == 0 {
                break;
            }
            if group.iter().any(|&i| graph.atoms[i].impacted) {
                continue;
            }
            if !polymer
                && group
                    .iter()
                    .any(|&i| !graph.atoms[i].is_hydrogen() && graph.atoms[i].coordinates.is_some())
            {
                continue;
            }
            if group.iter().any(|&i| {
                graph
                    .neighbors(i)
                    .any(|(j, _)| j != a && !group.contains(&j))
            }) {
                return Err(graph.error(a, "Leaving group has an unexpected external bond"));
            }
            graph.remove(&group, "leaving_group_substitution");
            budget -= 1;
        }
    }
    Ok(())
}

fn restore_lost_order(
    graph: &mut PreparedChemistry,
    linked: &[usize],
) -> Result<(), SemanticError> {
    for &a in linked {
        let lost: usize = graph.adjacency[a]
            .iter()
            .filter_map(|&k| {
                let b = &graph.bonds[k];
                (!b.active)
                    .then_some(b.template_order)
                    .flatten()
                    .map(usize::from)
            })
            .sum();
        let gained: usize = graph
            .neighbors(a)
            .filter(|&(b, _)| graph.atoms[b].residue != graph.atoms[a].residue)
            .map(|(_, o)| usize::from(o))
            .sum();
        let deficit = lost.saturating_sub(gained);
        for _ in 0..deficit {
            let mut candidates = Vec::new();
            for (b, order) in graph.neighbors(a) {
                let atom = &graph.atoms[b];
                if order != 1
                    || atom.residue != graph.atoms[a].residue
                    || !matches!(atom.element.as_str(), "O" | "S" | "N")
                {
                    continue;
                }
                if graph
                    .neighbors(b)
                    .filter(|&(n, _)| !graph.atoms[n].is_hydrogen())
                    .count()
                    != 1
                {
                    continue;
                }
                let has_h = !graph.hydrogens(b).is_empty();
                let charge = if atom.formal_charge == Some(-1) {
                    Some(0)
                } else {
                    atom.formal_charge
                };
                let degree = graph.degree(b) + 1 - usize::from(has_h);
                if (has_h || atom.formal_charge == Some(-1)) && valid(graph, b, degree, charge) {
                    candidates.push((has_h, atom.element.clone(), charge, degree, b));
                }
            }
            candidates.sort();
            let Some(&(has_h, ref element, charge, degree, b)) = candidates.first() else {
                return Err(graph.error(a, "Unsupported leaving-group bond-order rearrangement"));
            };
            if candidates
                .iter()
                .any(|c| c.1 != *element || c.2 != charge || c.3 != degree)
            {
                return Err(graph.error(a, "Ambiguous leaving-group bond-order rearrangement"));
            }
            if has_h {
                displace_h(graph, b, 1);
            }
            graph.change_order(a, b, 2, "leaving_bond_order_restoration");
            if let Some(charge) = charge {
                graph.set_charge(b, charge, "leaving_bond_order_restoration");
            }
        }
    }
    Ok(())
}

fn normalize_charges(graph: &mut PreparedChemistry) {
    for a in 0..graph.atoms.len() {
        if !graph.atoms[a].active || !graph.atoms[a].impacted {
            continue;
        }
        let key = (
            graph.atoms[a].element.as_str(),
            graph.atoms[a].formal_charge,
            graph.degree(a),
        );
        let charge = match key {
            ("B", Some(0), 4) => Some(-1),
            ("C", Some(1), 4) => Some(0),
            ("N", Some(0), 4) | ("P", Some(0), 4) | ("S" | "SE", Some(0), 3) => Some(1),
            ("O", Some(-1), 2) => Some(0),
            ("O", Some(0), 1) | ("SE", Some(0), 1) => Some(-1),
            ("O", Some(0), 3)
                if graph
                    .neighbors(a)
                    .all(|(b, _)| graph.atoms[b].element != "O") =>
            {
                Some(1)
            }
            _ => None,
        };
        if let Some(charge) = charge {
            graph.set_charge(a, charge, "linked_product_charge");
        }
    }
}

fn normalize_amides(graph: &mut PreparedChemistry) -> Result<(), SemanticError> {
    for a in 0..graph.atoms.len() {
        if !graph.atoms[a].active
            || !graph.atoms[a].impacted
            || graph.atoms[a].element != "N"
            || graph.atoms[a].formal_charge != Some(1)
        {
            continue;
        }
        let amide = graph.neighbors(a).any(|(carbon, order)| {
            order == 1
                && graph.atoms[carbon].element == "C"
                && graph
                    .neighbors(carbon)
                    .any(|(oxygen, o)| o == 2 && graph.atoms[oxygen].element == "O")
        });
        if amide {
            if displace_h(graph, a, 1) == 0 {
                return Err(graph.error(a, "Charged amide nitrogen has no removable hydrogen"));
            }
            graph.set_charge(a, 0, "amide_neutralization");
        }
    }
    Ok(())
}
