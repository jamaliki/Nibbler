use super::{
    PreparedChemistry, SemanticError,
    product::{displace_h, overvalent, valid},
};

pub(super) fn lower_pi(
    graph: &mut PreparedChemistry,
    a: usize,
    require_partner: bool,
) -> Result<bool, SemanticError> {
    let neighbors: Vec<_> = graph.neighbors(a).collect();
    for (b, order) in neighbors {
        if order < 2
            || !matches!(graph.atoms[b].element.as_str(), "C" | "N")
            || graph.atoms[a].residue != graph.atoms[b].residue
        {
            continue;
        }
        if !overvalent(graph, b) {
            if require_partner || (graph.atoms[b].element == "C" && order != 2) {
                continue;
            }
            let donor = graph.neighbors(a).find_map(|(d, o)| {
                (o == 1
                    && graph.atoms[d].residue != graph.atoms[a].residue
                    && matches!(graph.atoms[d].element.as_str(), "N" | "O" | "S")
                    && !graph.hydrogens(d).is_empty()
                    && valid(
                        graph,
                        d,
                        graph.degree(d).saturating_sub(1),
                        graph.atoms[d].formal_charge,
                    ))
                .then_some(d)
            });
            let Some(donor) = donor else { continue };
            let Some(h) = graph.hydrogens(donor).first().copied() else {
                continue;
            };
            let mut new_h = graph.atoms[h].clone();
            new_h.atom_id = format!("H@{}", graph.atoms[b].atom_id);
            new_h.template_atom_id = None;
            new_h.deposited_charge = None;
            new_h.residue = graph.atoms[b].residue;
            new_h.component_id = graph.atoms[b].component_id.clone();
            new_h.asym_id = graph.atoms[b].asym_id.clone();
            new_h.sequence_id = graph.atoms[b].sequence_id;
            new_h.auth_sequence_id = graph.atoms[b].auth_sequence_id.clone();
            new_h.insertion_code = graph.atoms[b].insertion_code.clone();
            new_h.alternate_id = graph.atoms[b].alternate_id.clone();
            new_h.source_row = None;
            new_h.coordinates = None;
            new_h.leaving = false;
            new_h.impacted = false;
            let index = graph.atoms.len();
            graph.atoms.push(new_h);
            graph.adjacency.push(Vec::new());
            graph.remove(&[h], "hydrogen_transfer");
            graph.add_bond(b, index, 1, super::BondSource::HydrogenTransfer)?;
            graph.atoms[donor].impacted = true;
        }
        graph.change_order(a, b, order - 1, "pi_bond_addition");
        return Ok(true);
    }
    Ok(false)
}

pub(super) fn carbon_product(graph: &mut PreparedChemistry, a: usize) -> Result<(), SemanticError> {
    if displace_h(graph, a, 1) > 0 {
        return Ok(());
    }
    let neighbors: Vec<_> = graph.neighbors(a).collect();
    let carbonyl = neighbors
        .iter()
        .any(|&(b, order)| order == 2 && graph.atoms[b].element == "O");
    if carbonyl {
        for &(b, order) in &neighbors {
            if order == 1
                && graph.atoms[b].element == "N"
                && graph.atoms[b].residue == graph.atoms[a].residue
                && graph.atoms[b].coordinates.is_none()
                && graph
                    .neighbors(b)
                    .filter(|&(n, _)| !graph.atoms[n].is_hydrogen())
                    .count()
                    == 1
            {
                let mut removed = graph.hydrogens(b);
                removed.push(b);
                graph.remove(&removed, "unresolved_amide_substitution");
                return Ok(());
            }
        }
    }
    for &(b, order) in &neighbors {
        if order == 2 && graph.atoms[b].element == "O" {
            graph.change_order(a, b, 1, "carbonyl_addition");
            graph.set_charge(b, -1, "carbonyl_addition");
            return Ok(());
        }
    }
    lower_pi(graph, a, false)?;
    Ok(())
}

pub(super) fn rebalance_terminals(
    graph: &mut PreparedChemistry,
    linked: &[usize],
) -> Result<(), SemanticError> {
    for &a in linked {
        if !overvalent(graph, a) || !matches!(graph.atoms[a].element.as_str(), "O" | "N" | "S") {
            continue;
        }
        let mut candidates = Vec::new();
        for (center, order) in graph.neighbors(a) {
            if order != 2 || graph.atoms[center].residue != graph.atoms[a].residue {
                continue;
            }
            for (b, order) in graph.neighbors(center) {
                if b == a
                    || order != 1
                    || graph.atoms[b].residue != graph.atoms[a].residue
                    || !matches!(graph.atoms[b].element.as_str(), "O" | "N" | "S")
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
                let charge = if graph.atoms[b].formal_charge == Some(-1) {
                    Some(0)
                } else {
                    graph.atoms[b].formal_charge
                };
                let degree = graph.degree(b) + 1 - usize::from(has_h);
                if (has_h || graph.atoms[b].formal_charge == Some(-1))
                    && valid(graph, b, degree, charge)
                {
                    let metal = graph.adjacency[b]
                        .iter()
                        .any(|&k| graph.bonds[k].active && graph.bonds[k].order == 0);
                    candidates.push((
                        metal,
                        graph.atoms[b].element.clone(),
                        charge,
                        degree,
                        has_h,
                        center,
                        b,
                    ));
                }
            }
        }
        candidates.sort();
        if let Some(&(metal, ref element, charge, degree, has_h, center, b)) = candidates.first() {
            if candidates
                .iter()
                .filter(|c| c.0 == metal)
                .any(|c| c.1 != *element || c.2 != charge || c.3 != degree)
            {
                continue;
            }
            if has_h {
                displace_h(graph, b, 1);
            }
            graph.change_order(a, center, 1, "terminal_resonance_shift");
            graph.change_order(center, b, 2, "terminal_resonance_shift");
            if let Some(charge) = charge {
                graph.set_charge(b, charge, "terminal_resonance_shift");
            }
        }
    }
    // Phosphoryl addition can lower a metal-coordinated P=O, but must not
    // neutralize an arbitrary overvalent phosphorus product.
    for &p in linked {
        if graph.atoms[p].element != "P" || graph.degree(p) != 6 {
            continue;
        }
        let oxygen = graph.neighbors(p).find_map(|(o, order)| {
            (order == 2
                && graph.atoms[o].element == "O"
                && graph.adjacency[o]
                    .iter()
                    .any(|&k| graph.bonds[k].active && graph.bonds[k].order == 0))
            .then_some(o)
        });
        if let Some(o) = oxygen {
            graph.change_order(p, o, 1, "coordinated_phosphoryl_addition");
            graph.set_charge(o, -1, "coordinated_phosphoryl_addition");
        }
    }
    Ok(())
}
