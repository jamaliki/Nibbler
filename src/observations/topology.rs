//! Deposited links and observed peptide continuity.

use super::residue::distance;
use super::source::Source;
use std::collections::HashMap;

pub(super) fn pdb(source: &Source) -> Vec<[i64; 3]> {
    let mut result = vec![[-1; 3]; source.groups.len()];
    let mut backbone = HashMap::new();
    for (group, rows) in source.groups.iter().enumerate() {
        let mut sites = [None; 3];
        for &row in rows {
            let atom = &source.atoms[row];
            match atom.name.as_str() {
                "N" => sites[0] = Some(atom.position),
                "CA" => sites[1] = Some(atom.position),
                "C" => sites[2] = Some(atom.position),
                "OXT" => result[group][1] = 1,
                _ => {}
            }
        }
        let first = &source.atoms[rows[0]];
        if first.sequence >= 1 && sites.iter().all(Option::is_some) {
            backbone.insert((&first.chain, first.sequence), (group, sites));
        }
    }
    for &(a, b) in &source.bonds {
        for (here, other) in [(a, b), (b, a)] {
            let atom = &source.atoms[here];
            match (atom.name.as_str(), source.atoms[other].name.as_str()) {
                ("SG", _) => result[atom.group][2] = 1,
                ("N", "C") if !super::continuity::distant(source, atom.group, source.atoms[other].group, true)
                    => result[atom.group][0] = 0,
                ("C", "N") if !super::continuity::distant(source, atom.group, source.atoms[other].group, true)
                    => result[atom.group][1] = 0,
                _ => {}
            }
        }
    }
    for ((chain, index), (group, sites)) in &backbone {
        if let Some((next, following)) = backbone.get(&(chain, index + 1))
            && !super::continuity::distant(source, *group, *next, true)
            && let (Some(c), Some(n)) = (sites[2], following[0])
        {
            let value = distance(c, n);
            if value > 0.8 && value < 1.8 {
                result[*group][1] = 0;
                result[*next][0] = 0;
            }
        }
    }
    for (group, rows) in source.groups.iter().enumerate() {
        if result[group][0] < 0 && source.atoms[rows[0]].sequence == 1 {
            result[group][0] = 1;
        }
    }
    result
}

pub(super) fn protein(source: &Source, policy: &str) -> Result<Vec<[i64; 3]>, String> {
    if !matches!(policy, "entity_sequence" | "observed_connectivity") {
        return Err("Unknown topology policy".into());
    }
    let mut result = vec![[-1; 3]; source.groups.len()];
    let mut backbone = HashMap::new();
    for (group, rows) in source.groups.iter().enumerate() {
        let mut sites = [None; 2];
        for &row in rows {
            let atom = &source.atoms[row];
            if atom.name == "N" {
                sites[0] = Some(atom.position);
            }
            if atom.name == "C" {
                sites[1] = Some(atom.position);
            }
        }
        let first = &source.atoms[rows[0]];
        backbone.insert((&first.chain, first.sequence), (sites, group));
        if policy == "entity_sequence" {
            result[group][0] = i64::from(first.sequence == 1);
            if let Some(&length) = source.lengths.get(&first.entity) {
                result[group][1] = i64::from(first.sequence == length);
            }
        }
    }
    if policy == "observed_connectivity" {
        for (group, rows) in source.groups.iter().enumerate() {
            let first = &source.atoms[rows[0]];
            let (here, _) = backbone[&(&first.chain, first.sequence)];
            for (end, step, a, b) in [(0, -1, 0, 1), (1, 1, 1, 0)] {
                if let Some((neighbour, other)) = backbone.get(&(&first.chain, first.sequence + step))
                    && let (Some(left), Some(right)) = (here[a], neighbour[b])
                {
                    let value = distance(left, right);
                    if value > 0.8 && value < 1.8
                        && !super::continuity::distant(source, group, *other, true) {
                        result[group][end] = 0;
                    }
                }
            }
        }
    }
    Ok(result)
}
