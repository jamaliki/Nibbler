use super::super::{
    PdbxModel,
    fields::{category_rows, optional_float, optional_text},
};
use super::input::order;
use super::{PreparedChemistry, SemanticError};
use crate::cif::CategoryView;
use std::collections::{BTreeMap, HashMap};

type EndpointIndex = HashMap<String, BTreeMap<Option<i64>, Vec<usize>>>;

pub(super) fn add(model: &PdbxModel, graph: &mut PreparedChemistry) -> Result<(), SemanticError> {
    let Some(block) = model.document.blocks().first() else {
        return Ok(());
    };
    let mut index: EndpointIndex = HashMap::new();
    for (r, residue) in graph.residues.iter().enumerate() {
        if let Some(&a) = residue.atoms.values().next() {
            index
                .entry(graph.atoms[a].asym_id.clone())
                .or_default()
                .entry(graph.atoms[a].sequence_id)
                .or_default()
                .push(r);
        }
    }
    for row in category_rows(
        CategoryView::new(block.entries()),
        "struct_conn",
        [
            "id",
            "conn_type_id",
            "ptnr1_label_asym_id",
            "ptnr1_label_comp_id",
            "ptnr1_label_seq_id",
            "ptnr1_label_atom_id",
            "ptnr2_label_asym_id",
            "ptnr2_label_comp_id",
            "ptnr2_label_seq_id",
            "ptnr2_label_atom_id",
            "pdbx_ptnr1_label_alt_id",
            "pdbx_ptnr2_label_alt_id",
            "ptnr1_symmetry",
            "ptnr2_symmetry",
            "pdbx_value_order",
            "pdbx_dist_value",
            "pdbx_role",
            "ptnr1_auth_seq_id",
            "ptnr2_auth_seq_id",
            "pdbx_ptnr1_PDB_ins_code",
            "pdbx_ptnr2_PDB_ins_code",
            "ptnr1_auth_asym_id",
            "ptnr2_auth_asym_id",
            "ptnr1_auth_comp_id",
            "ptnr2_auth_comp_id",
        ],
    ) {
        let values: Vec<_> = (0..25).map(|i| optional_text(&row, i)).collect();
        let kind = values[1].as_deref().unwrap_or("").to_ascii_lowercase();
        if !matches!(kind.as_str(), "covale" | "disulf" | "metalc") {
            continue;
        }
        let a = endpoint(graph, &index, &values, false)?;
        let b = endpoint(graph, &index, &values, true)?;
        let (Some(a), Some(b)) = (a, b) else {
            graph.diagnostics.push(format!(
                "Excluded conformer-specific link {}",
                values[0].as_deref().unwrap_or("?")
            ));
            continue;
        };
        let image = [&values[12], &values[13]]
            .iter()
            .any(|value| value.as_deref().is_some_and(|s| s != "1_555"));
        if image {
            if kind != "metalc" {
                return Err(graph.error(
                    a,
                    "Nonidentity covalent link requires explicit assembly coordinates",
                ));
            }
            graph.diagnostics.push(format!(
                "Excluded symmetry-image coordination {}",
                values[0].as_deref().unwrap_or("?")
            ));
            continue;
        }
        let distance = optional_float(&row, "struct_conn", "pdbx_dist_value", 15)?;
        let glyco = values[16].as_deref().is_some_and(|r| {
            matches!(
                r.to_ascii_lowercase().as_str(),
                "n-glycosylation" | "o-glycosylation" | "s-glycosylation" | "c-mannosylation"
            )
        }) || [a, b].iter().all(|&x| {
            graph.residues[graph.atoms[x].residue]
                .component_type
                .as_deref()
                .is_some_and(|t| t.to_ascii_lowercase().contains("saccharide"))
        });
        if kind != "metalc" && too_long(graph, a, b, distance, glyco) {
            continue;
        }
        let bond_order = if kind == "metalc" {
            0
        } else if let Some(value) = &values[14] {
            order(value)?
        } else {
            infer_order(graph, a, b)?
        };
        graph.add_bond(a, b, bond_order, super::BondSource::StructConn)?;
        if let Some(k) = graph.bond_index(a, b) {
            if let Some(id) = &values[0] {
                graph.bonds[k].connection_ids.push(id.clone());
            }
            if values[14].is_some() {
                graph.bonds[k].declared_order = Some(bond_order);
            }
        }
        if bond_order > 0 && graph.atoms[a].residue != graph.atoms[b].residue {
            graph.atoms[a].impacted = true;
            graph.atoms[b].impacted = true;
        }
    }
    let polymers = model
        .entities
        .iter()
        .filter(|e| e.kind == super::super::model::EntityKind::Polymer)
        .map(|e| e.id.as_str())
        .collect();
    add_polymer_links(graph, &polymers)
}

fn endpoint(
    graph: &PreparedChemistry,
    index: &EndpointIndex,
    values: &[Option<String>],
    second: bool,
) -> Result<Option<usize>, SemanticError> {
    let offset = if second { 6 } else { 2 };
    let alt = values[if second { 11 } else { 10 }].as_deref();
    let author = values[if second { 18 } else { 17 }].as_deref();
    let insertion = values[if second { 20 } else { 19 }].as_deref();
    let chain = values[offset].as_deref().and_then(|name| index.get(name));
    let number = values[offset + 2]
        .as_deref()
        .and_then(|s| s.parse::<i64>().ok());
    let residues: Vec<_> = if let Some(chain) = chain {
        if number.is_some() {
            chain.get(&number).into_iter().flatten().copied().collect()
        } else {
            chain.values().flatten().copied().collect()
        }
    } else {
        Vec::new()
    };
    let candidates: Vec<_> = residues
        .iter()
        .filter_map(|&r| {
            values[offset + 3]
                .as_ref()
                .and_then(|name| {
                    graph.residues[r]
                        .atoms
                        .get(name)
                        .or_else(|| graph.residues[r].aliases.get(name))
                })
                .map(|&a| (a, &graph.atoms[a]))
        })
        .filter(|(_, a)| {
            values[offset + 1]
                .as_deref()
                .is_none_or(|c| c == a.component_id)
                && author.is_none_or(|s| s == a.auth_sequence_id)
                && insertion.is_none_or(|s| a.insertion_code.as_deref() == Some(s))
                && values[if second { 22 } else { 21 }]
                    .as_deref()
                    .is_none_or(|s| graph.residues[a.residue].auth_chain == s)
                && values[if second { 24 } else { 23 }]
                    .as_deref()
                    .is_none_or(|s| graph.residues[a.residue].auth_component == s)
        })
        .collect();
    let selected: Vec<_> = candidates
        .iter()
        .filter(|(_, a)| {
            alt.is_none_or(|s| graph.residues[a.residue].selected_alt.as_deref() == Some(s))
        })
        .collect();
    if selected.is_empty() && alt.is_some() && !candidates.is_empty() {
        return Ok(None);
    }
    let [entry] = selected.as_slice() else {
        return Err(SemanticError::new(
            "PDBX_LINK_PRODUCT",
            "Authored link has an absent or ambiguous endpoint",
            vec!["struct_conn".into(), values[0].clone().unwrap_or_default()],
        ));
    };
    Ok(Some(entry.0))
}

fn too_long(
    graph: &mut PreparedChemistry,
    a: usize,
    b: usize,
    authored: Option<f64>,
    glyco: bool,
) -> bool {
    let expected = if [a, b]
        .iter()
        .all(|&x| matches!(graph.atoms[x].element.as_str(), "C" | "H" | "N" | "O"))
    {
        1.7
    } else if [a, b].iter().all(|&x| {
        matches!(
            graph.atoms[x].element.as_str(),
            "C" | "H" | "N" | "O" | "P" | "S"
        )
    }) {
        2.4
    } else {
        3.6
    };
    let limit = if glyco { 2.4 } else { expected };
    // Deposited distances describe the authored state; do not substitute a
    // different conformer's coordinates when no distance was supplied.
    if let Some(distance) = authored.filter(|&d| d > expected) {
        let kept = distance <= limit;
        graph.diagnostics.push(format!(
            "{} long authored link at {} ({} A)",
            if kept { "Kept" } else { "Excluded" },
            graph.atoms[a].atom_id,
            distance
        ));
        return !kept;
    }
    false
}

fn infer_order(graph: &PreparedChemistry, a: usize, b: usize) -> Result<u8, SemanticError> {
    if let Some(k) = graph.bond_index(a, b) {
        return Ok(graph.bonds[k].order);
    }
    let leaving_order = |x| {
        graph
            .neighbors(x)
            .filter(|&(n, _)| graph.atoms[n].leaving)
            .map(|(_, o)| o)
            .max()
            .unwrap_or(1)
    };
    let mut value = leaving_order(a).max(leaving_order(b));
    for x in [a, b] {
        if graph.atoms[x].hydrogen_count.is_none() {
            return Err(graph.error(x, "Missing linked component template"));
        }
        let replaceable: usize = graph
            .neighbors(x)
            .filter(|&(n, _)| graph.atoms[n].leaving || graph.atoms[n].is_hydrogen())
            .map(|(_, o)| usize::from(o))
            .sum();
        let free = super::valence::neutral_degree(&graph.atoms[x].element)
            .unwrap_or(graph.degree(x))
            .saturating_sub(graph.degree(x));
        if value > 1 {
            value = value.min(u8::try_from(replaceable + free).unwrap_or(u8::MAX));
        }
    }
    if value == 0 {
        return Err(graph.error(a, "Cannot infer a multiple bond supported by both partners"));
    }
    Ok(value)
}

fn add_polymer_links(
    graph: &mut PreparedChemistry,
    polymers: &std::collections::HashSet<&str>,
) -> Result<(), SemanticError> {
    let mut sequence = std::collections::BTreeMap::new();
    for (r, residue) in graph.residues.iter().enumerate() {
        if !polymers.contains(residue.entity_id.as_str()) {
            continue;
        }
        let Some(&first) = residue.atoms.values().next() else {
            continue;
        };
        let a = &graph.atoms[first];
        if let Some(number) = a.sequence_id {
            sequence.insert((a.asym_id.clone(), number), r);
        }
    }
    for ((chain, number), &left) in &sequence {
        let Some(next) = number.checked_add(1) else {
            continue;
        };
        let Some(&right) = sequence.get(&(chain.clone(), next)) else {
            continue;
        };
        if graph.residues[left].entity_id != graph.residues[right].entity_id {
            continue;
        }
        let kind = graph.residues[left]
            .component_type
            .as_deref()
            .unwrap_or("")
            .to_ascii_lowercase();
        let (out, incoming) = if kind.contains("peptide") {
            ("C", "N")
        } else if kind.contains("rna") || kind.contains("dna") {
            ("O3'", "P")
        } else {
            continue;
        };
        let (Some(&a), Some(&b)) = (
            graph.residues[left].atoms.get(out),
            graph.residues[right].atoms.get(incoming),
        ) else {
            continue;
        };
        if kind.contains("peptide") {
            let ca = |residue: usize| graph.residues[residue].atoms.get("CA")
                .and_then(|&atom| graph.atoms[atom].coordinates);
            if let (Some(p), Some(q)) = (ca(left), ca(right))
                && p.iter().zip(q).map(|(x, y)| (x - y).powi(2)).sum::<f64>() > 4_f64.powi(2)
            {
                graph.diagnostics.push(format!("Excluded distant CA sequence link {chain}:{number}"));
                continue;
            }
        }
        if graph.bond_index(a, b).is_some() {
            continue;
        }
        if let (Some(p), Some(q)) = (graph.atoms[a].coordinates, graph.atoms[b].coordinates)
            && p.iter().zip(q).map(|(x, y)| (x - y).powi(2)).sum::<f64>() > 2.4_f64.powi(2)
        {
            graph.diagnostics.push(format!(
                "Excluded distant sequence link {}:{}",
                chain, number
            ));
            continue;
        }
        if graph
            .neighbors(a)
            .any(|(n, _)| graph.atoms[n].residue == right)
            || graph
                .neighbors(b)
                .any(|(n, _)| graph.atoms[n].residue == left)
        {
            continue;
        }
        for x in [a, b] {
            let can_displace = graph
                .neighbors(x)
                .any(|(n, _)| graph.atoms[n].leaving || graph.atoms[n].is_hydrogen());
            let full = super::valence::neutral_degree(&graph.atoms[x].element)
                .is_some_and(|d| graph.degree(x) >= d);
            if full && !can_displace {
                return Err(graph.error(x, "Sequence does not justify a saturated backbone link"));
            }
        }
        graph.add_bond(a, b, 1, super::BondSource::PolymerSequence)?;
        graph.atoms[a].impacted = true;
        graph.atoms[b].impacted = true;
    }
    Ok(())
}
