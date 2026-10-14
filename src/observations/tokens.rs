//! Protein residue tokens and heavy ligand atom tokens.

use super::record::{Tokens, renumber_ligands};
use super::source::{Atom, Source};
use super::{HydrogenParents, LigandFeatures, Schema, residue, topology};
use std::collections::{BTreeSet, HashMap};

fn ligand_features(features: &LigandFeatures, atom: &Atom) -> [i16; 4] {
    features
        .get(&atom.component)
        .and_then(|component| component.get(&atom.name))
        .copied()
        .unwrap_or([0; 4])
}

fn hydrogen_overflow(
    source: &Source,
    features: &LigandFeatures,
    parent: usize,
    assigned: &[usize],
) -> String {
    let atom = &source.atoms[parent];
    let children: Vec<_> = assigned
        .iter()
        .map(|&row| {
            let child = &source.atoms[row];
            serde_json::json!({"name":child.name,"element":child.element,
            "position":child.position,"distance":residue::distance(atom.position,child.position)})
        })
        .collect();
    let ccd = features
        .get(&atom.component)
        .and_then(|component| component.get(&atom.name));
    let partners: Vec<_> = source
        .bonds
        .iter()
        .filter_map(|&(a, b)| {
            let other = if a == parent {
                b
            } else if b == parent {
                a
            } else {
                return None;
            };
            let other = &source.atoms[other];
            Some(
                serde_json::json!({"chain":other.chain,"component":other.component,
            "atom_name":other.name,"auth_sequence":other.auth_sequence,"position":other.position}),
            )
        })
        .collect();
    let detail = serde_json::json!({"parent_key":[atom.chain.clone(),atom.auth_sequence.clone(),
        atom.insertion.clone(),atom.component.clone(),atom.name.clone()],"parent_element":atom.element,
        "parent_position":atom.position,"capacity":3,"hydrogens":children,"covalent_partners":partners,
        "ccd_features":ccd.map(|values| serde_json::json!({"formal_charge":values[0],
            "hybridization":values[1],"total_hydrogen":values[2],"formal_charge_known":values[3] != 0}))});
    format!("Ligand parent has more than three hydrogens; unsupported Atom32 ligand: {detail}")
}

fn parent_error(atom: &Atom, parent: &str, reason: &str) -> String {
    let detail = serde_json::json!({"hydrogen_key":[atom.chain.clone(),atom.auth_sequence.clone(),
        atom.insertion.clone(),atom.component.clone(),atom.name.clone()],
        "position":atom.position,"parent":parent,"reason":reason});
    format!("Invalid ligand hydrogen parent: {detail}")
}

fn declared_parents(source: &Source) -> Result<HashMap<usize, Option<usize>>, String> {
    let mut result = HashMap::new();
    for &(child, parent) in &source.hydrogen_bonds {
        if let Some(previous) = result.insert(child, parent)
            && previous != parent
        {
            return Err(parent_error(
                &source.atoms[child],
                "struct_conn",
                "conflicting covalent parents",
            ));
        }
    }
    Ok(result)
}

fn ccd_parent<'a>(parents: &'a HashMap<String, String>, atom: &Atom) -> Option<&'a str> {
    if let Some(parent) = parents.get(&atom.name) {
        return Some(parent);
    }
    if matches!(atom.element.as_str(), "D" | "T") {
        let alias = atom.name.replacen(&atom.element, "H", 1);
        return parents.get(&alias).map(String::as_str);
    }
    None
}

fn ligand_parent(
    source: &Source,
    parents: &HydrogenParents,
    declared: &HashMap<usize, Option<usize>>,
    row: usize,
    rows: &[usize],
    heavy: &[(usize, [f32; 3])],
) -> Result<Option<usize>, String> {
    let atom = &source.atoms[row];
    let parent = if let Some(&parent) = declared.get(&row) {
        parent.ok_or_else(|| parent_error(atom, "struct_conn", "bonded parent is absent"))?
    } else if let Some(name) = parents
        .get(&atom.component)
        .and_then(|map| ccd_parent(map, atom))
    {
        rows.iter()
            .copied()
            .find(|&i| source.atoms[i].name == name && !residue::hydrogen(&source.atoms[i].element))
            .ok_or_else(|| parent_error(atom, name, "Dictionary parent is absent"))?
    } else {
        return Ok(residue::nearest(atom.position, heavy));
    };
    let target = &source.atoms[parent];
    let distance = residue::distance(atom.position, target.position);
    if !(residue::finite(atom.position)
        && residue::finite(target.position)
        && distance > 0.5
        && distance < 1.6
        && !target.protein)
    {
        return Err(parent_error(
            atom,
            &target.name,
            "bonded parent has invalid geometry or is a protein atom",
        ));
    }
    Ok(Some(parent))
}

fn add_proteins(
    source: &Source,
    schema: &Schema,
    complete: bool,
    protein_only: bool,
    policy: &str,
    tokens: &mut Tokens,
    mapping: &mut [usize],
) -> Result<(), String> {
    let topology = if protein_only {
        topology::protein(source, policy)?
    } else {
        topology::pdb(source)
    };
    let mut groups: Vec<_> = source
        .groups
        .iter()
        .enumerate()
        .filter_map(|(group, rows)| {
            rows.iter()
                .copied()
                .find(|&row| {
                    source.atoms[row].protein
                        && source.atoms[row].name == "N"
                        && residue::finite(source.atoms[row].position)
                })
                .map(|center| (center, group))
        })
        .collect();
    if protein_only {
        groups = source
            .groups
            .iter()
            .enumerate()
            .filter_map(|(group, rows)| {
                rows.iter()
                    .copied()
                    .find(|&row| source.atoms[row].protein)
                    .map(|first| (first, group))
            })
            .collect();
    }
    groups.sort_unstable();
    for (_, group) in groups {
        let rows = &source.groups[group];
        let first = &source.atoms[rows[0]];
        let aa = schema.identity(&first.component);
        if protein_only && aa >= 20 {
            continue;
        }
        let atoms: Vec<_> = rows.iter().map(|&row| &source.atoms[row]).collect();
        let mut packed = if aa < 20 {
            residue::pack(
                &schema.residues[aa],
                &atoms,
                topology[group],
                (!protein_only && topology[group][2] == 1).then_some(true),
                complete,
                true,
            )
        } else {
            residue::Packed::empty([-1; 3])
        };
        if aa >= 20 {
            for atom in &atoms {
                if let Some(slot) = ["N", "CA", "C", "O"]
                    .iter()
                    .position(|&name| name == atom.name)
                    && residue::finite(atom.position)
                {
                    packed.position[slot] = atom.position;
                    packed.mask[slot] = true;
                }
            }
        }
        if protein_only && !packed.mask[1] {
            continue;
        }
        if !protein_only && first.sequence < 1 {
            return Err("Protein token lacks a positive label_seq_id".into());
        }
        let token = tokens.identity.len();
        for &row in rows {
            let atom = &source.atoms[row];
            if !matches!(atom.element.as_str(), "H" | "D") && residue::finite(atom.position) {
                mapping[row] = token;
            }
        }
        tokens.push(first, aa, &packed, false, [0; 4]);
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn add_ligands(
    source: &Source,
    features: &LigandFeatures,
    parents: &HydrogenParents,
    elements: &[String],
    excluded: &[String],
    tokens: &mut Tokens,
    mapping: &mut [usize],
    complete: bool,
) -> Result<(), String> {
    let declared = declared_parents(source)?;
    let mut assigned: HashMap<usize, Vec<usize>> = HashMap::new();
    for (row, atom) in source.atoms.iter().enumerate() {
        if atom.protein
            || residue::hydrogen(&atom.element)
            || !residue::finite(atom.position)
            || excluded.contains(&atom.component)
        {
            continue;
        }
        mapping[row] = tokens.identity.len();
        let identity = elements
            .iter()
            .position(|element| element == &atom.element)
            .map_or(64, |index| 21 + index);
        let mut packed = residue::Packed::empty([-1; 3]);
        packed.position[0] = atom.position;
        packed.mask[0] = true;
        tokens.push(
            atom,
            identity,
            &packed,
            true,
            ligand_features(features, atom),
        );
        if complete {
            *tokens
                .hydrogen_known
                .last_mut()
                .ok_or("Missing ligand token")? = 1;
        }
    }
    for rows in &source.groups {
        let heavy: Vec<_> = rows
            .iter()
            .copied()
            .filter(|&row| !residue::hydrogen(&source.atoms[row].element))
            .map(|row| (row, source.atoms[row].position))
            .collect();
        for &row in rows.iter().filter(|&&row| {
            !source.atoms[row].protein
                && residue::hydrogen(&source.atoms[row].element)
                && !excluded.contains(&source.atoms[row].component)
        }) {
            let parent = ligand_parent(source, parents, &declared, row, rows, &heavy)?;
            if let Some(parent) = parent
                && mapping[parent] == usize::MAX
                && (declared.contains_key(&row)
                    || parents
                        .get(&source.atoms[row].component)
                        .and_then(|map| ccd_parent(map, &source.atoms[row]))
                        .is_some())
            {
                return Err(parent_error(
                    &source.atoms[row],
                    &source.atoms[parent].name,
                    "bonded parent has no ligand token",
                ));
            }
            let parent = parent.filter(|&parent| mapping[parent] != usize::MAX);
            let Some(parent) = parent else {
                for &other in rows {
                    if mapping[other] != usize::MAX {
                        tokens.hydrogen_known[mapping[other]] = 0;
                    }
                }
                continue;
            };
            let token = mapping[parent];
            let assigned = assigned.entry(parent).or_default();
            assigned.push(row);
            let count = tokens.masks[token * 32 + 15..token * 32 + 18]
                .iter()
                .filter(|&&value| value != 0)
                .count();
            if count >= 3 {
                return Err(hydrogen_overflow(source, features, parent, assigned));
            }
            let slot = 15 + count;
            tokens.masks[token * 32 + slot] = 1;
            tokens.positions[token * 96 + slot * 3..][..3]
                .copy_from_slice(&source.atoms[row].position);
        }
    }
    for (row, atom) in source.atoms.iter().enumerate() {
        let token = mapping[row];
        if token == usize::MAX || tokens.ligand[token] == 0 {
            continue;
        }
        let feature = ligand_features(features, atom);
        let count = tokens.masks[token * 32 + 15..(token + 1) * 32]
            .iter()
            .filter(|&&value| value != 0)
            .count();
        tokens.charge_known[token] = u8::from(
            tokens.hydrogen_known[token] != 0 && feature[3] != 0 && count == feature[2] as usize,
        );
    }
    Ok(())
}

pub(super) fn pack(
    source: &Source,
    schema: &Schema,
    features: &LigandFeatures,
    parents: &HydrogenParents,
    elements: &[String],
    excluded: &[String],
    options: &super::PackOptions<'_>,
) -> Result<Tokens, String> {
    let super::PackOptions {
        complete,
        topology_policy: policy,
        protein_only,
        ..
    } = *options;
    if protein_only && !source.has_protein_annotations {
        return Err("Protein export requires entity_poly protein annotations".into());
    }
    let mut result = Tokens::empty();
    let mut mapping = vec![usize::MAX; source.atoms.len()];
    add_proteins(
        source,
        schema,
        complete,
        protein_only,
        policy,
        &mut result,
        &mut mapping,
    )?;
    if !protein_only {
        add_ligands(
            source,
            features,
            parents,
            elements,
            excluded,
            &mut result,
            &mut mapping,
            complete,
        )?;
    }
    let chains: Vec<_> = if protein_only {
        let mut value = Vec::new();
        for key in &result.keys {
            if !value.contains(&key[0]) {
                value.push(key[0].clone());
            }
        }
        value
    } else {
        result
            .keys
            .iter()
            .map(|key| key[0].clone())
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect()
    };
    for (row, key) in result.keys.iter().enumerate() {
        result.chains[row] = chains
            .iter()
            .position(|chain| chain == &key[0])
            .ok_or("Missing token chain")? as i64;
    }
    if !protein_only {
        let mut pairs = BTreeSet::new();
        for &(a, b) in &source.bonds {
            let (left, right) = (mapping[a], mapping[b]);
            if left == usize::MAX
                || right == usize::MAX
                || source.atoms[a].protein == source.atoms[b].protein
            {
                continue;
            }
            pairs.insert(if source.atoms[a].protein {
                (right, left)
            } else {
                (left, right)
            });
        }
        result.attachments = pairs.into_iter().collect();
        for &(ligand, _) in &result.attachments {
            result.covalent[ligand] = 1;
        }
        renumber_ligands(&mut result);
    }
    if result.keys.iter().collect::<BTreeSet<_>>().len() != result.keys.len() {
        return Err("PDB tokens have duplicate mmCIF identities".into());
    }
    if protein_only && result.identity.is_empty() {
        return Err("No canonical protein residues with CA coordinates".into());
    }
    Ok(result)
}
