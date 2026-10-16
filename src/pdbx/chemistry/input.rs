use super::super::ComponentRegistry;
use super::super::model::{AtomSite, ComponentDefinition, PdbxModel};
use super::{PreparedChemistry, ProductAtom, SemanticError};
use std::collections::{BTreeMap, BTreeSet};

#[derive(Clone, Debug)]
pub(super) struct Residue {
    pub atoms: BTreeMap<String, usize>,
    pub aliases: BTreeMap<String, usize>,
    pub component_type: Option<String>,
    pub selected_alt: Option<String>,
    pub entity_id: String,
    pub auth_chain: String,
    pub auth_component: String,
}

type ResidueKey = (String, Option<i64>, String, Option<String>);

fn fail(message: impl Into<String>) -> SemanticError {
    SemanticError::new("PDBX_LINK_PRODUCT", message, vec!["atom_site".into()])
}

pub(super) fn build(
    model: &PdbxModel,
    registry: Option<&ComponentRegistry>,
    model_number: i64,
    alternate: Option<&str>,
    use_observed_hydrogens: bool,
) -> Result<PreparedChemistry, SemanticError> {
    if model_number < 1 {
        return Err(fail("Model number must be positive"));
    }
    let mut groups: BTreeMap<ResidueKey, Vec<(usize, &AtomSite)>> = BTreeMap::new();
    for (row, atom) in model.atom_sites.iter().enumerate().filter(|(_, a)| {
        a.model_number == model_number
            && (use_observed_hydrogens
                || !matches!(a.element.to_ascii_uppercase().as_str(), "H" | "D" | "T"))
    }) {
        groups
            .entry((
                atom.label_asym_id.clone(),
                atom.label_seq_id,
                atom.auth_seq_id.clone(),
                atom.insertion_code.clone(),
            ))
            .or_default()
            .push((row, atom));
    }
    if groups.is_empty() {
        return Err(fail(format!("Deposited model {model_number} has no atoms")));
    }
    let definitions: BTreeMap<_, _> = model
        .components
        .iter()
        .map(|c| (c.id.as_str(), c))
        .collect();
    let mut prepared_definitions = BTreeMap::new();
    for (&name, &embedded) in &definitions {
        let definition = supplement(embedded, registry.and_then(|r| r.get(name)))?;
        prepared_definitions.insert(name, definition);
    }
    let mut graph = PreparedChemistry {
        schema_version: "linked-product.v1",
        atoms: Vec::new(),
        bonds: Vec::new(),
        edits: Vec::new(),
        diagnostics: Vec::new(),
        adjacency: Vec::new(),
        residues: Vec::new(),
    };
    for rows in groups.values() {
        let selected_alt = select_alternate(rows, alternate);
        let selected: Vec<_> = rows
            .iter()
            .copied()
            .filter(|(_, a)| a.label_alt_id.is_none() || a.label_alt_id == selected_alt)
            .collect();
        let Some((_, first)) = selected.first() else {
            continue;
        };
        if selected
            .iter()
            .any(|(_, a)| a.label_comp_id != first.label_comp_id)
        {
            return Err(fail(
                "Selected conformer has conflicting component identities",
            ));
        }
        let definition = prepared_definitions.get(first.label_comp_id.as_str());
        add_residue(&mut graph, &selected, definition, selected_alt)?;
    }
    Ok(graph)
}

fn supplement(
    embedded: &ComponentDefinition,
    cached: Option<&ComponentDefinition>,
) -> Result<ComponentDefinition, SemanticError> {
    let Some(cached) = cached else {
        return Ok(embedded.clone());
    };
    if embedded.atoms.is_empty() {
        return Ok(cached.clone());
    }
    let mut result = embedded.clone();
    result.component_type = result.component_type.or(cached.component_type.clone());
    let original: BTreeSet<_> = embedded.atoms.iter().map(|a| a.atom_id.as_str()).collect();
    for a in &mut result.atoms {
        if let Some(c) = cached.atoms.iter().find(|c| c.atom_id == a.atom_id) {
            if !c.element.eq_ignore_ascii_case(&a.element) {
                return Err(fail("Component and CCD element definitions conflict"));
            }
            a.formal_charge = a.formal_charge.or(c.formal_charge);
            a.leaving = a.leaving.or(c.leaving);
            a.alternate_id = a.alternate_id.clone().or(c.alternate_id.clone());
        }
    }
    result.atoms.extend(
        cached
            .atoms
            .iter()
            .filter(|a| !original.contains(a.atom_id.as_str()))
            .cloned(),
    );
    if result.bonds.is_empty() {
        result.bonds = cached.bonds.clone();
    } else {
        result.bonds.extend(
            cached
                .bonds
                .iter()
                .filter(|b| {
                    !original.contains(b.first_atom_id.as_str())
                        || !original.contains(b.second_atom_id.as_str())
                })
                .cloned(),
        );
    }
    Ok(result)
}

fn select_alternate(rows: &[(usize, &AtomSite)], requested: Option<&str>) -> Option<String> {
    let mut scores: BTreeMap<&str, f64> = BTreeMap::new();
    for (_, a) in rows {
        if let Some(alt) = a.label_alt_id.as_deref() {
            *scores.entry(alt).or_default() += a.occupancy;
        }
    }
    if let Some(alt) = requested.filter(|alt| scores.contains_key(alt)) {
        return Some(alt.into());
    }
    scores
        .into_iter()
        .max_by(|(a, x), (b, y)| x.total_cmp(y).then_with(|| b.cmp(a)))
        .map(|(a, _)| a.into())
}

fn add_residue(
    graph: &mut PreparedChemistry,
    rows: &[(usize, &AtomSite)],
    definition: Option<&ComponentDefinition>,
    selected_alt: Option<String>,
) -> Result<(), SemanticError> {
    let Some((_, first)) = rows.first() else {
        return Ok(());
    };
    let residue = graph.residues.len();
    let mut indices = BTreeMap::new();
    let mut observed: BTreeMap<_, _> = rows
        .iter()
        .map(|&(r, a)| (a.label_atom_id.as_str(), (r, a)))
        .collect();
    if observed.len() != rows.len() {
        return Err(fail("Duplicate atom identity in the selected conformer"));
    }
    if let Some(comp) = definition.filter(|c| !c.atoms.is_empty()) {
        let names: BTreeSet<_> = comp.atoms.iter().map(|a| a.atom_id.as_str()).collect();
        for atom in &comp.atoms {
            if let Some(alias) = atom
                .alternate_id
                .as_deref()
                .filter(|&alias| !names.contains(alias))
                && let Some(site) = observed.remove(alias)
            {
                if observed.insert(atom.atom_id.as_str(), site).is_some() {
                    return Err(fail(
                        "Canonical and alternate CCD names select the same atom",
                    ));
                }
            }
        }
        if names.len() != comp.atoms.len() {
            return Err(fail("Duplicate atom in component definition"));
        }
        for a in &comp.atoms {
            let site = observed.get(a.atom_id.as_str()).copied();
            if site.is_some_and(|(_, s)| {
                !s.element.eq_ignore_ascii_case(&a.element)
                    && !(a.element.eq_ignore_ascii_case("H")
                        && matches!(s.element.to_ascii_uppercase().as_str(), "D" | "T"))
            }) {
                return Err(fail(format!(
                    "Element mismatch for {}/{}",
                    comp.id, a.atom_id
                )));
            }
            let index = graph.atoms.len();
            indices.insert(a.atom_id.clone(), index);
            graph.atoms.push(ProductAtom {
                model_number: first.model_number,
                asym_id: first.label_asym_id.clone(),
                sequence_id: first.label_seq_id,
                auth_sequence_id: first.auth_seq_id.clone(),
                component_id: first.label_comp_id.clone(),
                insertion_code: first.insertion_code.clone(),
                atom_id: a.atom_id.clone(),
                alternate_id: site
                    .and_then(|(_, s)| s.label_alt_id.clone())
                    .or(selected_alt.clone()),
                template_atom_id: Some(a.atom_id.clone()),
                deposited_charge: site.and_then(|(_, s)| s.formal_charge),
                source_row: site.map(|(r, _)| r),
                coordinates: site.map(|(_, s)| [s.x, s.y, s.z]),
                element: site.map_or_else(
                    || a.element.to_ascii_uppercase(),
                    |(_, s)| s.element.to_ascii_uppercase(),
                ),
                formal_charge: a.formal_charge,
                template_charge: a.formal_charge,
                hydrogen_count: Some(0),
                active: true,
                residue,
                leaving: a.leaving == Some(true),
                impacted: false,
            });
            graph.adjacency.push(Vec::new());
        }
        for bond in &comp.bonds {
            let (Some(&a), Some(&b)) = (
                indices.get(&bond.first_atom_id),
                indices.get(&bond.second_atom_id),
            ) else {
                return Err(fail(format!(
                    "Component {} bond has a missing endpoint",
                    comp.id
                )));
            };
            graph.add_bond(a, b, order(&bond.order)?, super::BondSource::Component)?;
        }
        map_unlisted_hydrogens(graph, &mut indices, &observed, &names, &comp.id)?;
    } else {
        for &(row, a) in rows {
            let index = graph.atoms.len();
            indices.insert(a.label_atom_id.clone(), index);
            graph.atoms.push(ProductAtom {
                model_number: a.model_number,
                asym_id: a.label_asym_id.clone(),
                sequence_id: a.label_seq_id,
                auth_sequence_id: a.auth_seq_id.clone(),
                component_id: a.label_comp_id.clone(),
                insertion_code: a.insertion_code.clone(),
                atom_id: a.label_atom_id.clone(),
                alternate_id: a.label_alt_id.clone(),
                source_row: Some(row),
                template_atom_id: None,
                deposited_charge: a.formal_charge,
                coordinates: Some([a.x, a.y, a.z]),
                element: a.element.to_ascii_uppercase(),
                formal_charge: a.formal_charge,
                template_charge: None,
                hydrogen_count: None,
                active: true,
                residue,
                leaving: false,
                impacted: false,
            });
            graph.adjacency.push(Vec::new());
        }
    }
    let mut aliases = BTreeMap::new();
    if let Some(comp) = definition {
        for atom in &comp.atoms {
            if let Some(alias) = atom
                .alternate_id
                .as_ref()
                .filter(|alias| !indices.contains_key(*alias))
                && let Some(&index) = indices.get(&atom.atom_id)
            {
                if aliases.insert(alias.clone(), index).is_some() {
                    return Err(fail("CCD alternate atom name has multiple targets"));
                }
            }
        }
    }
    graph.residues.push(Residue {
        atoms: indices,
        aliases,
        component_type: definition.and_then(|c| c.component_type.clone()),
        selected_alt,
        entity_id: first.label_entity_id.clone(),
        auth_chain: first.auth_asym_id.clone(),
        auth_component: first.auth_comp_id.clone(),
    });
    Ok(())
}

fn map_unlisted_hydrogens(
    graph: &mut PreparedChemistry,
    indices: &mut BTreeMap<String, usize>,
    observed: &BTreeMap<&str, (usize, &AtomSite)>,
    names: &BTreeSet<&str>,
    component: &str,
) -> Result<(), SemanticError> {
    for (&name, &(row, site)) in observed.iter().filter(|(name, _)| !names.contains(**name)) {
        if !matches!(site.element.to_ascii_uppercase().as_str(), "H" | "D" | "T") {
            return Err(fail(format!(
                "Observed heavy atom is absent from component {component}"
            )));
        }
        let point = [site.x, site.y, site.z];
        let mut parents: Vec<_> = indices
            .values()
            .filter_map(|&a| {
                let atom = &graph.atoms[a];
                if atom.is_hydrogen() {
                    return None;
                }
                atom.coordinates.map(|p| {
                    (
                        p.iter()
                            .zip(point)
                            .map(|(x, y)| (x - y).powi(2))
                            .sum::<f64>()
                            .sqrt(),
                        a,
                    )
                })
            })
            .collect();
        parents.sort_by(|a, b| a.0.total_cmp(&b.0));
        let Some(&(distance, parent)) = parents.first() else {
            return Err(fail("Unlisted hydrogen has no observed parent"));
        };
        if !(0.5..1.6).contains(&distance) || parents.get(1).is_some_and(|p| p.0 - distance < 0.2) {
            return Err(fail("Unlisted hydrogen has an ambiguous parent"));
        }
        let Some(h) = graph
            .hydrogens(parent)
            .into_iter()
            .find(|&h| graph.atoms[h].source_row.is_none())
        else {
            return Err(fail(
                "Unlisted hydrogen exceeds the component hydrogen count",
            ));
        };
        graph.atoms[h].atom_id = name.into();
        graph.atoms[h].source_row = Some(row);
        graph.atoms[h].coordinates = Some(point);
        graph.atoms[h].alternate_id = site.label_alt_id.clone();
        graph.atoms[h].element = site.element.to_ascii_uppercase();
        graph.atoms[h].deposited_charge = site.formal_charge;
        indices.insert(name.into(), h);
        graph.edits.push(super::ChemicalEdit {
            atom_indices: vec![h],
            reason: "hydrogen_name_mapping",
        });
    }
    Ok(())
}

pub(super) fn order(value: &str) -> Result<u8, SemanticError> {
    match value.to_ascii_lowercase().as_str() {
        "sing" | "single" => Ok(1),
        "doub" | "double" => Ok(2),
        "trip" | "triple" => Ok(3),
        _ => Err(fail(format!("Unsupported bond order {value:?}"))),
    }
}
