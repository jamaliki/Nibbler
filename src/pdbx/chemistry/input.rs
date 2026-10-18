use super::super::{ComponentRegistry, ShortText};
use super::super::model::{AtomSite, ComponentDefinition, PdbxModel};
use super::{PreparedChemistry, ProductAtom, SemanticError};
use crate::cif::FastMap;
use std::collections::{BTreeMap, BTreeSet};
use std::sync::Arc;

#[derive(Clone, Debug)]
pub(super) struct Residue {
    pub atoms: NameMap,
    pub aliases: NameMap,
    pub component_type: Option<ShortText>,
    pub selected_alt: Option<String>,
    pub entity_id: ShortText,
    pub auth_chain: ShortText,
    pub auth_component: ShortText,
}

/// Names in byte order, each with a position in a component definition.
#[derive(Debug, Default)]
pub(super) struct SortedNames(Vec<(String, usize)>);

impl SortedNames {
    fn position(&self, name: &str) -> Option<usize> {
        self.0
            .binary_search_by(|(key, _)| key.as_str().cmp(name))
            .ok()
            .map(|index| self.0[index].1)
    }
}

/// One residue's map from names to product atom indices.
///
/// Every residue of a component shares the template names, whose positions are offset by
/// `base`. `own` holds other names; an own name hides an equal template name.
#[derive(Clone, Debug, Default)]
pub(super) struct NameMap {
    shared: Option<Arc<SortedNames>>,
    base: usize,
    own: BTreeMap<String, usize>,
}

impl NameMap {
    fn shared(names: &Arc<SortedNames>, base: usize) -> Self {
        Self {
            shared: Some(Arc::clone(names)),
            base,
            own: BTreeMap::new(),
        }
    }

    pub(super) fn get(&self, name: &str) -> Option<usize> {
        if let Some(&index) = self.own.get(name) {
            return Some(index);
        }
        Some(self.base + self.shared.as_ref()?.position(name)?)
    }

    pub(super) fn contains_key(&self, name: &str) -> bool {
        self.get(name).is_some()
    }

    /// Add or replace a name; true when the name was already present.
    fn insert(&mut self, name: &str, index: usize) -> bool {
        let present = self.contains_key(name);
        self.own.insert(name.to_owned(), index);
        present
    }

    /// The index under the first name in byte order.
    pub(super) fn first(&self) -> Option<usize> {
        let shared = self.shared.as_ref().and_then(|names| names.0.first());
        match (self.own.iter().next(), shared) {
            (Some((own, &index)), Some((name, _))) if own <= name => Some(index),
            (_, Some(&(_, position))) => Some(self.base + position),
            (Some((_, &index)), None) => Some(index),
            (None, None) => None,
        }
    }

    /// Indices in name order.
    fn values(&self) -> Vec<usize> {
        let shared = self.shared.as_ref().map_or(&[][..], |names| &names.0[..]);
        let mut output = Vec::with_capacity(shared.len() + self.own.len());
        let mut own = self.own.iter().peekable();
        for (name, position) in shared {
            while let Some((_, &index)) = own.next_if(|(key, _)| *key < name) {
                output.push(index);
            }
            match own.next_if(|(key, _)| *key == name) {
                Some((_, &index)) => output.push(index),
                None => output.push(self.base + position),
            }
        }
        output.extend(own.map(|(_, &index)| index));
        output
    }
}

type ResidueKey<'a> = (&'a str, Option<i64>, &'a str, Option<&'a str>);

type Site<'a> = (usize, &'a AtomSite);

/// Observed atoms of one residue, sorted by atom name.
type Observed<'a> = Vec<(&'a str, Site<'a>)>;

/// Buffers that every residue reuses.
#[derive(Default)]
struct Scratch<'a> {
    selected: Vec<Site<'a>>,
    observed: Observed<'a>,
    sites: Vec<Option<Site<'a>>>,
    unlisted: Observed<'a>,
}

fn fail(message: impl Into<String>) -> SemanticError {
    SemanticError::new("PDBX_LINK_PRODUCT", message, vec!["atom_site".into()])
}

fn is_hydrogen_element(element: &str) -> bool {
    ["H", "D", "T"]
        .iter()
        .any(|name| element.eq_ignore_ascii_case(name))
}

/// Facts about one component definition that every residue of that component shares.
struct Template<'d> {
    /// Atom names in byte order with their definition positions.
    names: Arc<SortedNames>,
    /// Whether every definition atom has a distinct name.
    unique: bool,
    /// `(alternate, canonical)` names, in definition order, for alternates outside `names`.
    renames: Vec<(&'d str, &'d str)>,
    /// Bond endpoints as definition positions, with the parsed order; `None` for a bond
    /// with a missing endpoint and an order of `None` when it is unsupported.
    bonds: Vec<Option<(usize, usize, Option<u8>)>>,
    /// `(alternate, position)` for alternates outside `names`, in definition order.
    aliases: Vec<(&'d str, usize)>,
    /// `aliases` in byte order, for residues whose names are exactly the template names.
    alias_names: Arc<SortedNames>,
    /// Whether two entries of `aliases` share an alternate name.
    repeated_alias: bool,
}

impl<'d> Template<'d> {
    fn new(definition: &'d ComponentDefinition) -> Self {
        let atoms = &definition.atoms;
        let mut by_name: Vec<(&str, usize)> = atoms
            .iter()
            .enumerate()
            .map(|(position, atom)| (atom.atom_id.as_str(), position))
            .collect();
        by_name.sort_by(|left, right| left.0.cmp(right.0));
        let unique = by_name.windows(2).all(|pair| pair[0].0 != pair[1].0);
        let outside = |alias: &&str| {
            by_name
                .binary_search_by(|&(key, _)| key.cmp(alias))
                .is_err()
        };
        let renames = atoms
            .iter()
            .filter_map(|atom| {
                let alias = atom.alternate_id.as_deref().filter(outside)?;
                Some((alias, atom.atom_id.as_str()))
            })
            .collect();
        let aliases: Vec<(&str, usize)> = atoms
            .iter()
            .enumerate()
            .filter_map(|(position, atom)| {
                Some((atom.alternate_id.as_deref().filter(outside)?, position))
            })
            .collect();
        let mut alias_names = aliases.clone();
        alias_names.sort_by(|left, right| left.0.cmp(right.0));
        let repeated_alias = alias_names.windows(2).any(|pair| pair[0].0 == pair[1].0);
        let position = |name: &str| {
            by_name
                .binary_search_by(|&(key, _)| key.cmp(name))
                .ok()
                .map(|index| by_name[index].1)
        };
        let bonds = definition
            .bonds
            .iter()
            .map(|bond| {
                let first = position(&bond.first_atom_id)?;
                let second = position(&bond.second_atom_id)?;
                Some((first, second, order(&bond.order).ok()))
            })
            .collect();
        let owned = |names: &[(&str, usize)]| {
            Arc::new(SortedNames(
                names
                    .iter()
                    .map(|&(name, position)| (name.to_owned(), position))
                    .collect(),
            ))
        };
        Self {
            names: owned(&by_name),
            unique,
            renames,
            bonds,
            aliases,
            alias_names: owned(&alias_names),
            repeated_alias,
        }
    }
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
    // Atoms of one residue are usually adjacent, so a run of equal keys needs one lookup.
    // Sorting the groups afterwards gives the residue order a sorted map would.
    let mut lookup: FastMap<ResidueKey<'_>, usize> = FastMap::default();
    let mut groups: Vec<(ResidueKey<'_>, Vec<(usize, &AtomSite)>)> = Vec::new();
    let mut previous: Option<(ResidueKey<'_>, usize)> = None;
    for (row, atom) in model.atom_sites.iter().enumerate().filter(|(_, a)| {
        a.model_number == model_number && (use_observed_hydrogens || !is_hydrogen_element(&a.element))
    }) {
        let key = (
            atom.label_asym_id.as_str(),
            atom.label_seq_id,
            atom.auth_seq_id.as_str(),
            atom.insertion_code.as_deref(),
        );
        let slot = match previous {
            Some((last, slot)) if last == key => slot,
            _ => *lookup.entry(key).or_insert_with(|| {
                groups.push((key, Vec::new()));
                groups.len() - 1
            }),
        };
        previous = Some((key, slot));
        groups[slot].1.push((row, atom));
    }
    if groups.is_empty() {
        return Err(fail(format!("Deposited model {model_number} has no atoms")));
    }
    groups.sort_unstable_by(|left, right| left.0.cmp(&right.0));
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
    let templates: BTreeMap<&str, Template<'_>> = prepared_definitions
        .iter()
        .filter(|(_, definition)| !definition.atoms.is_empty())
        .map(|(&name, definition)| (name, Template::new(definition)))
        .collect();
    // Estimate the product size so the large vectors are not copied as they grow. Links
    // add at most one bond per residue plus one per authored connection.
    let (mut atom_count, mut bond_count) = (0, groups.len() + 256);
    for (_, rows) in &groups {
        let component = rows.first().map(|(_, a)| a.label_comp_id.as_str());
        match component.and_then(|component| prepared_definitions.get(component)) {
            Some(definition) if !definition.atoms.is_empty() => {
                atom_count += definition.atoms.len();
                bond_count += definition.bonds.len();
            }
            _ => atom_count += rows.len(),
        }
    }
    let mut graph = PreparedChemistry {
        schema_version: "linked-product.v1",
        atoms: Vec::with_capacity(atom_count),
        bonds: Vec::with_capacity(bond_count),
        edits: Vec::new(),
        diagnostics: Vec::new(),
        adjacency: Vec::with_capacity(atom_count),
        residues: Vec::with_capacity(groups.len()),
    };
    let mut scratch = Scratch::default();
    for (_, rows) in &groups {
        let selected_alt = select_alternate(rows, alternate);
        let mut selected = std::mem::take(&mut scratch.selected);
        selected.clear();
        selected.extend(rows.iter().copied().filter(|(_, a)| {
            a.label_alt_id.is_none() || a.label_alt_id.as_deref() == selected_alt.as_deref()
        }));
        let Some((_, first)) = selected.first() else {
            scratch.selected = selected;
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
        let component = first.label_comp_id.as_str();
        let definition = prepared_definitions.get(component);
        let template = templates.get(component);
        add_residue(&mut graph, &selected, definition, template, selected_alt, &mut scratch)?;
        scratch.selected = selected;
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

fn add_residue<'a>(
    graph: &mut PreparedChemistry,
    rows: &[Site<'a>],
    definition: Option<&ComponentDefinition>,
    template: Option<&Template<'a>>,
    selected_alt: Option<String>,
    scratch: &mut Scratch<'a>,
) -> Result<(), SemanticError> {
    let Some((_, first)) = rows.first() else {
        return Ok(());
    };
    let residue = graph.residues.len();
    let observed = &mut scratch.observed;
    observed.clear();
    observed.extend(rows.iter().map(|&(r, a)| (a.label_atom_id.as_str(), (r, a))));
    observed.sort_by(|left, right| left.0.cmp(right.0));
    if observed.windows(2).any(|pair| pair[0].0 == pair[1].0) {
        return Err(fail("Duplicate atom identity in the selected conformer"));
    }
    let find = |observed: &Observed<'a>, name: &str| {
        observed.binary_search_by(|&(key, _)| key.cmp(name))
    };
    let mut indices;
    let mut aliases = NameMap::default();
    if let (Some(comp), Some(template)) = (definition, template) {
        for &(alias, name) in &template.renames {
            if let Ok(index) = find(observed, alias) {
                let (_, site) = observed.remove(index);
                match find(observed, name) {
                    Ok(_) => {
                        return Err(fail(
                            "Canonical and alternate CCD names select the same atom",
                        ));
                    }
                    Err(index) => observed.insert(index, (name, site)),
                }
            }
        }
        if !template.unique {
            return Err(fail("Duplicate atom in component definition"));
        }
        // Both lists are sorted by name, so one pass pairs each definition atom with its
        // observed site and keeps the observed atoms the definition does not name.
        let sites = &mut scratch.sites;
        sites.clear();
        sites.resize(comp.atoms.len(), None);
        let unlisted = &mut scratch.unlisted;
        unlisted.clear();
        let names = &template.names.0;
        let mut next = 0;
        for &(name, site) in observed.iter() {
            while next < names.len() && names[next].0.as_str() < name {
                next += 1;
            }
            match names.get(next) {
                Some((key, position)) if key == name => {
                    sites[*position] = Some(site);
                    next += 1;
                }
                _ => unlisted.push((name, site)),
            }
        }
        let base = graph.atoms.len();
        for (a, &site) in comp.atoms.iter().zip(sites.iter()) {
            if site.is_some_and(|(_, s)| {
                !s.element.eq_ignore_ascii_case(&a.element)
                    && !(a.element.eq_ignore_ascii_case("H")
                        && (s.element.eq_ignore_ascii_case("D")
                            || s.element.eq_ignore_ascii_case("T")))
            }) {
                return Err(fail(format!(
                    "Element mismatch for {}/{}",
                    comp.id, a.atom_id
                )));
            }
            graph.atoms.push(ProductAtom {
                model_number: first.model_number,
                asym_id: first.label_asym_id.clone(),
                sequence_id: first.label_seq_id,
                auth_sequence_id: first.auth_seq_id.clone(),
                component_id: first.label_comp_id.clone(),
                insertion_code: first.insertion_code.clone(),
                atom_id: ShortText::new(&a.atom_id),
                alternate_id: site
                    .and_then(|(_, s)| s.label_alt_id.clone())
                    .or_else(|| selected_alt.as_deref().map(ShortText::new)),
                template_atom_id: Some(ShortText::new(&a.atom_id)),
                deposited_charge: site.and_then(|(_, s)| s.formal_charge),
                source_row: site.map(|(r, _)| r),
                coordinates: site.map(|(_, s)| [s.x, s.y, s.z]),
                element: site.map_or_else(
                    || ShortText::ascii_uppercase(&a.element),
                    |(_, s)| ShortText::ascii_uppercase(&s.element),
                ),
                formal_charge: a.formal_charge,
                template_charge: a.formal_charge,
                hydrogen_count: Some(0),
                active: true,
                residue,
                leaving: a.leaving == Some(true),
                impacted: false,
            });
            graph.adjacency.push(super::BondList::default());
        }
        indices = NameMap::shared(&template.names, base);
        for (bond, resolved) in comp.bonds.iter().zip(&template.bonds) {
            let Some((a, b, known)) = *resolved else {
                return Err(fail(format!(
                    "Component {} bond has a missing endpoint",
                    comp.id
                )));
            };
            let bond_order = match known {
                Some(value) => value,
                None => order(&bond.order)?,
            };
            graph.add_bond(base + a, base + b, bond_order, super::BondSource::Component)?;
        }
        map_unlisted_hydrogens(graph, &mut indices, unlisted, &comp.id)?;
        // Each canonical name keeps its own atom and an alternate that is itself a
        // canonical name is always present, so only the precomputed alternates can apply.
        // Without mapped hydrogens all of them apply.
        if indices.own.is_empty() && !template.repeated_alias {
            aliases = NameMap::shared(&template.alias_names, base);
        } else {
            for &(alias, position) in &template.aliases {
                if !indices.contains_key(alias) && aliases.insert(alias, base + position) {
                    return Err(fail("CCD alternate atom name has multiple targets"));
                }
            }
        }
    } else {
        indices = NameMap::default();
        for &(row, a) in rows {
            let index = graph.atoms.len();
            indices.insert(&a.label_atom_id, index);
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
                element: ShortText::ascii_uppercase(&a.element),
                formal_charge: a.formal_charge,
                template_charge: None,
                hydrogen_count: None,
                active: true,
                residue,
                leaving: false,
                impacted: false,
            });
            graph.adjacency.push(super::BondList::default());
        }
    }
    graph.residues.push(Residue {
        atoms: indices,
        aliases,
        component_type: definition.and_then(|c| c.component_type.as_deref().map(ShortText::new)),
        selected_alt,
        entity_id: first.label_entity_id.clone(),
        auth_chain: first.auth_asym_id.clone(),
        auth_component: first.auth_comp_id.clone(),
    });
    Ok(())
}

/// Give observed hydrogens that the definition does not name to template hydrogens.
fn map_unlisted_hydrogens(
    graph: &mut PreparedChemistry,
    indices: &mut NameMap,
    unlisted: &Observed<'_>,
    component: &str,
) -> Result<(), SemanticError> {
    for &(name, (row, site)) in unlisted {
        if !is_hydrogen_element(&site.element) {
            return Err(fail(format!(
                "Observed heavy atom is absent from component {component}"
            )));
        }
        let point = [site.x, site.y, site.z];
        let mut parents: Vec<_> = indices
            .values()
            .into_iter()
            .filter_map(|a| {
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
        graph.atoms[h].element = ShortText::ascii_uppercase(&site.element);
        graph.atoms[h].deposited_charge = site.formal_charge;
        indices.insert(name, h);
        graph.edits.push(super::ChemicalEdit {
            atom_indices: vec![h],
            reason: "hydrogen_name_mapping",
        });
    }
    Ok(())
}

pub(super) fn order(value: &str) -> Result<u8, SemanticError> {
    [
        ("sing", 1),
        ("single", 1),
        ("doub", 2),
        ("double", 2),
        ("trip", 3),
        ("triple", 3),
    ]
    .into_iter()
    .find(|(name, _)| value.eq_ignore_ascii_case(name))
    .map(|(_, order)| order)
    .ok_or_else(|| fail(format!("Unsupported bond order {value:?}")))
}
