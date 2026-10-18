//! Dictionary-backed metal hydrides that cannot use a single-parent H slot.

use super::source::Source;
use super::residue;
use crate::cif::{self, BlockCategory};
use reduce3::monlib::{Comp, MonLib};
use std::collections::{BTreeSet, HashMap};
use std::path::PathBuf;
use std::sync::{Arc, Mutex, OnceLock, PoisonError};

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
#[repr(u8)]
pub(crate) enum EdgeKind {
    MetalHydride = 1,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) struct Definition {
    pub metals: [String; 2],
    pub metal_elements: [String; 2],
    pub formal_charge: i16,
    pub kind: EdgeKind,
}

pub(crate) type Definitions = HashMap<String, Definition>;
type CacheKey = (PathBuf, PathBuf, String, Vec<String>, bool);
static DEFINITIONS: OnceLock<Mutex<HashMap<CacheKey, Arc<Definitions>>>> = OnceLock::new();

#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) struct Site {
    pub group: usize,
    pub hydrogen_name: String,
    pub hydrogen_row: Option<usize>,
    pub metal_rows: [usize; 2],
    pub formal_charge: i16,
    pub kind: EdgeKind,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) struct Edge {
    pub hydrogen: usize,
    pub metal: usize,
    pub kind: EdgeKind,
}

impl Site {
    /// Map both dependencies after the H token and its metal tokens exist.
    pub(crate) fn edges(&self, hydrogen: usize, mapping: &[usize]) -> Result<[Edge; 2], String> {
        let mut result = [Edge { hydrogen, metal: usize::MAX, kind: self.kind }; 2];
        for (edge, &row) in result.iter_mut().zip(&self.metal_rows) {
            edge.metal = *mapping.get(row).ok_or("Bridge metal row is outside the source")?;
            if edge.metal == usize::MAX || edge.metal == hydrogen {
                return Err("Bridge hydrogen requires both metal tokens".into());
            }
        }
        if result[0].metal == result[1].metal {
            return Err("Bridge hydrogen requires two distinct metal tokens".into());
        }
        Ok(result)
    }
}

fn hydrogen(element: &str) -> bool {
    ["H", "D", "T"].iter().any(|value| element.eq_ignore_ascii_case(value))
}

fn observed_name(name: &str, atom: &super::source::Atom) -> Option<bool> {
    if atom.name == name { return Some(hydrogen(&atom.element)); }
    for isotope in ["D", "T"] {
        if name.contains('H') && atom.name == name.replacen('H', isotope, 1) {
            return Some(atom.element.eq_ignore_ascii_case(isotope));
        }
    }
    None
}

fn metal(element: &str) -> bool {
    ["FE", "NI"].iter().any(|value| element.eq_ignore_ascii_case(value))
}

pub(super) fn possible_source(source: &Source) -> bool {
    source.groups.iter().any(|rows| rows.iter().filter(|&&row| {
        let atom = &source.atoms[row];
        !atom.protein && metal(&atom.element)
    }).take(2).count() == 2)
}

fn candidates(comp: &Comp) -> Definitions {
    let mut result = HashMap::new();
    for atom in comp.atoms.iter().filter(|atom| hydrogen(&atom.type_symbol)) {
        let mut neighbors = BTreeSet::new();
        let mut invalid = false;
        for bond in &comp.bonds {
            let other = if bond.a1 == atom.id { &bond.a2 }
                else if bond.a2 == atom.id { &bond.a1 } else { continue };
            let Some(other_atom) = comp.atom(other) else { invalid = true; break };
            if !metal(&other_atom.type_symbol)
                || !["single", "sing", "coval"].iter().any(|value| bond.type_.eq_ignore_ascii_case(value)) {
                invalid = true;
                break;
            }
            neighbors.insert((other.clone(), other_atom.type_symbol.to_ascii_uppercase()));
        }
        let values: Vec<_> = neighbors.into_iter().collect();
        if !invalid && let [(first, first_element), (second, second_element)] = &values[..] {
            result.insert(atom.id.clone(), Definition {
                metals: [first.clone(), second.clone()],
                metal_elements: [first_element.clone(), second_element.clone()],
                formal_charge: -1, kind: EdgeKind::MetalHydride,
            });
        }
    }
    result
}

fn text(table: &BlockCategory<'_>, row: usize, name: &str) -> Option<String> {
    table.column(name).and_then(|column| table.text(row, column)).map(|value| value.into_owned())
}

fn source_definition(comp: &Comp, name: &str, candidate: &Definition, document: &cif::CifDocument) -> bool {
    let Some(block) = document.blocks().iter().find(|block| {
        BlockCategory::find(block.entries(), "chem_comp_atom").is_some_and(|atoms| {
            (0..atoms.row_count()).any(|row| text(&atoms, row, "comp_id").as_deref() == Some(comp.id.as_str()))
        })
    }) else { return false };
    let Some(atoms) = BlockCategory::find(block.entries(), "chem_comp_atom") else { return false };
    let charged: Vec<_> = (0..atoms.row_count()).filter(|&row| {
        text(&atoms, row, "comp_id").as_deref() == Some(comp.id.as_str())
            && text(&atoms, row, "atom_id").as_deref() == Some(name)
    }).collect();
    let [row] = charged[..] else { return false };
    if text(&atoms, row, "charge").as_deref().and_then(cif::numeric::parse_integer) != Some(-1)
        || !text(&atoms, row, "type_symbol").is_some_and(|value| hydrogen(&value)) { return false; }
    let Some(bonds) = BlockCategory::find(block.entries(), "chem_comp_bond") else { return false };
    let mut neighbors = BTreeSet::new();
    for row in 0..bonds.row_count() {
        if text(&bonds, row, "comp_id").as_deref() != Some(comp.id.as_str()) { continue; }
        let (Some(first), Some(second)) = (text(&bonds, row, "atom_id_1"), text(&bonds, row, "atom_id_2")) else { return false };
        let other = if first == name { second } else if second == name { first } else { continue };
        if !text(&bonds, row, "value_order").is_some_and(|value| value.eq_ignore_ascii_case("sing"))
            || !text(&bonds, row, "pdbx_metal_coordination_flag").is_some_and(|value| value.eq_ignore_ascii_case("Y")) {
            return false;
        }
        neighbors.insert(other);
    }
    neighbors.into_iter().collect::<Vec<_>>() == candidate.metals
}

/// Read hydride charge and both metal bonds from the selected dictionary source.
pub(crate) fn dictionary(library: &MonLib, comp: &Comp) -> Result<Arc<Definitions>, String> {
    let key = (library.root.clone(), comp.source.clone(), comp.id.clone(), comp.applied_mods.clone(), comp.from_ccd);
    let cache = DEFINITIONS.get_or_init(Mutex::default);
    if let Some(value) = cache.lock().unwrap_or_else(PoisonError::into_inner).get(&key) {
        return Ok(Arc::clone(value));
    }
    let mut result = candidates(comp);
    if !result.is_empty() {
        let bytes = std::fs::read(&comp.source).map_err(|error| format!("Cannot read hydride dictionary {}: {error}", comp.id))?;
        let document = cif::parse(&bytes).map_err(|error| format!("Cannot parse hydride dictionary {}: {error}", comp.id))?;
        result.retain(|name, definition| source_definition(comp, name, definition, &document));
    }
    let result = Arc::new(result);
    cache.lock().unwrap_or_else(PoisonError::into_inner).insert(key, Arc::clone(&result));
    Ok(result)
}

fn site(source: &Source, group: usize, rows: &[usize], name: &str, definition: &Definition) -> Result<Site, String> {
    let mut metal_rows = [usize::MAX; 2];
    for (index, name) in definition.metals.iter().enumerate() {
        let found: Vec<_> = rows.iter().copied().filter(|&row| source.atoms[row].name == *name).collect();
        let [row] = found[..] else { return Err(format!("Bridge metal {name} is absent or ambiguous")) };
        let atom = &source.atoms[row];
        if atom.protein || !atom.element.eq_ignore_ascii_case(&definition.metal_elements[index])
            || !residue::finite(atom.position) { return Err(format!("Bridge metal {name} has invalid type or coordinates")); }
        metal_rows[index] = row;
    }
    let observed: Vec<_> = rows.iter().copied().filter(|&row| observed_name(name, &source.atoms[row]).is_some()).collect();
    if observed.len() > 1 { return Err(format!("Bridge hydrogen {name} is ambiguous")); }
    if observed.iter().any(|&row| {
        let atom = &source.atoms[row];
        atom.protein || observed_name(name, atom) != Some(true) || !residue::finite(atom.position)
    }) {
        return Err(format!("Bridge hydrogen {name} has invalid type or coordinates"));
    }
    if let Some(&row) = observed.first() {
        for &(_, parent) in source.hydrogen_bonds.iter().filter(|&&(child, _)| child == row) {
            if parent.is_none_or(|parent| !metal_rows.contains(&parent)) {
                return Err(format!("Bridge hydrogen {name} has an absent or conflicting declared metal"));
            }
        }
    }
    Ok(Site { group, hydrogen_name: name.into(), hydrogen_row: observed.first().copied(),
        metal_rows, formal_charge: definition.formal_charge, kind: definition.kind })
}

/// Resolve bridge sites even when the heavy source omits their H coordinate.
pub(super) fn for_source(library: &MonLib, source: &Source) -> Result<Vec<Site>, String> {
    let mut result = Vec::new();
    let mut selections = super::parents::Selections::new();
    for (group, rows) in source.groups.iter().enumerate() {
        let Some(&first) = rows.first() else { continue };
        if source.atoms[first].protein { continue; }
        let names = rows.iter().map(|&row| &source.atoms[row]).filter(|atom| !hydrogen(&atom.element))
            .map(|atom| &atom.name);
        let Some(definitions) = selections.get(library, source.atoms[first].component.as_str(), names,
            |comp| dictionary(library, &comp))? else { continue };
        let mut names: Vec<_> = definitions.keys().collect();
        names.sort_unstable();
        for name in names { result.push(site(source, group, rows, name, &definitions[name])?); }
    }
    Ok(result)
}

#[cfg(test)]
#[path = "bridges_tests.rs"]
mod tests;
