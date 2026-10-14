//! Resolve hydrogen dictionaries for each retained ligand residue.

use super::source::Source;
use super::{ResidueParents, residue};
use reduce3::monlib::{Comp, MonLib};
use std::collections::{BTreeSet, HashMap};
use std::path::PathBuf;
use std::sync::{Arc, Mutex, OnceLock, PoisonError};

pub(crate) type Parents = HashMap<String, String>;
type CacheKey = (PathBuf, PathBuf, String, Vec<String>, bool);
static PARENTS: OnceLock<Mutex<HashMap<CacheKey, Arc<Parents>>>> = OnceLock::new();

fn hydrogen(element: &str) -> bool {
    ["H", "D", "T"].iter().any(|value| element.eq_ignore_ascii_case(value))
}

pub(crate) fn dictionary(library: &MonLib, dictionary: &Comp) -> Result<Arc<Parents>, String> {
    let key = (library.root.clone(), dictionary.source.clone(), dictionary.id.clone(),
        dictionary.applied_mods.clone(), dictionary.from_ccd);
    let cache = PARENTS.get_or_init(Mutex::default);
    if let Some(parents) = cache.lock().unwrap_or_else(PoisonError::into_inner).get(&key) {
        return Ok(Arc::clone(parents));
    }
    let mut parents = HashMap::new();
    let bridges = super::bridges::dictionary(library, dictionary)?;
    for bond in &dictionary.bonds {
        let (Some(first), Some(second)) = (dictionary.atom(&bond.a1), dictionary.atom(&bond.a2)) else {
            return Err(format!("Reduce3 dictionary {} has a bond with an absent atom", dictionary.id));
        };
        let pair = match (hydrogen(&first.type_symbol), hydrogen(&second.type_symbol)) {
            (true, false) => Some((&bond.a1, &bond.a2)),
            (false, true) => Some((&bond.a2, &bond.a1)),
            _ => None,
        };
        if let Some((child, parent)) = pair
            && !bridges.contains_key(child)
            && let Some(previous) = parents.insert(child.clone(), parent.clone())
            && previous != *parent {
            return Err(format!("Reduce3 dictionary {} hydrogen {child} has conflicting parents", dictionary.id));
        }
    }
    parents.extend(reduce3::hplace::terminal::hydrogen_parents(dictionary));
    let parents = Arc::new(parents);
    cache.lock().unwrap_or_else(PoisonError::into_inner).insert(key, Arc::clone(&parents));
    Ok(parents)
}

fn selected(library: &MonLib, source: &Source, rows: &[usize]) -> Option<Arc<Comp>> {
    let first = &source.atoms[*rows.first()?];
    if first.protein {
        return None;
    }
    let names = rows.iter().filter_map(|&row| {
        let atom = &source.atoms[row];
        (!residue::hydrogen(&atom.element)).then(|| atom.name.clone())
    }).collect::<Vec<_>>();
    reduce3::hplace::hydrogen_dictionary(library, &first.component, &names, false)
}

pub(super) fn for_source(library: &MonLib, source: &Source) -> Result<ResidueParents, String> {
    let mut result = HashMap::new();
    for (group, rows) in source.groups.iter().enumerate() {
        if rows.iter().any(|&row| residue::hydrogen(&source.atoms[row].element))
            && let Some(comp) = selected(library, source, rows) {
            result.insert(group, dictionary(library, &comp)?);
        }
    }
    Ok(result)
}

pub(super) fn missing(library: &MonLib, source: &Source) -> Vec<String> {
    source.groups.iter().filter_map(|rows| {
        let first = &source.atoms[*rows.first()?];
        (!first.protein && selected(library, source, rows).is_none()).then(|| first.component.clone())
    }).collect::<BTreeSet<_>>().into_iter().collect()
}
