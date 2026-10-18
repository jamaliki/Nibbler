//! Resolve hydrogen dictionaries for each retained ligand residue.

use super::source::Source;
use super::{ResidueParents, residue};
use crate::hash::FastMap;
use crate::pdbx::ShortText;
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

/// Heavy-atom name lists seen for one component, each with its dictionary.
type Choices<'a, T> = Vec<(Vec<&'a ShortText>, Option<T>)>;

/// Hydrogen dictionaries chosen for one source, by component and heavy-atom names.
///
/// `hydrogen_dictionary` depends only on these inputs. RNA chains repeat a few
/// components with a few name lists, so most groups reuse an earlier choice. Each
/// component keeps its name lists, which compare as inline text without hashing.
pub(super) struct Selections<'a, T> {
    found: FastMap<&'a str, Choices<'a, T>>,
    names: Vec<&'a ShortText>,
}

impl<'a, T: Clone> Selections<'a, T> {
    pub fn new() -> Self {
        Self { found: FastMap::default(), names: Vec::new() }
    }

    /// The dictionary for `component` and `names`, passed through `resolve` once per distinct input.
    pub fn get(&mut self, library: &MonLib, component: &'a str, names: impl Iterator<Item = &'a ShortText>,
        resolve: impl FnOnce(Arc<Comp>) -> Result<T, String>) -> Result<Option<T>, String> {
        self.names.clear();
        self.names.extend(names);
        let known = self.found.entry(component).or_default();
        // Field compares inline; a slice compare calls a loop that is not inlined.
        let names = &self.names;
        if let Some((_, found)) = known.iter().find(|(known, _)| known.len() == names.len()
            && known.iter().zip(names).all(|(left, right)| **left == **right)) {
            return Ok(found.clone());
        }
        let owned: Vec<String> = self.names.iter().map(|name| name.as_str().to_owned()).collect();
        let found = reduce3::hplace::hydrogen_dictionary(library, component, &owned, false).map(resolve).transpose()?;
        known.push((self.names.clone(), found.clone()));
        Ok(found)
    }
}

/// The hydrogen dictionary of a nonprotein residue, chosen from its heavy-atom names.
fn selected<'a, T: Clone>(selections: &mut Selections<'a, T>, library: &MonLib, source: &'a Source, rows: &[usize],
    resolve: impl FnOnce(Arc<Comp>) -> Result<T, String>) -> Result<Option<T>, String> {
    let Some(&first) = rows.first() else { return Ok(None) };
    let first = &source.atoms[first];
    if first.protein {
        return Ok(None);
    }
    let names = rows.iter().map(|&row| &source.atoms[row]).filter(|atom| !residue::hydrogen(&atom.element))
        .map(|atom| &atom.name);
    selections.get(library, first.component.as_str(), names, resolve)
}

pub(super) fn for_source(library: &MonLib, source: &Source) -> Result<ResidueParents, String> {
    let mut result = HashMap::new();
    let mut selections = Selections::new();
    for (group, rows) in source.groups.iter().enumerate() {
        // `selected` returns nothing for a protein group, so skip its atoms before reading them.
        if rows.first().is_some_and(|&first| !source.atoms[first].protein)
            && rows.iter().any(|&row| residue::hydrogen(&source.atoms[row].element))
            && let Some(parents) = selected(&mut selections, library, source, rows,
                |comp| dictionary(library, &comp))? {
            result.insert(group, parents);
        }
    }
    Ok(result)
}

pub(super) fn missing(library: &MonLib, source: &Source) -> Vec<String> {
    let mut selections = Selections::new();
    let mut result = BTreeSet::new();
    for rows in &source.groups {
        let Some(&first) = rows.first() else { continue };
        let first = &source.atoms[first];
        if !first.protein && matches!(selected(&mut selections, library, source, rows, |_| Ok(())), Ok(None)) {
            result.insert(first.component.as_str());
        }
    }
    result.into_iter().map(str::to_owned).collect()
}
