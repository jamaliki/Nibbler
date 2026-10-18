//! Dense observations under a schema supplied by the caller.

pub(crate) mod afdb;
pub(crate) mod afdb_eligibility;
pub(crate) mod crops;
#[cfg(feature = "reduce3")]
pub(crate) mod bridges;
#[cfg(feature = "reduce3")]
pub(crate) mod parents;
mod record;
mod prepared;
mod residue;
mod schema;
pub(crate) mod source;
mod tokens;
mod topology;
mod continuity;
pub(crate) mod uncertainty;
pub(crate) mod product;
#[cfg(feature = "reduce3")]
pub(crate) mod molecular;

use crate::cif::CifDocument;
pub(crate) use record::Tokens;
pub(crate) use prepared::PreparedSource;
pub(crate) use schema::Schema;
use std::collections::HashMap;
#[cfg(feature = "reduce3")]
use crate::hash::{FastMap, FastSet};
#[cfg(feature = "reduce3")]
use std::borrow::Cow;
#[cfg(feature = "reduce3")]
use std::collections::BTreeSet;
#[cfg(feature = "reduce3")]
use std::fmt::Write;
use std::path::Path;
use std::sync::Arc;

pub(crate) type LigandFeatures = HashMap<String, HashMap<String, [i16; 4]>>;
pub(crate) type HydrogenParents = HashMap<String, HashMap<String, String>>;
pub(crate) type ResidueParents = HashMap<usize, Arc<HashMap<String, String>>>;

pub(crate) struct PackOptions<'a> {
    pub model: i64,
    pub complete: bool,
    pub topology_policy: &'a str,
    pub protein_only: bool,
    pub chem_data: Option<&'a Path>,
}

#[cfg(feature = "reduce3")]
pub(crate) fn missing_hydrogen_dictionaries(document: &CifDocument, model: i64, chem_data: &Path) -> Result<Vec<String>, String> {
    PreparedSource::read(document, model)?.missing_hydrogen_dictionaries(chem_data)
}

pub(crate) fn components(document: &CifDocument, model: i64) -> Result<Vec<String>, String> {
    Ok(PreparedSource::read(document, model)?.components())
}

#[cfg(feature = "reduce3")]
pub(crate) fn chemical_sources(document: &CifDocument, model: i64, chem_data: &Path)
    -> Result<Vec<std::path::PathBuf>, String> {
    let source = source::read(document, model)?;
    let library = crate::reduce::load_monomer_library(Some(chem_data)).map_err(|error| error.to_string())?;
    let mut files = BTreeSet::new();
    for name in ["geostd/list/mon_lib_list.cif", "geostd/list/geostd_list.cif",
        "mon_lib/list/mon_lib_list.cif", "geostd/ener_lib.cif"] {
        let file = library.root.join(name);
        if file.is_file() { files.insert(file); }
    }
    for rows in &source.groups {
        let Some(&row) = rows.first() else { continue; };
        let name = &source.atoms[row].component;
        let names = rows.iter().filter_map(|&row| {
            let atom = &source.atoms[row];
            (!residue::hydrogen(&atom.element)).then(|| atom.name.as_str().to_owned())
        }).collect::<Vec<_>>();
        if let Some(dictionary) = reduce3::hplace::hydrogen_dictionary(&library, name, &names, false) {
            files.insert(dictionary.source.clone());
        }
        if name.bytes().all(|value| value.is_ascii_alphanumeric()) && !name.is_empty() {
            let file = library.root.join("chemical_components").join(name[..1].to_ascii_lowercase())
                .join(format!("data_{name}.cif"));
            if file.is_file() { files.insert(file); }
        }
    }
    Ok(files.into_iter().collect())
}

#[cfg(feature = "reduce3")]
pub(crate) fn ligand_groups(
    document: &CifDocument,
    model: i64,
) -> Result<Vec<(String, Vec<String>)>, String> {
    let source = source::read(document, model)?;
    Ok(source
        .groups
        .iter()
        .filter_map(|rows| {
            let first = &source.atoms[*rows.first()?];
            (!first.protein).then(|| {
                (
                    first.component.as_str().to_owned(),
                    rows.iter()
                        .filter_map(|&row| {
                            let atom = &source.atoms[row];
                            (!residue::hydrogen(&atom.element)).then(|| atom.name.as_str().to_owned())
                        })
                        .collect(),
                )
            })
        })
        .collect())
}

#[cfg(feature = "reduce3")]
pub(crate) fn required_parent_sources(
    document: &CifDocument,
    expected_keys: &str,
    model: i64,
    chem_data: Option<&Path>,
) -> Result<Vec<u32>, String> {
    let keys = ParentKeys::parse(expected_keys, &[])?;
    let source = source::read(document, model)?;
    Ok(RetainedParents::new(&source, &keys, chem_data)?.required)
}

/// Retained hydrogen-parent keys, checked to be distinct.
#[cfg(feature = "reduce3")]
pub(crate) struct ParentKeys<'a>(Vec<[Cow<'a, str>; 5]>);

/// One key field, borrowed from the JSON text unless it holds an escape.
#[cfg(feature = "reduce3")]
#[derive(serde::Deserialize)]
struct KeyField<'a>(#[serde(borrow)] Cow<'a, str>);

#[cfg(feature = "reduce3")]
impl<'a> ParentKeys<'a> {
    /// Parse the key list and drop keys whose component (the fourth field) is in `ignored`.
    pub fn parse(text: &'a str, ignored: &[String]) -> Result<Self, String> {
        let keys = match compact_keys(text) {
            Some(keys) => keys,
            None => {
                let keys: Vec<[KeyField<'a>; 5]> = serde_json::from_str(text).map_err(|error| error.to_string())?;
                keys.into_iter().map(|key| key.map(|field| field.0)).collect()
            }
        };
        let keys: Vec<_> = keys.into_iter().filter(|key| !ignored.iter().any(|name| *name == key[3])).collect();
        let mut seen = FastSet::with_capacity_and_hasher(keys.len(), Default::default());
        if !keys.iter().all(|key| seen.insert(key)) {
            return Err("Retained hydrogen-parent keys contain duplicates".into());
        }
        Ok(Self(keys))
    }
}

/// Split key JSON in the compact form the index writes (`[["A","1","","ALA","protein"],...]`, with no
/// spaces, escapes or control characters) without serde. Any other text returns None, and serde then
/// parses or rejects it.
#[cfg(feature = "reduce3")]
fn compact_keys(text: &str) -> Option<Vec<[Cow<'_, str>; 5]>> {
    let mut rest = text.strip_prefix('[')?.strip_suffix(']')?;
    let mut keys = Vec::new();
    while !rest.is_empty() {
        if !keys.is_empty() {
            rest = rest.strip_prefix(',')?;
        }
        rest = rest.strip_prefix('[')?;
        let mut key: [Cow<'_, str>; 5] = Default::default();
        for (position, field) in key.iter_mut().enumerate() {
            if position > 0 {
                rest = rest.strip_prefix(',')?;
            }
            rest = rest.strip_prefix('"')?;
            let (value, tail) = rest.split_once('"')?;
            if value.bytes().any(|byte| byte == b'\\' || byte < 0x20) {
                return None;
            }
            *field = Cow::Borrowed(value);
            rest = tail;
        }
        rest = rest.strip_prefix(']')?;
        keys.push(key);
    }
    Some(keys)
}

/// The source rows of retained hydrogen parents, found in an already read source.
#[cfg(feature = "reduce3")]
pub(crate) struct RetainedParents<'a> {
    source: &'a source::Source,
    pub required: Vec<u32>,
}

#[cfg(feature = "reduce3")]
impl<'a> RetainedParents<'a> {
    pub fn new(source: &'a source::Source, keys: &ParentKeys<'_>, chem_data: Option<&Path>) -> Result<Self, String> {
        let required = required_sources(source, &keys.0, chem_data)?;
        Ok(Self { source, required })
    }

    pub fn ligand_sources(&self) -> Vec<u32> {
        // `required` is sorted and free of repeats.
        self.source.atoms.iter().filter(|atom| !atom.protein && !residue::hydrogen(&atom.element)
            && self.required.binary_search(&atom.source_row).is_ok())
            .map(|atom| atom.source_row).collect()
    }

    pub fn bridge_sources(&self, chem_data: Option<&Path>) -> Result<Vec<(String, [u32; 2])>, String> {
        if !bridges::possible_source(self.source) { return Ok(Vec::new()); }
        let library = crate::reduce::load_monomer_library(chem_data).map_err(|error| error.to_string())?;
        Ok(bridges::for_source(&library, self.source)?.into_iter().map(|site| {
            (site.hydrogen_name, site.metal_rows.map(|row| self.source.atoms[row].source_row))
        }).collect())
    }
}

/// Map each wanted key to its position in the parsed key list.
#[cfg(feature = "reduce3")]
fn key_positions<'k>(expected: &'k [[Cow<'_, str>; 5]]) -> FastMap<[&'k str; 5], usize> {
    let mut positions = FastMap::with_capacity_and_hasher(expected.len(), Default::default());
    for (index, key) in expected.iter().enumerate() {
        positions.insert(key.each_ref().map(|field| field.as_ref()), index);
    }
    positions
}

#[cfg(feature = "reduce3")]
fn required_sources(source: &source::Source, expected: &[[Cow<'_, str>; 5]],
    chem_data: Option<&Path>) -> Result<Vec<u32>, String> {
    let wanted = key_positions(expected);
    // The residue group each wanted key matched, by key position.
    let mut found: Vec<Option<usize>> = vec![None; expected.len()];
    let mut selected = Vec::new();
    let mut number = String::new();
    // The atoms of a protein residue share one key, so look it up once per run of such atoms.
    let mut last_protein: Option<(&source::Atom, Option<usize>)> = None;
    for atom in &source.atoms {
        if residue::hydrogen(&atom.element) || !residue::finite(atom.position) {
            continue;
        }
        let ligand = !atom.protein;
        let looked_up = match last_protein {
            Some((previous, index)) if !ligand && previous.sequence == atom.sequence
                && previous.chain.as_str() == atom.chain.as_str() && previous.insertion.as_str() == atom.insertion.as_str()
                && previous.component.as_str() == atom.component.as_str() => index,
            _ => {
                let sequence = if ligand {
                    atom.auth_sequence.as_str()
                } else {
                    number.clear();
                    // Writing to a String cannot fail.
                    let _ = write!(number, "{}", atom.sequence);
                    number.as_str()
                };
                let key = [atom.chain.as_str(), sequence, atom.insertion.as_str(), atom.component.as_str(),
                    if ligand { atom.name.as_str() } else { "protein" }];
                let index = wanted.get(&key).copied();
                if !ligand {
                    last_protein = Some((atom, index));
                }
                index
            }
        };
        let Some(index) = looked_up else {
            continue;
        };
        if let Some(group) = found[index].replace(atom.group)
            && group != atom.group
        {
            return Err(format!(
                "Retained hydrogen-parent key matches several residues: {:?}", expected[index]
            ));
        }
        selected.push(atom.source_row);
    }
    if bridges::possible_source(source)
        && expected.iter().zip(&found).any(|(key, group)| key[4] != "protein" && group.is_none()) {
        required_bridge_sources(source, &wanted, chem_data, &mut found, &mut selected)?;
    }
    if found.iter().any(Option::is_none) {
        let mut missing: Vec<_> = expected.iter().zip(&found)
            .filter(|(_, group)| group.is_none()).map(|(key, _)| key).collect();
        missing.sort_unstable();
        missing.truncate(3);
        return Err(format!(
            "Retained hydrogen-parent keys are absent from the selected model: {missing:?}"
        ));
    }
    selected.sort_unstable();
    selected.dedup();
    Ok(selected)
}

#[cfg(feature = "reduce3")]
fn required_bridge_sources(source: &source::Source, wanted: &FastMap<[&str; 5], usize>,
    chem_data: Option<&Path>, found: &mut [Option<usize>], selected: &mut Vec<u32>) -> Result<(), String> {
    let library = crate::reduce::load_monomer_library(chem_data).map_err(|error| error.to_string())?;
    for site in bridges::for_source(&library, source)? {
        let atom = &source.atoms[site.metal_rows[0]];
        let mut key = record::token_key(atom, true);
        key[4] = site.hydrogen_name;
        if let Some(&index) = wanted.get(&key.each_ref().map(String::as_str)) {
            if found[index].replace(site.group).is_some() {
                return Err(format!("Retained bridge key matches several residues: {key:?}"));
            }
            selected.extend(site.metal_rows.map(|row| source.atoms[row].source_row));
        }
    }
    Ok(())
}

pub(crate) fn pack(
    document: &CifDocument,
    schema: &Schema,
    features: &LigandFeatures,
    hydrogen_parents: &HydrogenParents,
    elements: &[String],
    excluded: &[String],
    options: &PackOptions<'_>,
) -> Result<Tokens, String> {
    let source = source::read(document, options.model)?;
    pack_source(&source, schema, features, hydrogen_parents, elements, excluded, options)
}

fn pack_source(
    source: &source::Source, schema: &Schema, features: &LigandFeatures,
    hydrogen_parents: &HydrogenParents, elements: &[String], excluded: &[String],
    options: &PackOptions<'_>,
) -> Result<Tokens, String> {
    #[cfg(feature = "reduce3")]
    let mut bridge_sites = Vec::new();
    #[cfg(feature = "reduce3")]
    let residue_parents = if options.chem_data.is_some()
        || (!options.protein_only && bridges::possible_source(source)) {
            let library = crate::reduce::load_monomer_library(options.chem_data).map_err(|error| error.to_string())?;
            if !options.protein_only {
                bridge_sites = bridges::for_source(&library, source)?;
                if !bridge_sites.is_empty() && !elements.iter().any(|element| element == "H") {
                    return Err("Metal hydride packing requires a ligand contract with H token 65".into());
                }
            }
            parents::for_source(&library, source)?
    } else {
        HashMap::new()
    };
    #[cfg(not(feature = "reduce3"))]
    let residue_parents = {
        if options.chem_data.is_some() {
            return Err("Hydrogen dictionary packing requires the Reduce3 feature".into());
        }
        HashMap::new()
    };
    tokens::pack(
        source,
        schema,
        features,
        hydrogen_parents,
        &residue_parents,
        #[cfg(feature = "reduce3")]
        &bridge_sites,
        elements,
        excluded,
        options,
    )
}
