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

use crate::cif::CifDocument;
pub(crate) use record::Tokens;
pub(crate) use prepared::PreparedSource;
pub(crate) use schema::Schema;
use std::collections::HashMap;
#[cfg(feature = "reduce3")]
use std::collections::BTreeSet;
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
            (!residue::hydrogen(&atom.element)).then(|| atom.name.clone())
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
                    first.component.clone(),
                    rows.iter()
                        .filter_map(|&row| {
                            let atom = &source.atoms[row];
                            (!residue::hydrogen(&atom.element)).then(|| atom.name.clone())
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
    Ok(RetainedParents::read(document, expected_keys, model, chem_data)?.required)
}

#[cfg(feature = "reduce3")]
pub(crate) struct RetainedParents {
    source: source::Source,
    pub required: Vec<u32>,
}

#[cfg(feature = "reduce3")]
impl RetainedParents {
    pub fn read(document: &CifDocument, expected_keys: &str, model: i64,
        chem_data: Option<&Path>) -> Result<Self, String> {
        let expected: Vec<[String; 5]> =
            serde_json::from_str(expected_keys).map_err(|error| error.to_string())?;
        let wanted: BTreeSet<_> = expected.iter().cloned().collect();
        if wanted.len() != expected.len() {
            return Err("Retained hydrogen-parent keys contain duplicates".into());
        }
        let source = source::read(document, model)?;
        let required = required_sources(&source, &wanted, chem_data)?;
        Ok(Self { source, required })
    }

    pub fn ligand_sources(&self) -> Vec<u32> {
        let required = self.required.iter().copied().collect::<BTreeSet<_>>();
        self.source.atoms.iter().filter(|atom|
            !atom.protein && !residue::hydrogen(&atom.element) && required.contains(&atom.source_row))
            .map(|atom| atom.source_row).collect()
    }

    pub fn bridge_sources(&self, chem_data: Option<&Path>) -> Result<Vec<(String, [u32; 2])>, String> {
        if !bridges::possible_source(&self.source) { return Ok(Vec::new()); }
        let library = crate::reduce::load_monomer_library(chem_data).map_err(|error| error.to_string())?;
        Ok(bridges::for_source(&library, &self.source)?.into_iter().map(|site| {
            (site.hydrogen_name, site.metal_rows.map(|row| self.source.atoms[row].source_row))
        }).collect())
    }
}

#[cfg(feature = "reduce3")]
fn required_sources(source: &source::Source, wanted: &BTreeSet<[String; 5]>,
    chem_data: Option<&Path>) -> Result<Vec<u32>, String> {
    let mut found = HashMap::new();
    let mut selected = Vec::new();
    for atom in &source.atoms {
        if residue::hydrogen(&atom.element) || !residue::finite(atom.position) {
            continue;
        }
        let key = record::token_key(atom, !atom.protein);
        if !wanted.contains(&key) {
            continue;
        }
        if let Some(group) = found.insert(key.clone(), atom.group)
            && group != atom.group
        {
            return Err(format!(
                "Retained hydrogen-parent key matches several residues: {key:?}"
            ));
        }
        selected.push(atom.source_row);
    }
    if bridges::possible_source(source)
        && wanted.iter().any(|key| key[4] != "protein" && !found.contains_key(key)) {
        required_bridge_sources(source, wanted, chem_data, &mut found, &mut selected)?;
    }
    if found.len() != wanted.len() {
        let missing: Vec<_> = wanted
            .iter()
            .filter(|key| !found.contains_key(*key))
            .take(3)
            .collect();
        return Err(format!(
            "Retained hydrogen-parent keys are absent from the selected model: {missing:?}"
        ));
    }
    selected.sort_unstable();
    selected.dedup();
    Ok(selected)
}

#[cfg(feature = "reduce3")]
fn required_bridge_sources(source: &source::Source, wanted: &BTreeSet<[String; 5]>,
    chem_data: Option<&Path>, found: &mut HashMap<[String; 5], usize>,
    selected: &mut Vec<u32>) -> Result<(), String> {
    let library = crate::reduce::load_monomer_library(chem_data).map_err(|error| error.to_string())?;
    for site in bridges::for_source(&library, source)? {
        let atom = &source.atoms[site.metal_rows[0]];
        let mut key = record::token_key(atom, true);
        key[4] = site.hydrogen_name;
        if wanted.contains(&key) {
            if found.insert(key.clone(), site.group).is_some() {
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
