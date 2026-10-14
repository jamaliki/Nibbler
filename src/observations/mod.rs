//! Dense observations under a schema supplied by the caller.

pub(crate) mod afdb;
pub(crate) mod afdb_eligibility;
pub(crate) mod crops;
mod record;
mod residue;
mod schema;
mod source;
mod tokens;
mod topology;

use crate::cif::CifDocument;
pub(crate) use record::Tokens;
pub(crate) use schema::Schema;
use std::collections::{BTreeSet, HashMap};

pub(crate) type LigandFeatures = HashMap<String, HashMap<String, [i16; 4]>>;
pub(crate) type HydrogenParents = HashMap<String, HashMap<String, String>>;

pub(crate) struct PackOptions<'a> {
    pub model: i64,
    pub complete: bool,
    pub topology_policy: &'a str,
    pub protein_only: bool,
}

pub(crate) fn components(document: &CifDocument, model: i64) -> Result<Vec<String>, String> {
    let source = source::read(document, model)?;
    Ok(source
        .atoms
        .iter()
        .filter(|atom| !atom.protein)
        .map(|atom| atom.component.clone())
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect())
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
) -> Result<Vec<u32>, String> {
    let expected: Vec<[String; 5]> =
        serde_json::from_str(expected_keys).map_err(|error| error.to_string())?;
    let wanted: BTreeSet<_> = expected.iter().cloned().collect();
    if wanted.len() != expected.len() {
        return Err("Retained hydrogen-parent keys contain duplicates".into());
    }
    let source = source::read(document, model)?;
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
    Ok(selected)
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
    tokens::pack(
        &source,
        schema,
        features,
        hydrogen_parents,
        elements,
        excluded,
        options,
    )
}
