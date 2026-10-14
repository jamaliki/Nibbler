//! Dense observations under a schema supplied by the caller.

pub(crate) mod afdb;
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
