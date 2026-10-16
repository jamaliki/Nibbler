//! One immutable source snapshot for component queries and target packing.

use super::{HydrogenParents, LigandFeatures, PackOptions, Schema, Tokens, source};
use crate::cif::CifDocument;
use std::collections::BTreeSet;
#[cfg(feature = "reduce3")]
use std::path::Path;

pub(crate) struct PreparedSource {
    pub(crate) source: source::Source,
    pub model: i64,
}

impl PreparedSource {
    pub(crate) fn read(document: &CifDocument, model: i64) -> Result<Self, String> {
        Ok(Self { source: source::read(document, model)?, model })
    }

    pub(crate) fn components(&self) -> Vec<String> {
        self.source.atoms.iter().filter(|atom| !atom.protein)
            .map(|atom| atom.component.clone()).collect::<BTreeSet<_>>().into_iter().collect()
    }

    #[cfg(feature = "reduce3")]
    pub(crate) fn missing_hydrogen_dictionaries(&self, chem_data: &Path) -> Result<Vec<String>, String> {
        let library = crate::reduce::load_monomer_library(Some(chem_data)).map_err(|error| error.to_string())?;
        Ok(super::parents::missing(&library, &self.source))
    }

    #[allow(clippy::too_many_arguments)]
    pub(crate) fn pack(&self, schema: &Schema, features: &LigandFeatures,
        hydrogen_parents: &HydrogenParents, elements: &[String], excluded: &[String],
        options: &PackOptions<'_>) -> Result<Tokens, String> {
        super::pack_source(&self.source, schema, features, hydrogen_parents, elements, excluded, options)
    }
}
