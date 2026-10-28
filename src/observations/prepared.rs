//! One immutable source snapshot for component queries and target packing.

use super::{HydrogenParents, LigandFeatures, PackOptions, Schema, TokenTable, Tokens, source};
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
        Ok(Self {
            source: source::read(document, model)?,
            model,
        })
    }

    pub(crate) fn components(&self) -> Vec<String> {
        // Atoms of one residue sit together, so skip repeats before the set lookup.
        let mut names = BTreeSet::new();
        let mut previous = None;
        for atom in self.source.atoms.iter().filter(|atom| !atom.protein) {
            if previous != Some(&atom.component) {
                names.insert(atom.component.as_str());
                previous = Some(&atom.component);
            }
        }
        names.into_iter().map(str::to_owned).collect()
    }

    #[cfg(feature = "reduce3")]
    pub(crate) fn missing_hydrogen_dictionaries(
        &self,
        chem_data: &Path,
    ) -> Result<Vec<String>, String> {
        let library = crate::reduce::load_monomer_library(Some(chem_data))
            .map_err(|error| error.to_string())?;
        Ok(super::parents::missing(&library, &self.source))
    }

    #[allow(clippy::too_many_arguments)]
    pub(crate) fn pack(
        &self,
        schema: &Schema,
        features: &LigandFeatures,
        hydrogen_parents: &HydrogenParents,
        table: &TokenTable,
        excluded: &[String],
        options: &PackOptions<'_>,
    ) -> Result<Tokens, String> {
        super::pack_source(
            &self.source,
            schema,
            features,
            hydrogen_parents,
            table,
            excluded,
            options,
        )
    }
}
