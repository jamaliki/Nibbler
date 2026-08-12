use std::error::Error;
use std::fmt::{self, Display, Formatter};
use std::str::FromStr;
use std::sync::OnceLock;

use regex::{Regex, RegexBuilder};

use super::dictionary::Dictionary;
use super::schema_artifact::decode_dictionary;

const PDBX_VERSION: &str = "5.416";
const PDBX_SHA256: &str = "78c038f20f69601b86f128c0d543576f4284da67c36284c344669da0f4f5e27a";
const MODELCIF_VERSION: &str = "1.4.9";
const MODELCIF_SHA256: &str = "aeb205bbf23d459089c694eec004978c78497e5f5d3dd5642120dd0d62e9b844";

static PDBX: OnceLock<Result<LoadedSchema, SchemaError>> = OnceLock::new();
static MODELCIF: OnceLock<Result<LoadedSchema, SchemaError>> = OnceLock::new();

/// One built-in, lock-pinned DDL2 schema.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum SchemaName {
    /// PDBx/mmCIF 5.416.
    Pdbx,
    /// ModelCIF 1.4.9, including its pinned PDBx base.
    ModelCif,
}

impl SchemaName {
    /// Return the stable Python/API selector.
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Pdbx => "pdbx",
            Self::ModelCif => "modelcif",
        }
    }
}

impl FromStr for SchemaName {
    type Err = SchemaError;

    fn from_str(value: &str) -> Result<Self, Self::Err> {
        match value.to_ascii_lowercase().as_str() {
            "pdbx" | "mmcif" => Ok(Self::Pdbx),
            "modelcif" => Ok(Self::ModelCif),
            _ => Err(SchemaError::new(
                "CIF_SCHEMA_UNKNOWN",
                format!("unknown schema {value:?}; expected 'pdbx' or 'modelcif'"),
            )),
        }
    }
}

/// One decoded schema and its precompiled DDL2 type expressions.
pub(crate) struct LoadedSchema {
    dictionary: Dictionary,
    patterns: std::collections::BTreeMap<String, Regex>,
}

impl LoadedSchema {
    /// Return the immutable compiled dictionary.
    #[must_use]
    pub(crate) const fn dictionary(&self) -> &Dictionary {
        &self.dictionary
    }

    pub(crate) fn pattern(&self, type_code: &str) -> Option<&Regex> {
        self.patterns.get(&type_code.to_ascii_lowercase())
    }
}

/// A schema selector, artifact, lock, or regular-expression failure.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct SchemaError {
    code: &'static str,
    message: String,
}

impl SchemaError {
    pub(crate) fn new(code: &'static str, message: impl Into<String>) -> Self {
        Self {
            code,
            message: message.into(),
        }
    }

    /// Return the stable machine-readable error code.
    #[must_use]
    pub const fn code(&self) -> &'static str {
        self.code
    }

    /// Return the human-readable failure.
    #[must_use]
    pub fn message(&self) -> &str {
        &self.message
    }
}

impl Display for SchemaError {
    fn fmt(&self, formatter: &mut Formatter<'_>) -> fmt::Result {
        write!(formatter, "{}: {}", self.code, self.message)
    }
}

impl Error for SchemaError {}

/// Return one lazily decoded, lock-verified built-in schema.
///
/// # Errors
///
/// Returns [`SchemaError`] if the embedded artifact is corrupt, disagrees with the
/// source lock, or contains a DDL2 expression unsupported by the regex engine.
pub(crate) fn loaded_schema(name: SchemaName) -> Result<&'static LoadedSchema, SchemaError> {
    let result = match name {
        SchemaName::Pdbx => PDBX.get_or_init(|| {
            load(
                include_bytes!("../../schemas/compiled/pdbx.nbs"),
                SchemaName::Pdbx,
                PDBX_VERSION,
                PDBX_SHA256,
            )
        }),
        SchemaName::ModelCif => MODELCIF.get_or_init(|| {
            load(
                include_bytes!("../../schemas/compiled/modelcif.nbs"),
                SchemaName::ModelCif,
                MODELCIF_VERSION,
                MODELCIF_SHA256,
            )
        }),
    };
    result.as_ref().map_err(Clone::clone)
}

fn load(
    artifact: &[u8],
    name: SchemaName,
    version: &str,
    sha256: &str,
) -> Result<LoadedSchema, SchemaError> {
    let dictionary = decode_dictionary(artifact)
        .map_err(|error| SchemaError::new("CIF_SCHEMA_ARTIFACT", error.to_string()))?;
    let metadata = dictionary.metadata();
    if metadata.schema_name() != name.as_str()
        || metadata.version() != version
        || metadata.source_sha256() != sha256
    {
        return Err(SchemaError::new(
            "CIF_SCHEMA_LOCK_MISMATCH",
            format!(
                "embedded {} artifact does not match its source lock",
                name.as_str()
            ),
        ));
    }

    let mut patterns = std::collections::BTreeMap::new();
    for item_type in dictionary.types.values() {
        let pattern = normalize_ddl2_pattern(item_type.pattern());
        let anchored = format!(r"\A(?:{pattern})\z");
        let pattern = RegexBuilder::new(&anchored)
            .dot_matches_new_line(true)
            .build()
            .map_err(|error| {
                SchemaError::new(
                    "CIF_SCHEMA_PATTERN",
                    format!("invalid DDL2 type {:?}: {error}", item_type.code()),
                )
            })?;
        patterns.insert(item_type.code().to_ascii_lowercase(), pattern);
    }
    Ok(LoadedSchema {
        dictionary,
        patterns,
    })
}

fn normalize_ddl2_pattern(pattern: &str) -> String {
    // DDL2 dictionaries put `]` and `[` first in a character class as `[][`. Rust's
    // regex syntax requires both characters to be escaped explicitly.
    pattern.replace("[][", r"[\]\[")
}

#[cfg(test)]
#[allow(clippy::expect_used)]
mod tests {
    use super::{SchemaName, loaded_schema};

    #[test]
    fn embedded_schemas_match_their_locks() {
        let pdbx = loaded_schema(SchemaName::Pdbx).expect("PDBx artifact is valid");
        assert_eq!(pdbx.dictionary().metadata().version(), "5.416");
        assert_eq!(pdbx.dictionary().categories().len(), 607);
        assert_eq!(pdbx.dictionary().items().len(), 6801);
        let formula_weight = pdbx
            .dictionary()
            .item("_chem_comp.formula_weight")
            .expect("known PDBx item");
        assert!(
            formula_weight
                .ranges()
                .iter()
                .any(|range| { range.minimum() == Some("1.0") && range.maximum().is_none() })
        );

        let modelcif = loaded_schema(SchemaName::ModelCif).expect("ModelCIF artifact is valid");
        assert_eq!(modelcif.dictionary().metadata().version(), "1.4.9");
        assert_eq!(modelcif.dictionary().categories().len(), 619);
        assert_eq!(modelcif.dictionary().items().len(), 6722);
    }
}
