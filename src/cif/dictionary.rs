use std::collections::BTreeMap;
use std::error::Error;
use std::fmt::{self, Display, Formatter};

/// Metadata that binds a compiled schema to one immutable dictionary input.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct DictionaryMetadata {
    pub(crate) schema_name: String,
    pub(crate) dictionary_name: String,
    pub(crate) version: String,
    pub(crate) source_sha256: String,
}

impl DictionaryMetadata {
    /// Return the short schema selector, such as `pdbx` or `modelcif`.
    #[must_use]
    pub fn schema_name(&self) -> &str {
        &self.schema_name
    }

    /// Return the DDL2 dictionary title.
    #[must_use]
    pub fn dictionary_name(&self) -> &str {
        &self.dictionary_name
    }

    /// Return the pinned dictionary version.
    #[must_use]
    pub fn version(&self) -> &str {
        &self.version
    }

    /// Return the SHA-256 digest of the source dictionary.
    #[must_use]
    pub fn source_sha256(&self) -> &str {
        &self.source_sha256
    }
}

/// One named DDL2 value type.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct TypeDefinition {
    pub(crate) code: String,
    pub(crate) primitive: String,
    pub(crate) pattern: String,
}

impl TypeDefinition {
    /// Return the dictionary-local type code.
    #[must_use]
    pub fn code(&self) -> &str {
        &self.code
    }

    /// Return the DDL2 primitive code.
    #[must_use]
    pub fn primitive(&self) -> &str {
        &self.primitive
    }

    /// Return the unanchored DDL2 regular-expression construct.
    #[must_use]
    pub fn pattern(&self) -> &str {
        &self.pattern
    }
}

/// One optional inclusive numeric interval from `_item_range`.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ItemRange {
    pub(crate) minimum: Option<String>,
    pub(crate) maximum: Option<String>,
}

impl ItemRange {
    /// Return the inclusive lower bound, if present.
    #[must_use]
    pub fn minimum(&self) -> Option<&str> {
        self.minimum.as_deref()
    }

    /// Return the inclusive upper bound, if present.
    #[must_use]
    pub fn maximum(&self) -> Option<&str> {
        self.maximum.as_deref()
    }
}

/// One DDL2 item definition.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ItemDefinition {
    pub(crate) name: String,
    pub(crate) category: String,
    pub(crate) mandatory: bool,
    pub(crate) type_code: String,
    pub(crate) enumerations: Vec<String>,
    pub(crate) ranges: Vec<ItemRange>,
    pub(crate) aliases: Vec<String>,
    pub(crate) parents: Vec<String>,
}

impl ItemDefinition {
    /// Return the canonical item name including its leading underscore.
    #[must_use]
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Return the owning category ID without a leading underscore.
    #[must_use]
    pub fn category(&self) -> &str {
        &self.category
    }

    /// Return whether DDL2 marks the item mandatory when its category occurs.
    #[must_use]
    pub const fn mandatory(&self) -> bool {
        self.mandatory
    }

    /// Return the dictionary-local type code.
    #[must_use]
    pub fn type_code(&self) -> &str {
        &self.type_code
    }

    /// Return allowed enumeration values.
    #[must_use]
    pub fn enumerations(&self) -> &[String] {
        &self.enumerations
    }

    /// Return allowed numeric intervals.
    #[must_use]
    pub fn ranges(&self) -> &[ItemRange] {
        &self.ranges
    }

    /// Return historical aliases accepted for this item.
    #[must_use]
    pub fn aliases(&self) -> &[String] {
        &self.aliases
    }

    /// Return parent item names referenced by DDL2 links.
    #[must_use]
    pub fn parents(&self) -> &[String] {
        &self.parents
    }
}

/// One DDL2 category definition.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct CategoryDefinition {
    pub(crate) id: String,
    pub(crate) mandatory: bool,
    pub(crate) keys: Vec<String>,
    pub(crate) groups: Vec<String>,
    pub(crate) canonical_items: Vec<String>,
}

impl CategoryDefinition {
    /// Return the category ID without a leading underscore.
    #[must_use]
    pub fn id(&self) -> &str {
        &self.id
    }

    /// Return whether DDL2 marks the whole category mandatory.
    #[must_use]
    pub const fn mandatory(&self) -> bool {
        self.mandatory
    }

    /// Return the composite key item names.
    #[must_use]
    pub fn keys(&self) -> &[String] {
        &self.keys
    }

    /// Return category group memberships.
    #[must_use]
    pub fn groups(&self) -> &[String] {
        &self.groups
    }

    /// Return item names in dictionary encounter order.
    #[must_use]
    pub fn canonical_items(&self) -> &[String] {
        &self.canonical_items
    }
}

/// A compact logical DDL2 schema before or after artifact serialization.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Dictionary {
    pub(crate) metadata: DictionaryMetadata,
    pub(crate) types: BTreeMap<String, TypeDefinition>,
    pub(crate) categories: BTreeMap<String, CategoryDefinition>,
    pub(crate) items: BTreeMap<String, ItemDefinition>,
    pub(crate) aliases: BTreeMap<String, Option<String>>,
}

impl Dictionary {
    /// Return immutable schema-lock metadata.
    #[must_use]
    pub const fn metadata(&self) -> &DictionaryMetadata {
        &self.metadata
    }

    /// Find a type definition by ASCII case-insensitive code.
    #[must_use]
    pub fn value_type(&self, code: &str) -> Option<&TypeDefinition> {
        self.types.get(&fold(code))
    }

    /// Find a category by ASCII case-insensitive ID.
    #[must_use]
    pub fn category(&self, id: &str) -> Option<&CategoryDefinition> {
        self.categories.get(&fold(id.trim_start_matches('_')))
    }

    /// Find an item or alias by ASCII case-insensitive full name.
    #[must_use]
    pub fn item(&self, name: &str) -> Option<&ItemDefinition> {
        let key = fold(name);
        match self.aliases.get(&key) {
            Some(Some(canonical)) => self.items.get(canonical),
            Some(None) => None,
            None => self.items.get(&key),
        }
    }

    /// Iterate over all categories in stable lexical order.
    pub fn categories(&self) -> impl ExactSizeIterator<Item = &CategoryDefinition> {
        self.categories.values()
    }

    /// Iterate over all canonical items in stable lexical order.
    pub fn items(&self) -> impl ExactSizeIterator<Item = &ItemDefinition> {
        self.items.values()
    }
}

/// A malformed or internally inconsistent DDL2 dictionary.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct DictionaryError(String);

impl DictionaryError {
    pub(super) fn new(message: impl Into<String>) -> Self {
        Self(message.into())
    }

    /// Return the human-readable compiler failure.
    #[must_use]
    pub fn message(&self) -> &str {
        &self.0
    }
}

impl Display for DictionaryError {
    fn fmt(&self, formatter: &mut Formatter<'_>) -> fmt::Result {
        formatter.write_str(&self.0)
    }
}

impl Error for DictionaryError {}

pub(super) fn fold(value: &str) -> String {
    value.to_ascii_lowercase()
}
