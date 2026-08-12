use std::collections::BTreeMap;
use std::error::Error;
use std::fmt::{self, Display, Formatter};

use super::dictionary::{
    CategoryDefinition, Dictionary, DictionaryMetadata, ItemDefinition, ItemRange, TypeDefinition,
};

const MAGIC: &[u8; 8] = b"NIBSCM01";
const FORMAT_VERSION: u32 = 1;
const MAX_RECORDS: usize = 10_000_000;

/// A malformed, unsupported, or excessively large compiled-schema artifact.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct SchemaArtifactError(String);

impl SchemaArtifactError {
    fn new(message: impl Into<String>) -> Self {
        Self(message.into())
    }

    /// Return the human-readable artifact failure.
    #[must_use]
    pub fn message(&self) -> &str {
        &self.0
    }
}

impl Display for SchemaArtifactError {
    fn fmt(&self, formatter: &mut Formatter<'_>) -> fmt::Result {
        formatter.write_str(&self.0)
    }
}

impl Error for SchemaArtifactError {}

/// Encode a compiled dictionary into Nibbler's deterministic binary schema format.
///
/// # Errors
///
/// Returns [`SchemaArtifactError`] if a collection or string cannot be represented by
/// the format's 32-bit length fields.
pub fn encode_dictionary(dictionary: &Dictionary) -> Result<Vec<u8>, SchemaArtifactError> {
    let mut writer = Writer::default();
    writer.bytes.extend_from_slice(MAGIC);
    writer.u32(FORMAT_VERSION);
    writer.string(&dictionary.metadata.schema_name)?;
    writer.string(&dictionary.metadata.dictionary_name)?;
    writer.string(&dictionary.metadata.version)?;
    writer.string(&dictionary.metadata.source_sha256)?;

    writer.count(dictionary.types.len())?;
    for definition in dictionary.types.values() {
        writer.string(&definition.code)?;
        writer.string(&definition.primitive)?;
        writer.string(&definition.pattern)?;
    }

    writer.count(dictionary.categories.len())?;
    for category in dictionary.categories.values() {
        writer.string(&category.id)?;
        writer.boolean(category.mandatory);
        writer.strings(&category.keys)?;
        writer.strings(&category.groups)?;
        writer.strings(&category.canonical_items)?;
    }

    writer.count(dictionary.items.len())?;
    for item in dictionary.items.values() {
        writer.string(&item.name)?;
        writer.string(&item.category)?;
        writer.boolean(item.mandatory);
        writer.string(&item.type_code)?;
        writer.strings(&item.enumerations)?;
        writer.count(item.ranges.len())?;
        for range in &item.ranges {
            writer.optional_string(range.minimum.as_deref())?;
            writer.optional_string(range.maximum.as_deref())?;
        }
        writer.strings(&item.aliases)?;
        writer.strings(&item.parents)?;
    }

    writer.count(dictionary.aliases.len())?;
    for (alias, canonical) in &dictionary.aliases {
        writer.string(alias)?;
        writer.optional_string(canonical.as_deref())?;
    }
    Ok(writer.bytes)
}

/// Decode one bounded Nibbler compiled-schema artifact.
///
/// # Errors
///
/// Returns [`SchemaArtifactError`] for a bad magic value, unsupported format version,
/// truncation, invalid UTF-8, duplicate map key, excessive count, or trailing bytes.
pub fn decode_dictionary(bytes: &[u8]) -> Result<Dictionary, SchemaArtifactError> {
    let mut reader = Reader::new(bytes);
    if reader.take(MAGIC.len())? != MAGIC {
        return Err(SchemaArtifactError::new("invalid compiled-schema magic"));
    }
    let version = reader.u32()?;
    if version != FORMAT_VERSION {
        return Err(SchemaArtifactError::new(format!(
            "unsupported compiled-schema format {version}"
        )));
    }
    let metadata = DictionaryMetadata {
        schema_name: reader.string()?,
        dictionary_name: reader.string()?,
        version: reader.string()?,
        source_sha256: reader.string()?,
    };

    let mut types = BTreeMap::new();
    for _ in 0..reader.count()? {
        let definition = TypeDefinition {
            code: reader.string()?,
            primitive: reader.string()?,
            pattern: reader.string()?,
        };
        insert_unique(&mut types, definition.code.to_ascii_lowercase(), definition)?;
    }

    let mut categories = BTreeMap::new();
    for _ in 0..reader.count()? {
        let category = CategoryDefinition {
            id: reader.string()?,
            mandatory: reader.boolean()?,
            keys: reader.strings()?,
            groups: reader.strings()?,
            canonical_items: reader.strings()?,
        };
        insert_unique(&mut categories, category.id.to_ascii_lowercase(), category)?;
    }

    let mut items = BTreeMap::new();
    for _ in 0..reader.count()? {
        let name = reader.string()?;
        let category = reader.string()?;
        let mandatory = reader.boolean()?;
        let type_code = reader.string()?;
        let enumerations = reader.strings()?;
        let mut ranges = Vec::new();
        for _ in 0..reader.count()? {
            ranges.push(ItemRange {
                minimum: reader.optional_string()?,
                maximum: reader.optional_string()?,
            });
        }
        let aliases = reader.strings()?;
        let parents = reader.strings()?;
        let item = ItemDefinition {
            name,
            category,
            mandatory,
            type_code,
            enumerations,
            ranges,
            aliases,
            parents,
        };
        insert_unique(&mut items, item.name.to_ascii_lowercase(), item)?;
    }

    let mut aliases = BTreeMap::new();
    for _ in 0..reader.count()? {
        let alias = reader.string()?;
        let canonical = reader.optional_string()?;
        insert_unique(&mut aliases, alias, canonical)?;
    }
    if !reader.is_finished() {
        return Err(SchemaArtifactError::new(
            "compiled schema contains trailing bytes",
        ));
    }
    Ok(Dictionary {
        metadata,
        types,
        categories,
        items,
        aliases,
    })
}

fn insert_unique<T>(
    map: &mut BTreeMap<String, T>,
    key: String,
    value: T,
) -> Result<(), SchemaArtifactError> {
    if map.insert(key.clone(), value).is_some() {
        return Err(SchemaArtifactError::new(format!(
            "duplicate compiled-schema key {key:?}"
        )));
    }
    Ok(())
}

#[derive(Default)]
struct Writer {
    bytes: Vec<u8>,
}

impl Writer {
    fn u32(&mut self, value: u32) {
        self.bytes.extend_from_slice(&value.to_le_bytes());
    }

    fn count(&mut self, value: usize) -> Result<(), SchemaArtifactError> {
        let value = u32::try_from(value)
            .map_err(|_| SchemaArtifactError::new("compiled-schema collection is too large"))?;
        self.u32(value);
        Ok(())
    }

    fn boolean(&mut self, value: bool) {
        self.bytes.push(u8::from(value));
    }

    fn string(&mut self, value: &str) -> Result<(), SchemaArtifactError> {
        self.count(value.len())?;
        self.bytes.extend_from_slice(value.as_bytes());
        Ok(())
    }

    fn strings(&mut self, values: &[String]) -> Result<(), SchemaArtifactError> {
        self.count(values.len())?;
        for value in values {
            self.string(value)?;
        }
        Ok(())
    }

    fn optional_string(&mut self, value: Option<&str>) -> Result<(), SchemaArtifactError> {
        self.boolean(value.is_some());
        if let Some(value) = value {
            self.string(value)?;
        }
        Ok(())
    }
}

struct Reader<'a> {
    bytes: &'a [u8],
    offset: usize,
}

impl<'a> Reader<'a> {
    const fn new(bytes: &'a [u8]) -> Self {
        Self { bytes, offset: 0 }
    }

    fn take(&mut self, count: usize) -> Result<&'a [u8], SchemaArtifactError> {
        let end = self
            .offset
            .checked_add(count)
            .ok_or_else(|| SchemaArtifactError::new("compiled-schema offset overflow"))?;
        let output = self
            .bytes
            .get(self.offset..end)
            .ok_or_else(|| SchemaArtifactError::new("truncated compiled schema"))?;
        self.offset = end;
        Ok(output)
    }

    fn u32(&mut self) -> Result<u32, SchemaArtifactError> {
        let bytes = self.take(4)?;
        let array = <[u8; 4]>::try_from(bytes)
            .map_err(|_| SchemaArtifactError::new("invalid compiled-schema integer"))?;
        Ok(u32::from_le_bytes(array))
    }

    fn count(&mut self) -> Result<usize, SchemaArtifactError> {
        let count = usize::try_from(self.u32()?)
            .map_err(|_| SchemaArtifactError::new("compiled-schema count is too large"))?;
        if count > MAX_RECORDS {
            return Err(SchemaArtifactError::new(format!(
                "compiled-schema count exceeds {MAX_RECORDS}"
            )));
        }
        Ok(count)
    }

    fn boolean(&mut self) -> Result<bool, SchemaArtifactError> {
        match self.take(1)? {
            [0] => Ok(false),
            [1] => Ok(true),
            _ => Err(SchemaArtifactError::new("invalid compiled-schema boolean")),
        }
    }

    fn string(&mut self) -> Result<String, SchemaArtifactError> {
        let count = self.count()?;
        let bytes = self.take(count)?;
        let text = std::str::from_utf8(bytes)
            .map_err(|_| SchemaArtifactError::new("compiled schema contains invalid UTF-8"))?;
        Ok(text.to_owned())
    }

    fn strings(&mut self) -> Result<Vec<String>, SchemaArtifactError> {
        let count = self.count()?;
        let mut output = Vec::with_capacity(count);
        for _ in 0..count {
            output.push(self.string()?);
        }
        Ok(output)
    }

    fn optional_string(&mut self) -> Result<Option<String>, SchemaArtifactError> {
        self.boolean()?.then(|| self.string()).transpose()
    }

    const fn is_finished(&self) -> bool {
        self.offset == self.bytes.len()
    }
}

#[cfg(test)]
#[allow(clippy::expect_used)]
mod tests {
    use super::{decode_dictionary, encode_dictionary};
    use crate::cif::{compile_dictionary, parse};

    const DICTIONARY: &[u8] = br#"
data_test.dic
_dictionary.title test.dic
_dictionary.version 1.0
loop_
_item_type_list.code
_item_type_list.primitive_code
_item_type_list.construct
code char '[A-Za-z]+'
save_example
_category.id example
_category.mandatory_code no
_category_key.name '_example.id'
save_
save__example.id
_item.name '_example.id'
_item.category_id example
_item.mandatory_code yes
_item_type.code code
save_
"#;

    #[test]
    fn artifact_round_trip_is_exact_and_rejects_trailing_bytes() {
        let document = parse(DICTIONARY).expect("test dictionary syntax is valid");
        let dictionary =
            compile_dictionary(&document, "test", "abc").expect("test dictionary is consistent");
        let encoded = encode_dictionary(&dictionary).expect("test schema fits artifact limits");
        assert_eq!(decode_dictionary(&encoded), Ok(dictionary));

        let mut invalid = encoded;
        invalid.push(0);
        assert!(decode_dictionary(&invalid).is_err());
    }
}
