use std::collections::BTreeMap;
use std::io::Cursor;

use serde::{Deserialize, Serialize};

use super::dictionary::{
    CategoryDefinition, Dictionary, DictionaryMetadata, ItemDefinition, TypeDefinition, fold,
};

const MAGIC: &[u8; 6] = b"NIBSCM";
const FORMAT_VERSION: u16 = 2;
const HEADER_SIZE: usize = MAGIC.len() + size_of::<u16>();
const MAX_ARTIFACT_BYTES: usize = 16 * 1024 * 1024;

#[derive(Serialize)]
struct Artifact<'a> {
    metadata: &'a DictionaryMetadata,
    types: Vec<&'a TypeDefinition>,
    categories: Vec<&'a CategoryDefinition>,
    items: Vec<&'a ItemDefinition>,
    aliases: Vec<(&'a str, Option<&'a str>)>,
}

#[derive(Deserialize)]
struct DecodedArtifact {
    metadata: DictionaryMetadata,
    types: Vec<TypeDefinition>,
    categories: Vec<CategoryDefinition>,
    items: Vec<ItemDefinition>,
    aliases: Vec<(String, Option<String>)>,
}

/// Encode a compiled dictionary into Nibbler's deterministic schema format.
///
/// # Errors
///
/// Returns an error if MessagePack cannot represent the dictionary.
pub fn encode_dictionary(dictionary: &Dictionary) -> Result<Vec<u8>, rmp_serde::encode::Error> {
    let body = rmp_serde::to_vec(&Artifact {
        metadata: &dictionary.metadata,
        types: dictionary.types.values().collect(),
        categories: dictionary.categories.values().collect(),
        items: dictionary.items.values().collect(),
        aliases: dictionary
            .aliases
            .iter()
            .map(|(alias, canonical)| (alias.as_str(), canonical.as_deref()))
            .collect(),
    })?;
    let mut artifact = Vec::with_capacity(HEADER_SIZE + body.len());
    artifact.extend_from_slice(MAGIC);
    artifact.extend_from_slice(&FORMAT_VERSION.to_le_bytes());
    artifact.extend_from_slice(&body);
    Ok(artifact)
}

pub(super) fn decode_dictionary(bytes: &[u8]) -> Result<Dictionary, String> {
    if bytes.len() > MAX_ARTIFACT_BYTES {
        return Err(format!(
            "compiled schema exceeds the {MAX_ARTIFACT_BYTES}-byte limit"
        ));
    }
    let header = bytes
        .get(..HEADER_SIZE)
        .ok_or_else(|| "truncated compiled-schema header".to_owned())?;
    if &header[..MAGIC.len()] != MAGIC {
        return Err("invalid compiled-schema magic".to_owned());
    }
    let version = u16::from_le_bytes([header[MAGIC.len()], header[MAGIC.len() + 1]]);
    if version != FORMAT_VERSION {
        return Err(format!("unsupported compiled-schema format {version}"));
    }

    let body = &bytes[HEADER_SIZE..];
    let mut decoder = rmp_serde::Deserializer::new(Cursor::new(body));
    let artifact = DecodedArtifact::deserialize(&mut decoder)
        .map_err(|error| format!("invalid compiled schema: {error}"))?;
    if decoder.position() != body.len() as u64 {
        return Err("compiled schema contains trailing bytes".to_owned());
    }
    Ok(Dictionary {
        metadata: artifact.metadata,
        types: collect_unique(artifact.types, |value| &value.code)?,
        categories: collect_unique(artifact.categories, |value| &value.id)?,
        items: collect_unique(artifact.items, |value| &value.name)?,
        aliases: collect_aliases(artifact.aliases)?,
    })
}

fn collect_unique<T>(
    values: Vec<T>,
    key: impl Fn(&T) -> &str,
) -> Result<BTreeMap<String, T>, String> {
    let mut output = BTreeMap::new();
    for value in values {
        let key = fold(key(&value));
        if output.insert(key.clone(), value).is_some() {
            return Err(format!("duplicate compiled-schema key {key:?}"));
        }
    }
    Ok(output)
}

fn collect_aliases(
    values: Vec<(String, Option<String>)>,
) -> Result<BTreeMap<String, Option<String>>, String> {
    let mut output = BTreeMap::new();
    for (alias, canonical) in values {
        let alias = fold(&alias);
        if output.insert(alias.clone(), canonical).is_some() {
            return Err(format!("duplicate compiled-schema key {alias:?}"));
        }
    }
    Ok(output)
}

#[cfg(test)]
#[allow(clippy::expect_used)]
mod tests {
    use super::{
        FORMAT_VERSION, HEADER_SIZE, MAX_ARTIFACT_BYTES, decode_dictionary, encode_dictionary,
    };
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
    fn artifact_round_trip_is_exact_and_bounded() {
        let document = parse(DICTIONARY).expect("test dictionary syntax is valid");
        let dictionary =
            compile_dictionary(&document, "test", "abc").expect("test dictionary is consistent");
        let encoded = encode_dictionary(&dictionary).expect("test dictionary is serializable");
        assert_eq!(decode_dictionary(&encoded), Ok(dictionary));

        let mut invalid = encoded.clone();
        invalid.push(0);
        assert!(decode_dictionary(&invalid).is_err());

        invalid = encoded;
        invalid[HEADER_SIZE - 2..HEADER_SIZE].copy_from_slice(&(FORMAT_VERSION + 1).to_le_bytes());
        assert!(decode_dictionary(&invalid).is_err());
        assert!(decode_dictionary(&vec![0; MAX_ARTIFACT_BYTES + 1]).is_err());
    }
}
