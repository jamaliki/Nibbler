//! Compilation of normalized dictionaries from parsed DDL2 documents.

use std::collections::{BTreeMap, btree_map::Entry};

use super::dictionary::{
    CategoryDefinition, Dictionary, DictionaryError, DictionaryMetadata, ItemDefinition, ItemRange,
    TypeDefinition, fold,
};
use super::dictionary_source::{
    item_records, loop_records, range_records, required_scalar, scalar, values,
};
use super::document::{CifDocument, CifEntry};

/// Compile one parsed DDL2 document into a normalized schema.
///
/// The source document must contain one data block with `_dictionary.title`,
/// `_dictionary.version`, an `_item_type_list` loop, and category/item save frames.
///
/// # Errors
///
/// Returns [`DictionaryError`] for missing metadata, malformed DDL2 records, duplicate
/// incompatible definitions, or references to undefined categories, items, or types.
pub fn compile_dictionary(
    document: &CifDocument,
    schema_name: &str,
    source_sha256: &str,
) -> Result<Dictionary, DictionaryError> {
    let [block] = document.blocks() else {
        return Err(DictionaryError::new(
            "a DDL2 dictionary must contain exactly one data block",
        ));
    };
    let dictionary_name = required_scalar(block.entries(), "_dictionary.title")?;
    let version = required_scalar(block.entries(), "_dictionary.version")?;
    let mut dictionary = Dictionary {
        metadata: DictionaryMetadata {
            schema_name: schema_name.to_owned(),
            dictionary_name: dictionary_name.to_owned(),
            version: version.to_owned(),
            source_sha256: source_sha256.to_owned(),
        },
        types: BTreeMap::new(),
        categories: BTreeMap::new(),
        items: BTreeMap::new(),
        aliases: BTreeMap::new(),
    };
    compile_types(block.entries(), &mut dictionary)?;

    let frames = block
        .entries()
        .iter()
        .filter_map(|entry| match entry {
            CifEntry::Frame(frame) => Some(frame),
            CifEntry::Item(_) | CifEntry::Loop(_) => None,
        })
        .collect::<Vec<_>>();
    for frame in &frames {
        compile_category(frame.entries(), &mut dictionary)
            .map_err(|error| DictionaryError::new(format!("save_{}: {error}", frame.code())))?;
        register_item_order(frame.entries(), &mut dictionary);
    }

    let mut pending_aliases = Vec::new();
    let mut pending_links = Vec::new();
    for has_type in [true, false] {
        for frame in frames
            .iter()
            .filter(|frame| scalar(frame.entries(), "_item_type.code").is_some() == has_type)
        {
            compile_items(
                frame.entries(),
                &mut dictionary,
                &mut pending_aliases,
                &mut pending_links,
            )
            .map_err(|error| DictionaryError::new(format!("save_{}: {error}", frame.code())))?;
        }
    }
    apply_aliases(&mut dictionary, pending_aliases)?;
    apply_links(&mut dictionary, pending_links)?;
    resolve_inherited_constraints(&mut dictionary);
    validate_compiled_dictionary(&dictionary)?;
    Ok(dictionary)
}

fn register_item_order(entries: &[CifEntry], dictionary: &mut Dictionary) {
    for record in item_records(entries) {
        let [name, category, _] = record.as_slice() else {
            continue;
        };
        let item_key = fold(name);
        let category = dictionary
            .categories
            .entry(fold(category))
            .or_insert_with(|| CategoryDefinition {
                id: (*category).to_owned(),
                mandatory: false,
                keys: Vec::new(),
                groups: Vec::new(),
                canonical_items: Vec::new(),
            });
        if !category.canonical_items.contains(&item_key) {
            category.canonical_items.push(item_key);
        }
    }
}

fn compile_types(entries: &[CifEntry], dictionary: &mut Dictionary) -> Result<(), DictionaryError> {
    for record in loop_records(
        entries,
        &[
            "_item_type_list.code",
            "_item_type_list.primitive_code",
            "_item_type_list.construct",
        ],
    ) {
        let [code, primitive, pattern] = record.as_slice() else {
            return Err(DictionaryError::new("invalid item type record width"));
        };
        let definition = TypeDefinition {
            code: (*code).to_owned(),
            primitive: (*primitive).to_owned(),
            pattern: (*pattern).to_owned(),
        };
        if dictionary.types.insert(fold(code), definition).is_some() {
            return Err(DictionaryError::new(format!(
                "duplicate DDL2 type code {code:?}"
            )));
        }
    }
    if dictionary.types.is_empty() {
        return Err(DictionaryError::new(
            "dictionary has no item type definitions",
        ));
    }
    Ok(())
}

fn compile_category(
    entries: &[CifEntry],
    dictionary: &mut Dictionary,
) -> Result<(), DictionaryError> {
    let Some(id) = scalar(entries, "_category.id") else {
        return Ok(());
    };
    let key = fold(id);
    let definition = CategoryDefinition {
        id: id.to_owned(),
        mandatory: scalar(entries, "_category.mandatory_code").is_some_and(is_mandatory),
        keys: values(entries, "_category_key.name")
            .into_iter()
            .map(fold)
            .collect(),
        groups: values(entries, "_category_group.id")
            .into_iter()
            .map(str::to_owned)
            .collect(),
        canonical_items: Vec::new(),
    };
    match dictionary.categories.entry(key) {
        Entry::Vacant(entry) => {
            entry.insert(definition);
        }
        Entry::Occupied(mut entry) => {
            let existing = entry.get_mut();
            if !existing.id.eq_ignore_ascii_case(id) {
                return Err(DictionaryError::new(format!(
                    "incompatible duplicate DDL2 category {id:?}"
                )));
            }
            existing.mandatory |= definition.mandatory;
            merge_unique(&mut existing.keys, definition.keys);
            merge_unique(&mut existing.groups, definition.groups);
        }
    }
    Ok(())
}

fn compile_items(
    entries: &[CifEntry],
    dictionary: &mut Dictionary,
    pending_aliases: &mut Vec<(String, String)>,
    pending_links: &mut Vec<(String, String)>,
) -> Result<(), DictionaryError> {
    let item_records = item_records(entries);
    if item_records.is_empty() {
        return Ok(());
    }
    let type_code = scalar(entries, "_item_type.code");
    let enumerations = values(entries, "_item_enumeration.value")
        .into_iter()
        .map(str::to_owned)
        .collect::<Vec<_>>();
    let ranges = range_records(entries);
    let aliases = values(entries, "_item_aliases.alias_name")
        .into_iter()
        .map(str::to_owned)
        .collect::<Vec<_>>();

    for record in &item_records {
        let [name, category, mandatory] = record.as_slice() else {
            return Err(DictionaryError::new("invalid item record width"));
        };
        let key = fold(name);
        let category_key = fold(category);
        match dictionary.items.entry(key.clone()) {
            Entry::Vacant(entry) => {
                let item = ItemDefinition {
                    name: (*name).to_owned(),
                    category: (*category).to_owned(),
                    mandatory: is_mandatory(mandatory),
                    type_code: type_code.unwrap_or_default().to_owned(),
                    enumerations: enumerations.clone(),
                    ranges: ranges.clone(),
                    aliases: Vec::new(),
                    parents: Vec::new(),
                };
                entry.insert(item);
                dictionary
                    .categories
                    .entry(category_key)
                    .or_insert_with(|| CategoryDefinition {
                        id: (*category).to_owned(),
                        mandatory: false,
                        keys: Vec::new(),
                        groups: Vec::new(),
                        canonical_items: vec![key.clone()],
                    });
            }
            Entry::Occupied(mut entry) => merge_item_refinement(
                entry.get_mut(),
                category,
                mandatory,
                type_code,
                &enumerations,
                &ranges,
            )?,
        }
    }

    if let Some(primary) = item_records.first().and_then(|record| record.first()) {
        pending_aliases.extend(
            aliases
                .into_iter()
                .map(|alias| (alias, (*primary).to_owned())),
        );
    }
    pending_links.extend(
        loop_records(
            entries,
            &["_item_linked.child_name", "_item_linked.parent_name"],
        )
        .into_iter()
        .filter_map(|record| match record.as_slice() {
            [child, parent] => Some(((*child).to_owned(), (*parent).to_owned())),
            _ => None,
        }),
    );
    Ok(())
}

fn merge_item_refinement(
    existing: &mut ItemDefinition,
    category: &str,
    mandatory: &str,
    type_code: Option<&str>,
    enumerations: &[String],
    ranges: &[ItemRange],
) -> Result<(), DictionaryError> {
    let type_conflicts =
        type_code.is_some_and(|code| !existing.type_code.eq_ignore_ascii_case(code));
    let complete_definition_conflicts = type_code.is_some()
        && ((!existing.enumerations.is_empty()
            && !enumerations.is_empty()
            && existing.enumerations != enumerations)
            || (!existing.ranges.is_empty() && !ranges.is_empty() && existing.ranges != ranges));
    if !existing.category.eq_ignore_ascii_case(category)
        || type_conflicts
        || complete_definition_conflicts
    {
        return Err(DictionaryError::new(format!(
            "incompatible duplicate item definition {:?}",
            existing.name
        )));
    }
    existing.mandatory |= is_mandatory(mandatory);
    if !enumerations.is_empty() && (type_code.is_none() || existing.enumerations.is_empty()) {
        existing.enumerations = enumerations.to_vec();
    }
    if !ranges.is_empty() && (type_code.is_none() || existing.ranges.is_empty()) {
        existing.ranges = ranges.to_vec();
    }
    Ok(())
}

fn apply_aliases(
    dictionary: &mut Dictionary,
    aliases: Vec<(String, String)>,
) -> Result<(), DictionaryError> {
    for (alias, canonical) in aliases {
        let canonical_key = fold(&canonical);
        let Some(item) = dictionary.items.get_mut(&canonical_key) else {
            return Err(DictionaryError::new(format!(
                "alias {alias:?} references undefined item {canonical:?}"
            )));
        };
        let alias_key = fold(&alias);
        match dictionary.aliases.entry(alias_key) {
            Entry::Vacant(entry) => {
                entry.insert(Some(canonical_key.clone()));
            }
            Entry::Occupied(mut entry) => {
                if entry
                    .get()
                    .as_ref()
                    .is_some_and(|previous| previous != &canonical_key)
                {
                    entry.insert(None);
                }
            }
        }
        if !item
            .aliases
            .iter()
            .any(|existing| existing.eq_ignore_ascii_case(&alias))
        {
            item.aliases.push(alias);
        }
    }
    Ok(())
}

fn apply_links(
    dictionary: &mut Dictionary,
    links: Vec<(String, String)>,
) -> Result<(), DictionaryError> {
    for (child, parent) in links {
        let parent_key = fold(&parent);
        if !dictionary.items.contains_key(&parent_key) {
            return Err(DictionaryError::new(format!(
                "item link references undefined parent {parent:?}"
            )));
        }
        let child_key = fold(&child);
        let Some(item) = dictionary.items.get_mut(&child_key) else {
            return Err(DictionaryError::new(format!(
                "item link references undefined child {child:?}"
            )));
        };
        if !item.parents.contains(&parent_key) {
            item.parents.push(parent_key);
        }
    }
    Ok(())
}

fn resolve_inherited_constraints(dictionary: &mut Dictionary) {
    loop {
        let inherited = dictionary
            .items
            .iter()
            .filter(|(_, item)| item.type_code.is_empty())
            .filter_map(|(key, item)| {
                item.parents.iter().find_map(|parent| {
                    dictionary.items.get(parent).and_then(|definition| {
                        (!definition.type_code.is_empty()).then(|| {
                            (
                                key.clone(),
                                definition.type_code.clone(),
                                definition.enumerations.clone(),
                                definition.ranges.clone(),
                            )
                        })
                    })
                })
            })
            .collect::<Vec<_>>();
        if inherited.is_empty() {
            break;
        }
        for (key, type_code, enumerations, ranges) in inherited {
            let Some(item) = dictionary.items.get_mut(&key) else {
                continue;
            };
            item.type_code = type_code;
            if item.enumerations.is_empty() {
                item.enumerations = enumerations;
            }
            if item.ranges.is_empty() {
                item.ranges = ranges;
            }
        }
    }
    if dictionary.types.contains_key("any") {
        for item in dictionary
            .items
            .values_mut()
            .filter(|item| item.type_code.is_empty())
        {
            item.type_code = "any".to_owned();
        }
    }
}

fn validate_compiled_dictionary(dictionary: &Dictionary) -> Result<(), DictionaryError> {
    for item in dictionary.items.values() {
        if !dictionary.types.contains_key(&fold(&item.type_code)) {
            return Err(DictionaryError::new(format!(
                "item {:?} references undefined type {:?}",
                item.name, item.type_code
            )));
        }
        if !dictionary.categories.contains_key(&fold(&item.category)) {
            return Err(DictionaryError::new(format!(
                "item {:?} references undefined category {:?}",
                item.name, item.category
            )));
        }
    }
    for category in dictionary.categories.values() {
        for key in &category.keys {
            if !dictionary.items.contains_key(key) {
                return Err(DictionaryError::new(format!(
                    "category {:?} references undefined key {key:?}",
                    category.id
                )));
            }
        }
    }
    Ok(())
}

fn is_mandatory(value: &str) -> bool {
    value.eq_ignore_ascii_case("yes")
}

fn merge_unique(target: &mut Vec<String>, values: Vec<String>) {
    for value in values {
        if !target
            .iter()
            .any(|existing| existing.eq_ignore_ascii_case(&value))
        {
            target.push(value);
        }
    }
}
