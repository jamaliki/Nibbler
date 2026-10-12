use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex, OnceLock};

use reduce3::monlib::MonLib;

use crate::cif::{CifBlock, CifDocument, CifEntry, write_canonical};

type Overlays = HashMap<(PathBuf, Vec<String>), Arc<MonLib>>;
static OVERLAYS: OnceLock<Mutex<Overlays>> = OnceLock::new();

pub(super) fn overlay(
    base: Arc<MonLib>,
    documents: &[Arc<CifDocument>],
) -> Result<Arc<MonLib>, String> {
    if documents.is_empty() {
        return Ok(base);
    }
    let texts: Vec<_> = documents
        .iter()
        .map(|d| dictionary_text(d))
        .collect::<Result<_, _>>()?;
    let key = (base.root.clone(), texts.clone());
    let mut cache = OVERLAYS
        .get_or_init(|| Mutex::new(HashMap::new()))
        .lock()
        .map_err(|_| "Dictionary cache is poisoned")?;
    if let Some(library) = cache.get(&key) {
        return Ok(Arc::clone(library));
    }
    let inputs: Vec<_> = texts
        .iter()
        .map(|text| (text.as_str(), Path::new("<validation-dictionary>")))
        .collect();
    let library = Arc::new(base.with_validation_cif(&inputs)?);
    cache.insert(key, Arc::clone(&library));
    Ok(library)
}

pub(super) fn dictionary_text(document: &CifDocument) -> Result<String, String> {
    if document.blocks().iter().any(|block| block.entries().iter().any(|entry| {
        matches!(entry, CifEntry::Loop(table) if table.tags().iter().any(|tag| dictionary_tag(tag)) && !table.tags().iter().all(|tag| dictionary_tag(tag)))
    })) {
        return Err("Restraint and coordinate tags cannot share a dictionary loop".into());
    }
    if document.blocks().iter().any(|block| {
        block
            .entries()
            .iter()
            .any(|entry| matches!(entry, CifEntry::Frame(_)))
    }) {
        return Err("Save-frame restraint dictionaries need reference interpretation".into());
    }
    let blocks: Vec<_> = document
        .blocks()
        .iter()
        .filter_map(|block| {
            let entries: Vec<_> = block
                .entries()
                .iter()
                .filter(|entry| match entry {
                    CifEntry::Item(item) => dictionary_tag(item.tag()),
                    CifEntry::Loop(table) => table.tags().iter().all(|tag| dictionary_tag(tag)),
                    CifEntry::Frame(_) => false,
                })
                .cloned()
                .collect();
            (!entries.is_empty())
                .then(|| CifBlock::data(block.code().unwrap_or("dictionary").to_owned(), entries))
        })
        .collect();
    if blocks.is_empty() {
        return Err("The supplied dictionary has no chemical definitions".into());
    }
    write_canonical(&CifDocument::new(blocks)).map_err(|error| error.to_string())
}

fn dictionary_tag(tag: &str) -> bool {
    let tag = tag.to_ascii_lowercase();
    tag.starts_with("_chem_") || tag.starts_with("_lib_")
}
