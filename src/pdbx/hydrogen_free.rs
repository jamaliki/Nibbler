//! Remove deposited hydrogen sites while keeping every other category as authored.

use super::SemanticError;
use super::selected_document::check_tables_for;
use crate::cif::{
    BlockCategory, CifBlock, CifDocument, CifEntry, CifLoop, CifValueRef, in_category,
};
use std::borrow::Cow;
use std::collections::HashSet;

/// Element symbols that name a hydrogen isotope in `_atom_site.type_symbol`.
const HYDROGEN_TYPES: [&str; 3] = ["H", "D", "T"];

/// An atom identifier, keeping `?` and `.` distinct from present text.
#[derive(Debug, Eq, Hash, PartialEq)]
enum AtomId<'a> {
    Text(Cow<'a, str>),
    Unknown,
    NotApplicable,
}

impl<'a> AtomId<'a> {
    fn of(value: Option<CifValueRef<'a>>) -> Option<Self> {
        Some(match value? {
            CifValueRef::Unknown => Self::Unknown,
            CifValueRef::NotApplicable => Self::NotApplicable,
            CifValueRef::Text(text) => Self::Text(Cow::Borrowed(text.as_str())),
            CifValueRef::Integer(number, original) => Self::Text(original.map_or_else(
                || Cow::Owned(number.to_string()),
                |text| Cow::Borrowed(text.as_str()),
            )),
            CifValueRef::Float(number, _, original) => Self::Text(original.map_or_else(
                || Cow::Owned(number.to_string()),
                |text| Cow::Borrowed(text.as_str()),
            )),
        })
    }
}

/// Copy a document without its H, D and T `atom_site` rows.
///
/// A row is a hydrogen when its `type_symbol` is `H`, `D` or `T` in any case; a missing
/// symbol is kept. `atom_site_anisotrop` keeps the rows whose `id` names a retained
/// atom, so rows of removed or absent atoms go too. A category left with no rows is
/// dropped, and one that keeps rows keeps its loop or scalar form. Every other
/// category, `struct_conn` included, stays as authored, and a document that loses no
/// row is returned as a shared clone.
///
/// # Errors
///
/// Rejects documents without exactly one block, without one `atom_site` table with a
/// `type_symbol` column, and with mixed-category coordinate loops.
pub fn hydrogen_free_document(source: &CifDocument) -> Result<CifDocument, SemanticError> {
    let [block] = source.blocks() else {
        return Err(error("Hydrogen removal requires one source block"));
    };
    check_tables_for(block.entries(), "Hydrogen removal", error)?;
    let atoms = BlockCategory::find(block.entries(), "atom_site")
        .ok_or_else(|| error("Source has no atom_site table"))?;
    let symbol = atoms
        .column("type_symbol")
        .ok_or_else(|| error("Source atom_site has no type_symbol column"))?;
    let kept: Vec<bool> = (0..atoms.row_count())
        .map(|row| {
            atoms.text(row, symbol).is_none_or(|text| {
                !HYDROGEN_TYPES
                    .iter()
                    .any(|hydrogen| text.eq_ignore_ascii_case(hydrogen))
            })
        })
        .collect();
    let id = atoms.column("id");
    let ids: HashSet<AtomId<'_>> = (0..atoms.row_count())
        .filter(|&row| kept[row])
        .filter_map(|row| AtomId::of(id.and_then(|column| atoms.value(row, column))))
        .collect();
    let anisotrop_kept =
        |value: Option<CifValueRef<'_>>| AtomId::of(value).is_some_and(|id| ids.contains(&id));
    let anisotrop_id = |tags: &[String]| {
        tags.iter()
            .position(|tag| tag.eq_ignore_ascii_case("_atom_site_anisotrop.id"))
    };
    // Each anisotropic row is judged on its own id, in every loop or scalar form it takes.
    let scalar_anisotrop = block.entries().iter().find_map(|entry| match entry {
        CifEntry::Item(item) if item.tag().eq_ignore_ascii_case("_atom_site_anisotrop.id") => {
            Some(item.value().as_ref())
        }
        _ => None,
    });
    let scalar_anisotrop_kept = anisotrop_kept(scalar_anisotrop);
    let rows_of = |values: &CifLoop| -> Option<Vec<bool>> {
        let tag = values.tags().first()?;
        if in_category(tag, "atom_site") {
            Some(kept.clone())
        } else if in_category(tag, "atom_site_anisotrop") {
            let column = anisotrop_id(values.tags());
            Some(
                (0..values.row_count())
                    .map(|row| anisotrop_kept(column.and_then(|c| values.value(row, c))))
                    .collect(),
            )
        } else {
            None
        }
    };
    let item_kept = |tag: &str| {
        if in_category(tag, "atom_site") {
            kept[0]
        } else if in_category(tag, "atom_site_anisotrop") {
            scalar_anisotrop_kept
        } else {
            true
        }
    };
    let changed = block.entries().iter().any(|entry| match entry {
        CifEntry::Item(item) => !item_kept(item.tag()),
        CifEntry::Loop(values) => rows_of(values).is_some_and(|rows| rows.contains(&false)),
        CifEntry::Frame(_) => false,
    });
    if !changed {
        return Ok(source.clone());
    }
    let mut entries = Vec::with_capacity(block.entries().len());
    for entry in block.entries() {
        match entry {
            CifEntry::Item(item) => {
                if item_kept(item.tag()) {
                    entries.push(entry.clone());
                }
            }
            CifEntry::Loop(values) => {
                let Some(rows) = rows_of(values) else {
                    entries.push(entry.clone());
                    continue;
                };
                let mut cells = Vec::new();
                for row in (0..values.row_count()).filter(|&row| rows[row]) {
                    for column in 0..values.column_count() {
                        let value = values
                            .value(row, column)
                            .ok_or_else(|| error("Invalid source loop"))?;
                        cells.push(value.to_owned());
                    }
                }
                if !cells.is_empty() {
                    entries.push(CifEntry::Loop(CifLoop::new(values.tags().to_vec(), cells)));
                }
            }
            CifEntry::Frame(_) => entries.push(entry.clone()),
        }
    }
    Ok(CifDocument::new(vec![CifBlock::data(
        block.code().unwrap_or_default().to_owned(),
        entries,
    )]))
}

fn error(message: &str) -> SemanticError {
    SemanticError::new("PDBX_HYDROGEN_REMOVAL", message, Vec::new())
}
