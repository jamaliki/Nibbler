//! Explicit row selection preserves source identifiers and matching authored links.

use super::SemanticError;
use super::fields::{category_rows, optional_integer, optional_text};
use super::model::AtomSite;
use crate::cif::{
    BlockCategory, CategoryView, CifBlock, CifDocument, CifEntry, CifLoop, in_category,
};
use std::collections::{HashMap, HashSet};

type SiteIndex<'a> = HashMap<(&'a str, &'a str), Vec<(usize, &'a AtomSite)>>;

/// Copy selected atom_site rows without renumbering deposited atom identifiers.
///
/// Keep authored links to unresolved atoms, but remove links whose observed
/// endpoints belong only to discarded models or alternatives. Other categories
/// remain shared. Selecting every source row returns a shared document clone.
///
/// # Errors
///
/// Rejects invalid row indices, repeated atom tables, and mixed-category
/// coordinate loops, which cannot represent this row view.
pub fn selected_atom_document(
    source: &CifDocument,
    selected: &[usize],
) -> Result<CifDocument, SemanticError> {
    let [block] = source.blocks() else {
        return Err(error("Selection requires one source block"));
    };
    check_tables(block.entries())?;
    let table = BlockCategory::find(block.entries(), "atom_site")
        .ok_or_else(|| error("Source has no atom_site table"))?;
    if selected.iter().any(|&r| r >= table.row_count()) {
        return Err(error("Selected atom row is out of range"));
    }
    if selected.iter().copied().eq(0..table.row_count()) {
        return Ok(source.clone());
    }
    let rows: HashSet<_> = selected.iter().copied().collect();
    let connection_rows = selected_connections(CategoryView::new(block.entries()), &rows)?;
    let mut ids = HashSet::new();
    let has_anisotrop = BlockCategory::find(block.entries(), "atom_site_anisotrop").is_some();
    if has_anisotrop && let Some(col) = table.column("id") {
        for &r in selected {
            if let Some(id) = table.text(r, col) {
                ids.insert(id.into_owned());
            }
        }
    }
    let mut entries = Vec::with_capacity(block.entries().len());
    let mut connection_offset = 0;
    let anisotrop_keep = BlockCategory::find(block.entries(), "atom_site_anisotrop")
        .and_then(|table| table.column("id").and_then(|col| table.text(0, col)))
        .is_some_and(|id| ids.contains(id.as_ref()));
    for entry in block.entries() {
        if let CifEntry::Item(item) = entry {
            let keep = if in_category(item.tag(), "atom_site") {
                rows.contains(&0)
            } else if in_category(item.tag(), "atom_site_anisotrop") {
                anisotrop_keep
            } else if in_category(item.tag(), "struct_conn") {
                connection_rows.contains(&0)
            } else {
                true
            };
            if keep {
                entries.push(entry.clone());
            }
            continue;
        }
        if let CifEntry::Loop(values) = entry {
            let atom = values
                .tags()
                .iter()
                .all(|tag| in_category(tag, "atom_site"));
            let anisotrop = values
                .tags()
                .iter()
                .all(|tag| in_category(tag, "atom_site_anisotrop"));
            let connection = values
                .tags()
                .iter()
                .all(|tag| in_category(tag, "struct_conn"));
            if atom || anisotrop || connection {
                let id_col = values
                    .tags()
                    .iter()
                    .position(|tag| tag.eq_ignore_ascii_case("_atom_site_anisotrop.id"));
                let mut cells = Vec::new();
                for r in 0..values.row_count() {
                    let keep = if atom {
                        rows.contains(&r)
                    } else if connection {
                        connection_rows.contains(&(r + connection_offset))
                    } else {
                        id_col
                            .and_then(|c| values.value(r, c))
                            .and_then(|v| v.as_text())
                            .is_some_and(|s| ids.contains(s))
                    };
                    if keep {
                        for c in 0..values.tags().len() {
                            let value = values
                                .value(r, c)
                                .ok_or_else(|| error("Invalid source loop"))?;
                            cells.push(value.to_owned());
                        }
                    }
                }
                if connection {
                    connection_offset += values.row_count();
                }
                if !cells.is_empty() {
                    entries.push(CifEntry::Loop(CifLoop::new(values.tags().to_vec(), cells)));
                }
                continue;
            }
        }
        entries.push(entry.clone());
    }
    Ok(CifDocument::new(vec![CifBlock::data(
        block.code().unwrap_or_default().to_owned(),
        entries,
    )]))
}

pub(super) fn check_tables(entries: &[CifEntry]) -> Result<(), SemanticError> {
    let mut atom_loops = 0;
    let mut scalar = [false; 3];
    let mut loops = [0; 3];
    for entry in entries {
        for (index, category) in ["atom_site", "atom_site_anisotrop", "struct_conn"]
            .into_iter()
            .enumerate()
        {
            match entry {
                CifEntry::Item(item) if in_category(item.tag(), category) => {
                    scalar[index] = true;
                }
                CifEntry::Loop(values)
                    if values.tags().iter().any(|tag| in_category(tag, category)) =>
                {
                    if !values.tags().iter().all(|tag| in_category(tag, category)) {
                        return Err(error(
                            "Selection cannot filter a mixed-category coordinate loop",
                        ));
                    }
                    if category == "atom_site" {
                        atom_loops += 1;
                    }
                    loops[index] += 1;
                }
                _ => {}
            }
        }
    }
    if scalar
        .iter()
        .zip(loops)
        .any(|(&items, count)| items && count > 0)
    {
        return Err(error(
            "Selection cannot mix scalar and loop rows in one coordinate category",
        ));
    }
    if atom_loops + usize::from(scalar[0]) != 1 {
        return Err(error("Selection requires one atom_site table"));
    }
    Ok(())
}

fn selected_connections(
    categories: CategoryView<'_>,
    selected: &HashSet<usize>,
) -> Result<HashSet<usize>, SemanticError> {
    let fields = [
        "ptnr1_label_asym_id",
        "ptnr1_label_atom_id",
        "ptnr1_label_comp_id",
        "ptnr1_label_seq_id",
        "ptnr1_auth_seq_id",
        "pdbx_ptnr1_label_alt_id",
        "pdbx_ptnr1_PDB_ins_code",
        "ptnr2_label_asym_id",
        "ptnr2_label_atom_id",
        "ptnr2_label_comp_id",
        "ptnr2_label_seq_id",
        "ptnr2_auth_seq_id",
        "pdbx_ptnr2_label_alt_id",
        "pdbx_ptnr2_PDB_ins_code",
    ];
    let mut connections = category_rows(categories, "struct_conn", fields).peekable();
    if connections.peek().is_none() {
        return Ok(HashSet::new());
    }
    let atoms = super::source::atom_sites(categories)?;
    let mut index: SiteIndex<'_> = HashMap::new();
    for (row, atom) in atoms.iter().enumerate() {
        index
            .entry((&atom.label_asym_id, &atom.label_atom_id))
            .or_default()
            .push((row, atom));
    }
    let mut retained = HashSet::new();
    for row in connections {
        let mut keep = true;
        for offset in [0, 7] {
            let values = [0, 1, 2, 4, 5, 6].map(|i| optional_text(&row, offset + i));
            let seq = optional_integer(&row, "struct_conn", fields[offset + 3], offset + 3)?;
            let key = (
                values[0].as_deref().unwrap_or_default(),
                values[1].as_deref().unwrap_or_default(),
            );
            let mut observed = false;
            let mut chosen = false;
            for &(r, atom) in index.get(&key).into_iter().flatten() {
                if values[2]
                    .as_deref()
                    .is_some_and(|s| s != atom.label_comp_id)
                    || seq.is_some_and(|s| Some(s) != atom.label_seq_id)
                    || values[3].as_deref().is_some_and(|s| s != atom.auth_seq_id)
                    || values[4]
                        .as_deref()
                        .is_some_and(|s| atom.label_alt_id.as_deref().is_some_and(|a| a != s))
                    || values[5]
                        .as_deref()
                        .is_some_and(|s| atom.insertion_code.as_deref() != Some(s))
                {
                    continue;
                }
                observed = true;
                chosen |= selected.contains(&r);
            }
            keep &= !observed || chosen;
        }
        if keep {
            retained.insert(row.row_index());
        }
    }
    Ok(retained)
}

fn error(message: &str) -> SemanticError {
    SemanticError::new("PDBX_SELECTION", message, Vec::new())
}
