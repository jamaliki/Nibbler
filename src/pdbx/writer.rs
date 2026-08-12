//! Deterministic PDBx profile ordering over the shared CIF serializer.

use crate::cif::{BlockKind, CifBlock, CifDocument, CifEntry, CifLoop, CifValue, split_tag};

use super::PdbxModel;
use super::model::{ComponentAtom, ComponentBond, ComponentDefinition, ComponentResolution};

const CATEGORY_ORDER: &[&str] = &[
    "entry",
    "audit_conform",
    "entity",
    "entity_poly",
    "entity_poly_seq",
    "pdbx_entity_nonpoly",
    "pdbx_entity_branch",
    "pdbx_entity_branch_list",
    "pdbx_entity_branch_link",
    "chem_comp",
    "chem_comp_atom",
    "chem_comp_bond",
    "struct_asym",
    "pdbx_poly_seq_scheme",
    "pdbx_nonpoly_scheme",
    "pdbx_branch_scheme",
    "atom_type",
    "atom_site",
    "struct_conn_type",
    "struct_conn",
];

const ATOM_SITE_ORDER: &[&str] = &[
    "group_pdb",
    "id",
    "type_symbol",
    "label_atom_id",
    "label_alt_id",
    "label_comp_id",
    "label_asym_id",
    "label_entity_id",
    "label_seq_id",
    "pdbx_pdb_ins_code",
    "cartn_x",
    "cartn_y",
    "cartn_z",
    "occupancy",
    "b_iso_or_equiv",
    "pdbx_formal_charge",
    "auth_seq_id",
    "auth_comp_id",
    "auth_asym_id",
    "auth_atom_id",
    "pdbx_pdb_model_num",
];

/// Construct a canonical-order document without duplicating CIF formatting logic.
#[must_use]
pub fn canonical_document(model: &PdbxModel) -> CifDocument {
    canonical_document_ordered(model, CATEGORY_ORDER)
}

pub(crate) fn canonical_document_ordered(
    model: &PdbxModel,
    category_order: &[&str],
) -> CifDocument {
    let blocks = model
        .source_document()
        .blocks()
        .iter()
        .map(|block| reorder_block(block, model, category_order))
        .collect();
    CifDocument::new(blocks)
}

fn reorder_block(block: &CifBlock, model: &PdbxModel, category_order: &[&str]) -> CifBlock {
    let mut source_entries = block.entries().to_vec();
    append_resolved_components(&mut source_entries, model);
    let mut entries = source_entries
        .iter()
        .enumerate()
        .map(|(source_index, entry)| {
            let entry = reorder_columns(entry);
            (category_rank(&entry, category_order), source_index, entry)
        })
        .collect::<Vec<_>>();
    entries.sort_by_key(|(rank, source_index, _)| (*rank, *source_index));
    let entries = entries.into_iter().map(|(_, _, entry)| entry).collect();
    match block.kind() {
        BlockKind::Data => CifBlock::data(block.code().unwrap_or_default().to_owned(), entries),
        BlockKind::Global => CifBlock::global(entries),
    }
}

fn append_resolved_components(entries: &mut Vec<CifEntry>, model: &PdbxModel) {
    let missing = model
        .components
        .iter()
        .filter(|component| {
            component.resolution != ComponentResolution::Embedded
                && component.resolution != ComponentResolution::Unresolved
        })
        .collect::<Vec<_>>();
    if missing.is_empty() {
        return;
    }
    if !merge_category_rows(entries, "chem_comp", &missing, component_value) {
        entries.push(component_loop(&missing));
    }
    let atoms = missing
        .iter()
        .flat_map(|component| component.atoms.iter().map(move |atom| (*component, atom)))
        .collect::<Vec<_>>();
    if !atoms.is_empty()
        && !merge_category_rows(entries, "chem_comp_atom", &atoms, component_atom_value)
    {
        entries.push(component_atom_loop(&atoms));
    }
    let bonds = missing
        .iter()
        .flat_map(|component| component.bonds.iter().map(move |bond| (*component, bond)))
        .collect::<Vec<_>>();
    if !bonds.is_empty()
        && !merge_category_rows(entries, "chem_comp_bond", &bonds, component_bond_value)
    {
        entries.push(component_bond_loop(&bonds));
    }
}

fn merge_category_rows<T>(
    entries: &mut Vec<CifEntry>,
    selected_category: &str,
    records: &[T],
    value: impl Fn(&T, &str) -> CifValue,
) -> bool {
    let loop_index = entries.iter().position(|entry| {
        matches!(entry, CifEntry::Loop(cif_loop) if cif_loop.tags().first().and_then(|tag| split_tag(tag).map(|(category, _)| category)).is_some_and(|name| name.eq_ignore_ascii_case(selected_category)))
    });
    if let Some(index) = loop_index {
        let CifEntry::Loop(cif_loop) = &entries[index] else {
            return false;
        };
        let tags = cif_loop.tags().to_vec();
        let mut values = cif_loop
            .values()
            .map(|value| value.to_owned())
            .collect::<Vec<_>>();
        for record in records {
            values.extend(
                tags.iter()
                    .map(|tag| value(record, split_tag(tag).map(|(_, item)| item).unwrap_or(""))),
            );
        }
        entries[index] = CifEntry::Loop(CifLoop::new(tags, values));
        return true;
    }
    merge_scalar_category(entries, selected_category, records, value)
}

fn merge_scalar_category<T>(
    entries: &mut Vec<CifEntry>,
    selected_category: &str,
    records: &[T],
    value: impl Fn(&T, &str) -> CifValue,
) -> bool {
    let indices = entries
        .iter()
        .enumerate()
        .filter_map(|(index, entry)| match entry {
            CifEntry::Item(item)
                if split_tag(item.tag())
                    .map(|(category, _)| category)
                    .is_some_and(|name| name.eq_ignore_ascii_case(selected_category)) =>
            {
                Some(index)
            }
            _ => None,
        })
        .collect::<Vec<_>>();
    let Some(&insert_at) = indices.first() else {
        return false;
    };
    let mut tags = Vec::with_capacity(indices.len());
    let mut values = Vec::with_capacity((records.len() + 1) * indices.len());
    for &index in &indices {
        let CifEntry::Item(item) = &entries[index] else {
            continue;
        };
        tags.push(item.tag().to_owned());
        values.push(item.value().clone());
    }
    for record in records {
        values.extend(
            tags.iter()
                .map(|tag| value(record, split_tag(tag).map(|(_, item)| item).unwrap_or(""))),
        );
    }
    for &index in indices.iter().rev() {
        entries.remove(index);
    }
    entries.insert(insert_at, CifEntry::Loop(CifLoop::new(tags, values)));
    true
}

fn component_loop(components: &[&ComponentDefinition]) -> CifEntry {
    let tags = [
        "_chem_comp.id",
        "_chem_comp.type",
        "_chem_comp.name",
        "_chem_comp.formula",
        "_chem_comp.formula_weight",
        "_chem_comp.pdbx_formal_charge",
    ];
    let values = components
        .iter()
        .flat_map(|component| {
            tags.map(|tag| {
                component_value(
                    component,
                    split_tag(tag).map(|(_, item)| item).unwrap_or(""),
                )
            })
        })
        .collect();
    CifEntry::Loop(CifLoop::new(tags.map(str::to_owned).to_vec(), values))
}

fn component_atom_loop(atoms: &[(&ComponentDefinition, &ComponentAtom)]) -> CifEntry {
    let tags = [
        "_chem_comp_atom.comp_id",
        "_chem_comp_atom.atom_id",
        "_chem_comp_atom.type_symbol",
        "_chem_comp_atom.charge",
    ];
    let values = atoms
        .iter()
        .flat_map(|record| {
            tags.map(|tag| {
                component_atom_value(record, split_tag(tag).map(|(_, item)| item).unwrap_or(""))
            })
        })
        .collect();
    CifEntry::Loop(CifLoop::new(tags.map(str::to_owned).to_vec(), values))
}

fn component_bond_loop(bonds: &[(&ComponentDefinition, &ComponentBond)]) -> CifEntry {
    let tags = [
        "_chem_comp_bond.comp_id",
        "_chem_comp_bond.atom_id_1",
        "_chem_comp_bond.atom_id_2",
        "_chem_comp_bond.value_order",
        "_chem_comp_bond.pdbx_aromatic_flag",
    ];
    let values = bonds
        .iter()
        .flat_map(|record| {
            tags.map(|tag| {
                component_bond_value(record, split_tag(tag).map(|(_, item)| item).unwrap_or(""))
            })
        })
        .collect();
    CifEntry::Loop(CifLoop::new(tags.map(str::to_owned).to_vec(), values))
}

fn component_value(component: &&ComponentDefinition, item: &str) -> CifValue {
    match item.to_ascii_lowercase().as_str() {
        "id" => CifValue::text(&component.id),
        "type" => optional_text(component.component_type.as_deref()),
        "name" => optional_text(component.name.as_deref()),
        "formula" => optional_text(component.formula.as_deref()),
        "formula_weight" => optional_float(component.formula_weight),
        "pdbx_formal_charge" => optional_integer(component.formal_charge.map(i64::from)),
        _ => CifValue::Unknown,
    }
}

fn component_atom_value(
    (component, atom): &(&ComponentDefinition, &ComponentAtom),
    item: &str,
) -> CifValue {
    match item.to_ascii_lowercase().as_str() {
        "comp_id" => CifValue::text(&component.id),
        "atom_id" => CifValue::text(&atom.atom_id),
        "type_symbol" => CifValue::text(&atom.element),
        "charge" => optional_integer(atom.formal_charge.map(i64::from)),
        _ => CifValue::Unknown,
    }
}

fn component_bond_value(
    (component, bond): &(&ComponentDefinition, &ComponentBond),
    item: &str,
) -> CifValue {
    match item.to_ascii_lowercase().as_str() {
        "comp_id" => CifValue::text(&component.id),
        "atom_id_1" => CifValue::text(&bond.first_atom_id),
        "atom_id_2" => CifValue::text(&bond.second_atom_id),
        "value_order" => CifValue::text(&bond.order),
        "pdbx_aromatic_flag" => bond.aromatic.map_or(CifValue::Unknown, |value| {
            CifValue::text(if value { "y" } else { "n" })
        }),
        _ => CifValue::Unknown,
    }
}

fn optional_text(value: Option<&str>) -> CifValue {
    value.map_or(CifValue::Unknown, CifValue::text)
}

fn optional_float(value: Option<f64>) -> CifValue {
    value.map_or(CifValue::Unknown, |value| {
        CifValue::Float(value, None, None)
    })
}

fn optional_integer(value: Option<i64>) -> CifValue {
    value.map_or(CifValue::Unknown, |value| CifValue::Integer(value, None))
}

fn reorder_columns(entry: &CifEntry) -> CifEntry {
    let CifEntry::Loop(cif_loop) = entry else {
        return entry.clone();
    };
    if split_tag(cif_loop.tags().first().map_or("", String::as_str)).map(|(category, _)| category)
        != Some("atom_site")
    {
        return entry.clone();
    }
    let mut columns = (0..cif_loop.column_count()).collect::<Vec<_>>();
    columns.sort_by_key(|column| {
        let item = split_tag(&cif_loop.tags()[*column])
            .map(|(_, item)| item)
            .unwrap_or("");
        ATOM_SITE_ORDER
            .iter()
            .position(|expected| item.eq_ignore_ascii_case(expected))
            .unwrap_or(ATOM_SITE_ORDER.len())
    });
    let tags = columns
        .iter()
        .map(|column| cif_loop.tags()[*column].clone())
        .collect();
    let mut values = Vec::with_capacity(cif_loop.value_count());
    for row_index in 0..cif_loop.row_count() {
        if let Some(row) = cif_loop.row(row_index) {
            values.extend(
                columns
                    .iter()
                    .filter_map(|column| row.get(*column))
                    .map(|value| value.to_owned()),
            );
        }
    }
    CifEntry::Loop(CifLoop::new(tags, values))
}

fn category_rank(entry: &CifEntry, category_order: &[&str]) -> usize {
    let tag = match entry {
        CifEntry::Item(item) => item.tag(),
        CifEntry::Loop(cif_loop) => cif_loop.tags().first().map_or("", String::as_str),
        CifEntry::Frame(_) => return category_order.len() + 1,
    };
    let category = split_tag(tag).map(|(category, _)| category).unwrap_or("");
    category_order
        .iter()
        .position(|expected| category.eq_ignore_ascii_case(expected))
        .unwrap_or(category_order.len())
}
