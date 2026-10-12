use reduce3::cifsource::{CifSource, CifTable};
use reduce3::interp::{FlatAtoms, validation::ValidationGraph};
use reduce3::model::Structure;

use crate::cif::CifBlock;
use crate::reduce::BlockSource;
mod images;

/// Preserve annotated bonds that the preparation hierarchy does not import.
pub(super) fn uncovered(
    block: &CifBlock,
    structure: &Structure,
    flat: &FlatAtoms,
    graph: &ValidationGraph,
) -> Vec<String> {
    let source = BlockSource::new(block);
    let (Some(connections), Some(atoms)) = (source.table("struct_conn"), source.table("atom_site"))
    else {
        return Vec::new();
    };
    let mut reasons = Vec::new();
    for row in 0..connections.row_count() {
        let kind = connections
            .column("conn_type_id")
            .map(|column| connections.cell(row, column));
        if kind.is_some_and(|value| {
            let kind = value.trim();
            !kind.starts_with("covale") && !matches!(kind, "disulf" | "metalc" | "modres")
        }) {
            continue;
        }
        let image = match images::relative(&source, &connections, row, structure) {
            Ok(image) => image,
            Err(reason) => {
                reasons.push(format!("Connection row {row}: {reason}"));
                continue;
            }
        };
        let find = |side| {
            let fields = endpoint_fields(&connections, row, side);
            if !fields.iter().any(|(name, _)| name.ends_with("atom_id")) {
                return Vec::new();
            }
            let Some(fields) = fields
                .into_iter()
                .map(|(name, value)| atoms.column(name).map(|column| (column, value)))
                .collect::<Option<Vec<_>>>()
            else {
                return Vec::new();
            };
            flat.path
                .iter()
                .enumerate()
                .filter_map(|(index, &p)| {
                    let source_row = structure.atom(p).src as usize;
                    fields
                        .iter()
                        .all(|(column, value)| atoms.cell(source_row, *column).trim() == value)
                        .then_some(index as u32)
                })
                .collect::<Vec<_>>()
        };
        let (first, second) = (find(1), find(2));
        if first.len() != 1 || second.len() != 1 {
            reasons.push(format!(
                "Connection row {row} has missing or ambiguous atom identities"
            ));
        } else if !(graph.symmetry_bonds.iter().any(|bond| {
            bond.indices == [first[0], second[0]] && images::same(image, bond.image)
                || bond.indices == [second[0], first[0]]
                    && image
                        .inverse()
                        .is_some_and(|reverse| images::same(reverse, bond.image))
        }) || image.is_identity()
            && graph
                .bonds
                .iter()
                .any(|bond| bond.indices == [first[0].min(second[0]), first[0].max(second[0])]))
        {
            reasons.push(format!(
                "Connection row {row} lacks a weighted chemical link definition"
            ));
        }
    }
    reasons
}

fn endpoint_fields(table: &impl CifTable, row: usize, side: usize) -> Vec<(&'static str, String)> {
    let mut result = Vec::new();
    for name in [
        "label_asym_id",
        "auth_asym_id",
        "label_seq_id",
        "auth_seq_id",
        "label_comp_id",
        "auth_comp_id",
        "label_atom_id",
        "auth_atom_id",
    ] {
        let Some(column) = table.column(&format!("ptnr{side}_{name}")) else {
            continue;
        };
        if table.missing(row, column).is_none() {
            result.push((name, table.cell(row, column).trim().to_owned()));
        }
    }
    for (partner, atom) in [
        (format!("pdbx_ptnr{side}_label_alt_id"), "label_alt_id"),
        (format!("pdbx_ptnr{side}_PDB_ins_code"), "pdbx_PDB_ins_code"),
    ] {
        let Some(column) = table.column(&partner) else {
            continue;
        };
        if table.missing(row, column).is_none() {
            result.push((atom, table.cell(row, column).trim().to_owned()));
        }
    }
    result
}
