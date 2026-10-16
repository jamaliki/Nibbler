//! Read selected atom rows and the covalent links used by training targets.

use std::collections::{BTreeSet, HashMap, HashSet};

use crate::cif::{BlockCategory, CifDocument, CifValueRef, numeric};

pub(super) struct Atom {
    pub source_row: u32,
    pub name: String,
    pub element: String,
    pub chain: String,
    pub component: String,
    pub entity: String,
    pub auth_sequence: String,
    pub auth_chain: String,
    pub auth_component: String,
    pub insertion: String,
    pub alternate: String,
    pub sequence: i64,
    pub position: [f32; 3],
    pub protein: bool,
    pub group: usize,
}

pub(super) struct Source {
    pub atoms: Vec<Atom>,
    pub groups: Vec<Vec<usize>>,
    pub bonds: Vec<(usize, usize)>,
    pub hydrogen_bonds: Vec<(usize, Option<usize>)>,
    pub has_protein_annotations: bool,
    pub lengths: HashMap<String, i64>,
    pub entity_types: HashMap<String, String>,
}

type Occurrence = (String, i64, String, String, String, String);
type Residue = (String, i64, String, String);
type Endpoint = (String, i64, String, String);
type Connections = (Vec<(usize, usize)>, Vec<(usize, Option<usize>)>);

fn text(table: &BlockCategory<'_>, row: usize, column: Option<usize>) -> String {
    match column.and_then(|c| table.value(row, c)) {
        Some(CifValueRef::Unknown) => "?".to_owned(),
        Some(CifValueRef::NotApplicable) | None => ".".to_owned(),
        _ => column
            .and_then(|c| table.text(row, c))
            .map_or_else(|| ".".to_owned(), |s| s.into_owned()),
    }
}

fn normal(value: &str) -> &str {
    let value = value.trim();
    if matches!(value, "" | "." | "?") {
        ""
    } else {
        value
    }
}

fn number(table: &BlockCategory<'_>, row: usize, column: Option<usize>) -> Option<f64> {
    match column.and_then(|c| table.value(row, c))? {
        CifValueRef::Integer(value, _) => Some(value as f64),
        CifValueRef::Float(value, _, _) => value.is_finite().then_some(value),
        CifValueRef::Text(value) => numeric::parse_float(value.as_str()),
        CifValueRef::Unknown | CifValueRef::NotApplicable => None,
    }
}

fn integer(table: &BlockCategory<'_>, row: usize, column: Option<usize>) -> Option<i64> {
    match column.and_then(|c| table.value(row, c))? {
        CifValueRef::Integer(value, _) => Some(value),
        CifValueRef::Text(value) => numeric::parse_integer(value.as_str()),
        _ => None,
    }
}

fn checked_number(
    table: &BlockCategory<'_>,
    row: usize,
    column: Option<usize>,
    name: &str,
) -> Result<Option<f64>, String> {
    let number = number(table, row, column);
    let present = column
        .and_then(|c| table.value(row, c))
        .is_some_and(|value| !matches!(value, CifValueRef::Unknown | CifValueRef::NotApplicable));
    if present && number.is_none() {
        return Err(format!(
            "Invalid {name} number in atom_site row {}",
            row + 1
        ));
    }
    Ok(number)
}

fn entity_types(entries: &[crate::cif::CifEntry]) -> HashMap<String, String> {
    let mut result = HashMap::new();
    if let Some(table) = BlockCategory::find(entries, "entity_poly") {
        let entity = table.column("entity_id");
        let kind = table.column("type");
        for row in 0..table.row_count() {
            let entity = text(&table, row, entity);
            let kind = text(&table, row, kind);
            if !normal(&entity).is_empty() && !normal(&kind).is_empty() {
                result.insert(entity, kind);
            }
        }
    }
    result
}

fn entity_lengths(entries: &[crate::cif::CifEntry]) -> HashMap<String, i64> {
    let mut result = HashMap::new();
    if let Some(table) = BlockCategory::find(entries, "entity_poly_seq") {
        let entity = table.column("entity_id");
        let num = table.column("num");
        for row in 0..table.row_count() {
            if let Some(number) = integer(&table, row, num) {
                let entity = text(&table, row, entity);
                result
                    .entry(entity)
                    .and_modify(|n: &mut i64| *n = (*n).max(number))
                    .or_insert(number.max(0));
            }
        }
    }
    result
}

fn endpoint(
    table: &BlockCategory<'_>,
    row: usize,
    partner: usize,
    atoms: &[Atom],
    lookup: &HashMap<Endpoint, Vec<usize>>,
) -> Result<Option<usize>, String> {
    let field = |name: &str| text(table, row, table.column(&format!("ptnr{partner}_{name}")));
    let chain = field("label_asym_id");
    let component = field("label_comp_id");
    let name = field("label_atom_id");
    if [chain.as_str(), component.as_str(), name.as_str()]
        .iter()
        .any(|s| normal(s).is_empty())
    {
        return Ok(None);
    }
    let sequence = numeric::parse_integer(normal(&field("label_seq_id"))).unwrap_or(-1);
    let alternate = text(
        table,
        row,
        table.column(&format!("pdbx_ptnr{partner}_label_alt_id")),
    );
    let author = field("auth_seq_id");
    let author_chain = field("auth_asym_id");
    let author_component = field("auth_comp_id");
    let insertion = text(table, row, table.column(&format!("pdbx_ptnr{partner}_PDB_ins_code")));
    let candidates = lookup.get(&(
            normal(&chain).to_owned(),
            sequence,
            normal(&component).to_owned(),
            normal(&name).to_owned(),
        ));
    let selected: Vec<_> = candidates.into_iter().flatten()
        .copied()
        .filter(|&i| {
            (normal(&alternate).is_empty() || normal(&atoms[i].alternate) == normal(&alternate))
                && (normal(&author).is_empty() || normal(&atoms[i].auth_sequence) == normal(&author))
                && (normal(&author_chain).is_empty() || normal(&atoms[i].auth_chain) == normal(&author_chain))
                && (normal(&author_component).is_empty() || normal(&atoms[i].auth_component) == normal(&author_component))
                && (normal(&insertion).is_empty() || normal(&atoms[i].insertion) == normal(&insertion))
        }).collect();
    if selected.len() > 1 {
        return Err(format!("Connection row {} partner {partner} has ambiguous atom identity", row + 1));
    }
    Ok(selected.first().copied())
}

fn connections(entries: &[crate::cif::CifEntry], atoms: &[Atom]) -> Result<Connections, String> {
    let mut pairs = BTreeSet::new();
    let mut hydrogen_pairs = BTreeSet::new();
    let Some(table) = BlockCategory::find(entries, "struct_conn") else {
        return Ok((Vec::new(), Vec::new()));
    };
    let mut lookup: HashMap<Endpoint, Vec<usize>> = HashMap::new();
    for (index, atom) in atoms.iter().enumerate() {
        lookup
            .entry((
                normal(&atom.chain).to_owned(),
                atom.sequence.max(-1),
                normal(&atom.component).to_owned(),
                normal(&atom.name).to_owned(),
            ))
            .or_default()
            .push(index);
    }
    let kind = table.column("conn_type_id");
    let symmetries = [
        table.column("ptnr1_symmetry"),
        table.column("ptnr2_symmetry"),
    ];
    for row in 0..table.row_count() {
        let kind = normal(&text(&table, row, kind)).to_ascii_lowercase();
        if !(kind.starts_with("covale") || kind.starts_with("disulf")) {
            continue;
        }
        if symmetries.iter().any(|&c| {
            let value = text(&table, row, c);
            !matches!(normal(&value), "" | "1_555")
        }) {
            continue;
        }
        let first = endpoint(&table, row, 1, atoms, &lookup)?;
        let second = endpoint(&table, row, 2, atoms, &lookup)?;
        let hydrogen = |row: usize| matches!(atoms[row].element.as_str(), "H" | "D" | "T");
        if let (Some(a), Some(b)) = (first, second) {
            if hydrogen(a) != hydrogen(b) {
                hydrogen_pairs.insert(if hydrogen(a) {
                    (a, Some(b))
                } else {
                    (b, Some(a))
                });
            }
            if atoms[a].group != atoms[b].group {
                pairs.insert((a.min(b), a.max(b)));
            }
        } else if let Some(child) = first.or(second).filter(|&row| hydrogen(row)) {
            hydrogen_pairs.insert((child, None));
        }
    }
    Ok((
        pairs.into_iter().collect(),
        hydrogen_pairs.into_iter().collect(),
    ))
}

pub(super) fn read(document: &CifDocument, model: i64) -> Result<Source, String> {
    let block = document
        .blocks()
        .iter()
        .find(|block| BlockCategory::find(block.entries(), "atom_site").is_some())
        .ok_or_else(|| "Source has no atom_site table".to_owned())?;
    if BlockCategory::find(block.entries(), "pdbx_struct_assembly_gen")
        .is_some_and(|table| table.column("assembly_id").is_some() && table.row_count() > 0)
    {
        return Err("Expand PDB assemblies before tokenization and protonation".to_owned());
    }
    let table = BlockCategory::find(block.entries(), "atom_site")
        .ok_or_else(|| "Source has no atom_site table".to_owned())?;
    let types = entity_types(block.entries());
    let columns = [
        "label_atom_id",
        "type_symbol",
        "label_asym_id",
        "label_comp_id",
        "label_entity_id",
        "auth_seq_id",
        "pdbx_PDB_ins_code",
        "label_alt_id",
        "group_PDB",
        "auth_asym_id",
        "auth_comp_id",
    ]
    .map(|name| table.column(name));
    let model_column = table.column("pdbx_PDB_model_num");
    let sequence_column = table.column("label_seq_id");
    let coordinates = ["Cartn_x", "Cartn_y", "Cartn_z"].map(|name| table.column(name));
    let occupancy_column = table.column("occupancy");
    let mut atoms = Vec::new();
    let mut occupancy = Vec::new();
    let mut groups = Vec::new();
    let mut occurrences = HashMap::<Occurrence, usize>::new();
    let mut alternatives = HashMap::<Residue, HashSet<String>>::new();
    for row in 0..table.row_count() {
        if integer(&table, row, model_column).unwrap_or(-1) != model {
            continue;
        }
        let fields = columns.map(|column| text(&table, row, column));
        let sequence = integer(&table, row, sequence_column).unwrap_or(-1);
        let occurrence = (
            normal(&fields[2]).to_owned(),
            sequence.max(-1),
            normal(&fields[3]).to_owned(),
            normal(&fields[4]).to_owned(),
            normal(&fields[6]).to_owned(),
            normal(&fields[5]).to_owned(),
        );
        let residue = (
            occurrence.0.clone(),
            occurrence.1,
            occurrence.4.clone(),
            occurrence.5.clone(),
        );
        if !normal(&fields[7]).is_empty() {
            let choices = alternatives.entry(residue).or_default();
            choices.insert(normal(&fields[7]).to_owned());
            if choices.len() > 1 {
                return Err("Resolve alternate conformers before native target packing".to_owned());
            }
        }
        let group = *occurrences.entry(occurrence).or_insert_with(|| {
            groups.push(Vec::new());
            groups.len() - 1
        });
        let protein = if types.is_empty() {
            normal(&fields[8]).eq_ignore_ascii_case("ATOM")
        } else {
            types
                .get(normal(&fields[4]))
                .is_some_and(|t| t.to_ascii_lowercase().starts_with("polypeptide"))
        };
        let mut position = [f32::NAN; 3];
        for (axis, column) in coordinates.iter().copied().enumerate() {
            position[axis] =
                checked_number(&table, row, column, "coordinate")?.map_or(f32::NAN, |n| n as f32);
        }
        groups[group].push(atoms.len());
        atoms.push(Atom {
            source_row: u32::try_from(row)
                .map_err(|_| "atom_site row exceeds the source index limit")?,
            name: fields[0].clone(),
            element: fields[1].clone(),
            chain: fields[2].clone(),
            component: fields[3].clone(),
            entity: fields[4].clone(),
            auth_sequence: fields[5].clone(),
            auth_chain: fields[9].clone(),
            auth_component: fields[10].clone(),
            insertion: fields[6].clone(),
            alternate: fields[7].clone(),
            sequence,
            position,
            protein,
            group,
        });
        occupancy.push(
            checked_number(&table, row, occupancy_column, "occupancy")?.map_or(0.0, |n| n as f32),
        );
    }
    if atoms.is_empty() {
        return Err(format!("model_id {model} is not present in the source"));
    }
    // A selected document can retain shared and selected alternate rows for one name.
    let mut selected = vec![false; atoms.len()];
    for rows in &groups {
        let mut by_name = HashMap::<String, usize>::new();
        for &row in rows {
            let name = if normal(&atoms[row].name).is_empty() {
                format!("__row_{row}")
            } else {
                normal(&atoms[row].name).to_owned()
            };
            let previous = by_name.get(&name).copied();
            let score = |i: usize| {
                (
                    !normal(&atoms[i].alternate).is_empty(),
                    if occupancy[i].is_finite() {
                        -occupancy[i]
                    } else {
                        0.0
                    },
                    i,
                )
            };
            if previous.is_none_or(|i| score(row) < score(i)) {
                by_name.insert(name, row);
            }
        }
        for row in by_name.into_values() {
            selected[row] = true;
        }
    }
    let atoms: Vec<_> = atoms
        .into_iter()
        .enumerate()
        .filter_map(|(row, atom)| selected[row].then_some(atom))
        .collect();
    for rows in &mut groups {
        rows.clear();
    }
    for (row, atom) in atoms.iter().enumerate() {
        groups[atom.group].push(row);
    }
    let (bonds, hydrogen_bonds) = connections(block.entries(), &atoms)?;
    Ok(Source {
        atoms,
        groups,
        bonds,
        hydrogen_bonds,
        has_protein_annotations: !types.is_empty(),
        lengths: entity_lengths(block.entries()),
        entity_types: types,
    })
}

#[cfg(test)]
#[path = "source_tests.rs"]
mod tests;
