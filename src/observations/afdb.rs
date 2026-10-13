//! Construct a protein document from dense heavy-atom arrays.

use super::{Schema, residue::distance};
use crate::cif::{CifBlock, CifDocument, CifEntry, CifItem, CifLoop, CifValue};

fn table(category: &str, names: &[&str], values: Vec<CifValue>) -> CifEntry {
    CifEntry::Loop(CifLoop::new(
        names
            .iter()
            .map(|name| format!("_{category}.{name}"))
            .collect(),
        values,
    ))
}

pub(crate) fn document(
    schema: &Schema,
    positions: &[u8],
    masks: &[u8],
    identities: &[u8],
    identifier: &str,
) -> Result<(CifDocument, usize, String), String> {
    if identifier.is_empty()
        || identifier
            .chars()
            .any(|value| value.is_whitespace() || value == '\0')
    {
        return Err(
            "Protein document identifier must be nonempty and contain no whitespace".into(),
        );
    }
    if !identities.len().is_multiple_of(8) {
        return Err("Residue identities require int64 buffers".into());
    }
    let count = identities.len() / 8;
    if count == 0
        || positions.len() != count * 14 * 3 * 4
        || masks.len() != count * 14
        || masks.iter().any(|&value| value > 1)
    {
        return Err("Expected nonempty float32 Atom14 positions and boolean masks".into());
    }
    let aa: Vec<_> = identities
        .chunks_exact(8)
        .map(|value| i64::from_ne_bytes(value.try_into().unwrap_or([0; 8])))
        .collect();
    if aa.iter().any(|&value| !(0..20).contains(&value)) {
        return Err("AFDB requires canonical residue identities".into());
    }
    let xyz: Vec<_> = positions
        .chunks_exact(4)
        .map(|value| f32::from_ne_bytes(value.try_into().unwrap_or([0; 4])))
        .collect();
    let mut atom_values = Vec::new();
    let mut chains = Vec::<String>::new();
    let mut sequences = Vec::<String>::new();
    let mut chain = 0;
    let mut offset = 0;
    let mut atom_id = 0;
    let letters = b"ARNDCQEGHILKMFPSTWYV";
    let mut sequence = String::with_capacity(count);
    for row in 0..count {
        let continuous = if row == 0 {
            false
        } else {
            let c = xyz[(row - 1) * 42 + 6..][..3]
                .try_into()
                .unwrap_or([f32::NAN; 3]);
            let n = xyz[row * 42..][..3].try_into().unwrap_or([f32::NAN; 3]);
            let value = distance(c, n);
            masks[(row - 1) * 14 + 2] != 0 && masks[row * 14] != 0 && value > 0.8 && value < 1.8
        };
        if !continuous {
            chain += 1;
            offset = row;
            chains.push(format!("A{chain}"));
            sequences.push(String::new());
        }
        let identity = aa[row] as usize;
        let residue = &schema.residues[identity];
        let letter = char::from(letters[identity]);
        sequence.push(letter);
        sequences[chain - 1].push(letter);
        for slot in 0..14 {
            if masks[row * 14 + slot] == 0 {
                continue;
            }
            let name = &residue.names[slot];
            if name.is_empty() {
                return Err("AFDB observes an unnamed heavy-atom slot".into());
            }
            atom_id += 1;
            let seq = (row - offset + 1) as i64;
            let text = |value: &str| CifValue::text(value);
            atom_values.extend([
                text("ATOM"),
                CifValue::Integer(atom_id, None),
                text(
                    &name
                        .chars()
                        .next()
                        .ok_or("Empty AFDB atom name")?
                        .to_string(),
                ),
                text(name),
                CifValue::NotApplicable,
                text(&residue.residue),
                text(&chains[chain - 1]),
                CifValue::Integer(chain as i64, None),
                CifValue::Integer(seq, None),
                text(name),
                text(&residue.residue),
                text(&chains[chain - 1]),
                CifValue::Integer(seq, None),
                CifValue::NotApplicable,
            ]);
            for axis in 0..3 {
                let value = xyz[row * 42 + slot * 3 + axis];
                if !value.is_finite() {
                    return Err("AFDB observed coordinates must be finite".into());
                }
                atom_values.push(CifValue::Float(f64::from(value), None, None));
            }
            atom_values.extend([text("1"), text("0"), text("1")]);
        }
    }
    let mut entries = vec![CifEntry::Item(CifItem::new(
        "_entry.id".into(),
        CifValue::text(identifier),
    ))];
    let mut entity = Vec::new();
    let mut poly = Vec::new();
    let mut asym = Vec::new();
    for (index, name) in chains.iter().enumerate() {
        let id = (index + 1).to_string();
        entity.extend([CifValue::text(&id), CifValue::text("polymer")]);
        poly.extend([
            CifValue::text(&id),
            CifValue::text("polypeptide(L)"),
            CifValue::text(&sequences[index]),
        ]);
        asym.extend([CifValue::text(name), CifValue::text(&id)]);
    }
    entries.push(table("entity", &["id", "type"], entity));
    entries.push(table(
        "entity_poly",
        &["entity_id", "type", "pdbx_seq_one_letter_code_can"],
        poly,
    ));
    entries.push(table("struct_asym", &["id", "entity_id"], asym));
    entries.push(table(
        "atom_site",
        &[
            "group_PDB",
            "id",
            "type_symbol",
            "label_atom_id",
            "label_alt_id",
            "label_comp_id",
            "label_asym_id",
            "label_entity_id",
            "label_seq_id",
            "auth_atom_id",
            "auth_comp_id",
            "auth_asym_id",
            "auth_seq_id",
            "pdbx_PDB_ins_code",
            "Cartn_x",
            "Cartn_y",
            "Cartn_z",
            "occupancy",
            "B_iso_or_equiv",
            "pdbx_PDB_model_num",
        ],
        atom_values,
    ));
    Ok((
        CifDocument::new(vec![CifBlock::data(identifier.into(), entries)]),
        chain,
        sequence,
    ))
}
