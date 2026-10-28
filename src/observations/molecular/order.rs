//! Restore the indexed token order after protonation reorders source atoms.

use super::Failure;
use super::json;
use super::pack::Packed;
use crate::hash::{FastMap, FastSet};
use std::collections::{HashMap, HashSet};

fn listing(keys: &[&[&str]]) -> String {
    let items: Vec<String> = keys
        .iter()
        .take(3)
        .map(|key| json::repr_sequence(key.iter().copied(), true))
        .collect();
    format!("[{}]", items.join(", "))
}

/// Token keys as compact JSON, written as `json.dumps` writes them.
pub(super) fn key_text(keys: &[[&str; 5]]) -> String {
    let mut out = String::with_capacity(keys.len() * 40);
    out.push('[');
    for (index, key) in keys.iter().enumerate() {
        if index > 0 {
            out.push(',');
        }
        json::strings(&mut out, key.iter().copied());
    }
    out.push(']');
    out
}

/// Reorder token rows and their edges into the indexed order, given as compact JSON.
///
/// Returns the JSON of the final keys. When protonation keeps the order, the keys
/// written as JSON equal the indexed text and nothing needs parsing.
pub(super) fn indexed(packed: &mut Packed<'_>, expected: &str) -> Result<String, Failure> {
    let text = key_text(&packed.keys);
    if text == expected {
        let mut seen = FastSet::with_capacity_and_hasher(packed.keys.len(), Default::default());
        if !packed.keys.iter().all(|key| seen.insert(key)) {
            return Err("Molecular token identities must be unique before alignment".into());
        }
        return Ok(text);
    }
    if let Some(selected) = permutation(&packed.keys, expected) {
        if selected
            .iter()
            .enumerate()
            .all(|(index, &old)| index == old)
        {
            return Ok(text);
        }
        reorder(packed, &selected);
        return Ok(key_text(&packed.keys));
    }
    let expected: Vec<Vec<String>> =
        serde_json::from_str(expected).map_err(|error| error.to_string())?;
    let wanted_keys: Vec<Vec<&str>> = expected
        .iter()
        .map(|key| key.iter().map(String::as_str).collect())
        .collect();
    let actual: HashSet<&[&str]> = packed.keys.iter().map(|key| key.as_slice()).collect();
    let wanted: HashSet<&[&str]> = wanted_keys.iter().map(Vec::as_slice).collect();
    if actual.len() != packed.keys.len() || wanted.len() != expected.len() {
        return Err("Molecular token identities must be unique before alignment".into());
    }
    if actual != wanted {
        let mut missing: Vec<&[&str]> = wanted.difference(&actual).copied().collect();
        let mut added: Vec<&[&str]> = actual.difference(&wanted).copied().collect();
        missing.sort();
        added.sort();
        return Err(Failure::Value(format!(
            "Protonation changed molecular token identities: missing={}, added={}",
            listing(&missing),
            listing(&added)
        )));
    }
    if packed
        .keys
        .iter()
        .map(|key| key.as_slice())
        .eq(wanted_keys.iter().map(Vec::as_slice))
    {
        return Ok(text);
    }
    let lookup: HashMap<&[&str], usize> = packed
        .keys
        .iter()
        .enumerate()
        .map(|(index, key)| (key.as_slice(), index))
        .collect();
    let selected: Vec<usize> = wanted_keys
        .iter()
        .map(|key| lookup.get(key.as_slice()).copied().unwrap_or(0))
        .collect();
    reorder(packed, &selected);
    Ok(key_text(&packed.keys))
}

/// For each indexed key, the row that holds it, when the keys are unique and both lists hold
/// the same keys.
///
/// Returns `None` for any other case, and for JSON with escapes or keys that are not five
/// strings, which the general path then handles or reports.
fn permutation(keys: &[[&str; 5]], expected: &str) -> Option<Vec<usize>> {
    let wanted: Vec<[&str; 5]> = serde_json::from_str(expected).ok()?;
    if wanted.len() != keys.len() {
        return None;
    }
    let mut lookup = FastMap::with_capacity_and_hasher(keys.len(), Default::default());
    for (index, key) in keys.iter().enumerate() {
        if lookup.insert(key, index).is_some() {
            return None;
        }
    }
    let mut used = vec![false; keys.len()];
    let mut selected = Vec::with_capacity(keys.len());
    for key in &wanted {
        let index = *lookup.get(key)?;
        if std::mem::replace(&mut used[index], true) {
            return None;
        }
        selected.push(index);
    }
    Some(selected)
}

/// The decimal digits of `value`, left-aligned and padded with zero bytes.
///
/// These arrays sort as the decimal text sorts.
fn decimal(value: usize) -> [u8; 20] {
    let mut digits = [0_u8; 20];
    let mut length = 0;
    let mut rest = value;
    loop {
        digits[length] = b'0' + (rest % 10) as u8;
        length += 1;
        rest /= 10;
        if rest == 0 {
            break;
        }
    }
    digits[..length].reverse();
    digits
}

/// Put `values[selected[index]]` at `index` for every index, where `selected` is a permutation.
///
/// Each cycle of the permutation moves by swaps in place, so rows that keep their place cost
/// nothing and no second copy of the rows is allocated.
fn permute<T>(values: &mut [T], selected: &[usize], done: &mut [bool]) {
    done.fill(false);
    for start in 0..selected.len() {
        let mut current = start;
        while !done[current] {
            done[current] = true;
            let next = selected[current];
            if next == start {
                break;
            }
            values.swap(current, next);
            current = next;
        }
    }
}

/// Move row `selected[index]` to `index` and renumber every token reference to match.
fn reorder(packed: &mut Packed<'_>, selected: &[usize]) {
    let mut inverse = vec![0_i64; selected.len()];
    for (index, &old) in selected.iter().enumerate() {
        inverse[old] = index as i64;
    }
    let mut done = vec![false; selected.len()];
    permute(&mut packed.rows, selected, &mut done);
    if !packed.atom_sasa.is_empty() {
        permute(&mut packed.atom_sasa, selected, &mut done);
        permute(&mut packed.atom_sasa_mask, selected, &mut done);
    }
    let remap = |value: i64| {
        if value >= 0 {
            inverse[value as usize]
        } else {
            -1
        }
    };
    for value in packed
        .legacy_token_index
        .iter_mut()
        .chain(packed.owner.iter_mut())
    {
        *value = remap(*value);
    }
    for row in &mut packed.rows {
        row.ncaa_parent_token = remap(row.ncaa_parent_token);
    }
    for pair in packed
        .protein_attachment_index
        .iter_mut()
        .chain(packed.bridge_index.iter_mut())
        .chain(packed.ligand_bonds.iter_mut())
    {
        *pair = pair.map(|token| inverse[token as usize]);
    }
    for pair in &mut packed.attachment_pairs {
        pair[0] = inverse[pair[0] as usize];
        pair[1] = inverse[pair[1] as usize];
    }
    for (token, _) in &mut packed.names {
        *token = inverse[*token] as usize;
    }
    packed
        .names
        .sort_by_cached_key(|(token, _)| decimal(*token));
    // The keys now hold the indexed values in the indexed order.
    permute(&mut packed.keys, selected, &mut done);
}
