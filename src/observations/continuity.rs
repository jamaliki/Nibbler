//! Preserve sequence gaps and add gaps at observed polymer breaks.

use super::{record::Tokens, residue, source::Source};
use std::collections::BTreeMap;

pub(super) const INDEX_SEPARATION: i64 = 65;

fn site(source: &Source, group: usize, name: &str) -> Option<[f32; 3]> {
    source.groups[group].iter().find_map(|&row| {
        let atom = &source.atoms[row];
        ((atom.name == name || name == "O3'" && atom.name == "O3*")
            && residue::finite(atom.position))
            .then_some(atom.position)
    })
}

pub(super) fn distant(source: &Source, left: usize, right: usize, protein: bool) -> bool {
    let (a, b, limit) = if protein { ("CA", "CA", 4.) } else { ("O3'", "P", 2.4) };
    match (site(source, left, a), site(source, right, b)) {
        (Some(a), Some(b)) => residue::distance(a, b) > limit,
        // Missing nucleotide link atoms cannot establish backbone continuity.
        _ => !protein,
    }
}

pub(super) fn apply(source: &Source, mapping: &[usize], tokens: &mut Tokens) {
    let mut chains = BTreeMap::<_, BTreeMap<_, usize>>::new();
    for (group, rows) in source.groups.iter().enumerate() {
        let atom = &source.atoms[rows[0]];
        let kind = source.entity_types.get(&atom.entity).map(|s| s.to_ascii_lowercase());
        let nucleic = kind.is_some_and(|s| s.contains("ribonucleotide"));
        if atom.sequence > 0 && (atom.protein || nucleic) {
            chains.entry((&atom.chain, &atom.entity)).or_default().insert(atom.sequence, group);
        }
    }
    for groups in chains.values() {
        let mut previous: Option<(i64, usize)> = None;
        let (mut offset, mut segment) = (0, 0);
        for (&number, &group) in groups {
            let atom = &source.atoms[source.groups[group][0]];
            if let Some((before, left)) = previous {
                let gap = number > before + 1;
                let broken = distant(source, left, group, atom.protein);
                if gap || broken {
                    segment += 1;
                    if broken { offset += (INDEX_SEPARATION - (number - before)).max(0); }
                }
            }
            for &row in &source.groups[group] {
                let token = mapping[row];
                if token == usize::MAX { continue; }
                if atom.protein {
                    tokens.sequence[token] = number + offset;
                } else {
                    tokens.polymer_segments[token] = segment;
                }
            }
            previous = Some((number, group));
        }
    }
}

#[cfg(test)]
mod tests;
