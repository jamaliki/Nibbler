//! Preserve sequence gaps and add gaps at observed polymer breaks.

use super::{record::Tokens, residue, source::Source};
use std::collections::BTreeMap;

pub(super) const INDEX_SEPARATION: i64 = 65;

/// The longest distance, in Å, between the atoms that link consecutive residues (C to N
/// for a peptide, O3' to P for a nucleotide) that still counts as a covalent link.
/// Peptide C-N distances peak at 1.3 to 1.4 Å and clear breaks start near 2.7 Å; the same
/// 2.4 Å tolerance decides which sequence links the chemistry graph keeps.
const LINK_LIMIT: f32 = 2.4;

/// The CA-CA distance, in Å, that a peptide link cannot exceed when C or N is unobserved.
const CA_LIMIT: f32 = 4.;

#[derive(Clone, Copy, Default)]
struct Backbone {
    n: Option<[f32; 3]>,
    ca: Option<[f32; 3]>,
    c: Option<[f32; 3]>,
}

impl Backbone {
    /// Whether `right`, the next residue in sequence, is cut off from `self`.
    ///
    /// The peptide bond decides: C of `self` against N of `right`, since a stretched CA-CA
    /// distance does not break a bonded pair. Only when either atom is unobserved does
    /// CA-CA stand in, so a residue without a backbone N or C keeps the behaviour it had
    /// before the peptide bond was read. Without CA as well, nothing shows a break.
    fn cut_from(self, right: Self) -> bool {
        let (a, b, limit) = match (self.c, right.n) {
            (Some(c), Some(n)) => (c, n, LINK_LIMIT),
            _ => match (self.ca, right.ca) {
                (Some(a), Some(b)) => (a, b, CA_LIMIT),
                _ => return false,
            },
        };
        residue::distance(a, b) > limit
    }
}

fn site(source: &Source, group: usize, name: &str) -> Option<[f32; 3]> {
    source.groups[group].iter().find_map(|&row| {
        let atom = &source.atoms[row];
        ((atom.name == name || name == "O3'" && atom.name == "O3*")
            && residue::finite(atom.position))
        .then_some(atom.position)
    })
}

/// Whether the residue `right`, which follows `left` in sequence, is cut off from it.
pub(super) fn distant(source: &Source, left: usize, right: usize, protein: bool) -> bool {
    if protein {
        let backbone = |group| Backbone {
            n: site(source, group, "N"),
            ca: site(source, group, "CA"),
            c: site(source, group, "C"),
        };
        return backbone(left).cut_from(backbone(right));
    }
    match (site(source, left, "O3'"), site(source, right, "P")) {
        (Some(a), Some(b)) => residue::distance(a, b) > LINK_LIMIT,
        // Missing nucleotide link atoms cannot establish backbone continuity.
        _ => true,
    }
}

pub(super) fn apply(source: &Source, mapping: &[usize], tokens: &mut Tokens) {
    let mut chains = BTreeMap::<_, BTreeMap<_, usize>>::new();
    for (group, rows) in source.groups.iter().enumerate() {
        let atom = &source.atoms[rows[0]];
        let kind = source
            .entity_types
            .get(atom.entity.as_str())
            .map(|s| s.to_ascii_lowercase());
        let nucleic = kind.is_some_and(|s| s.contains("ribonucleotide"));
        if atom.sequence > 0 && !atom.protein && nucleic {
            chains
                .entry((&atom.chain, &atom.entity))
                .or_default()
                .insert(atom.sequence, group);
        }
    }
    for groups in chains.values() {
        let mut previous: Option<(i64, usize)> = None;
        let mut segment = 0;
        for (&number, &group) in groups {
            if let Some((before, left)) = previous
                && number == before + 1
                && distant(source, left, group, false)
            {
                segment += 1;
            }
            for &row in &source.groups[group] {
                let token = mapping[row];
                if token == usize::MAX {
                    continue;
                }
                tokens.polymer_segments[token] = segment;
            }
            previous = Some((number, group));
        }
    }
}

fn observed(tokens: &Tokens, row: usize, slot: usize) -> Option<[f32; 3]> {
    let position = tokens.positions[row * 96 + slot * 3..][..3]
        .try_into()
        .unwrap_or([f32::NAN; 3]);
    (tokens.masks[row * 32 + slot] != 0 && residue::finite(position)).then_some(position)
}

pub(super) fn retained_proteins(tokens: &mut Tokens) {
    let mut chains = BTreeMap::<_, Vec<usize>>::new();
    for row in 0..tokens.identity.len() {
        if tokens.ligand[row] == 0 && tokens.masks[row * 32 + 1] != 0 {
            chains.entry(tokens.chains[row]).or_default().push(row);
        }
    }
    for rows in chains.values_mut() {
        rows.sort_unstable_by_key(|&row| tokens.sequence[row]);
        let mut previous: Option<(i64, Backbone)> = None;
        let mut offset = 0;
        for &row in rows.iter() {
            let number = tokens.sequence[row];
            let backbone = Backbone {
                n: observed(tokens, row, 0),
                ca: observed(tokens, row, 1),
                c: observed(tokens, row, 2),
            };
            if let Some((before, left)) = previous
                && number == before + 1
                && left.cut_from(backbone)
            {
                offset += (INDEX_SEPARATION - (number - before)).max(0);
            }
            tokens.sequence[row] = number + offset;
            previous = Some((number, backbone));
        }
    }
}

#[cfg(test)]
#[allow(clippy::expect_used)]
mod tests;
