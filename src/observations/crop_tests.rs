#![allow(clippy::unwrap_used)]

use super::{
    Cells, CropPlans, SWEEP_MIN_TOKENS, Tally, View, compare, finish, permutation, plan,
    plan_with_chain_policy,
};

fn fixture(count: usize, ligand: &[u8], chains: &[i64]) -> super::CropPlans {
    let positions: Vec<u8> = (0..count)
        .flat_map(|index| [index as f32, 0., 0.])
        .flat_map(|value| value.to_ne_bytes())
        .collect();
    let chains: Vec<u8> = chains
        .iter()
        .flat_map(|value| value.to_ne_bytes())
        .collect();
    plan(&positions, &chains, ligand, &vec![1; count], 50, None).unwrap()
}

#[test]
fn stable_ties_and_short_chain_filters() {
    let result = fixture(60, &vec![0; 60], &vec![0; 60]);
    assert_eq!(result.centers.len(), 60);
    assert_eq!(result.ties[25], 0);
    assert_eq!(result.radii[25], 625.);
    assert_eq!(result.bits, vec![1; 60]);
    assert!((result.mass - 1.).abs() < 1e-14);
    let result = fixture(
        60,
        &vec![0; 60],
        &(0..60).map(|index| (index / 30) as i64).collect::<Vec<_>>(),
    );
    assert!(result.centers.is_empty());
}

#[test]
fn empty_and_protein_free_views_have_no_mass() {
    assert_eq!(plan(&[], &[], &[], &[], 50, None).unwrap().mass, 0.);
    assert_eq!(fixture(60, &vec![1; 60], &vec![0; 60]).mass, 0.);
}

#[test]
fn rejects_invalid_buffers() {
    assert!(plan(&[], &[], &[0], &[1], 50, None).is_err());
    assert!(plan(&[0; 12], &(-1_i64).to_ne_bytes(), &[0], &[1], 50, None).is_err());
    assert!(plan(&[0; 12], &(0_i64).to_ne_bytes(), &[2], &[1], 50, None).is_err());
    assert!(plan(&[0; 12], &(0_i64).to_ne_bytes(), &[0], &[1], 49, None).is_err());
    assert!(plan(&[0; 12], &(0_i64).to_ne_bytes(), &[0], &[1], 50, Some(&[])).is_err());
    assert!(plan(&[0; 12], &(0_i64).to_ne_bytes(), &[0], &[1], 50, Some(&[2])).is_err());
}

#[test]
fn sparse_chain_ids_survive_complete_residue_filtering() {
    let result = fixture(60, &vec![0; 60], &vec![130; 60]);
    assert_eq!(result.centers.len(), 60);
    assert_eq!(result.words, 3);
    assert_eq!(result.bits, [0, 0, 4].repeat(60));
    assert!((result.mass - 1.).abs() < 1e-14);
    assert_eq!(fixture(17, &vec![0; 17], &vec![41; 17]).mass, 0.);
    assert_eq!(fixture(2, &[1, 1], &[1, 2]).mass, 0.);
}

#[test]
fn masked_centers_stay_in_the_neighbor_budget_and_cohort() {
    let positions: Vec<u8> = (0..53).flat_map(|_| [0_f32; 3]).flat_map(f32::to_ne_bytes).collect();
    let chains: Vec<u8> = (0..53).map(|row| if row < 50 { 0_i64 } else { 1 }).flat_map(i64::to_ne_bytes).collect();
    let ligand: Vec<_> = (0..53).map(|row| u8::from(row >= 50)).collect();
    let mut eligible = vec![1; 53];
    eligible[52] = 0;
    let result = plan(&positions, &chains, &ligand, &vec![1; 53], 53, Some(&eligible)).unwrap();
    assert_eq!(result.centers, (0..52).collect::<Vec<_>>());
    assert_eq!(result.bits, vec![3; 52]);
    assert!((result.mass - 1.).abs() < 1e-14);
    assert!((result.probabilities[50] - (0.01 / 52. + 0.99 / 2.)).abs() < 1e-14);
    let limited = plan(&positions, &chains, &ligand, &vec![1; 53], 52, Some(&eligible)).unwrap();
    assert_eq!(limited.bits, vec![1; 52]);
    assert_eq!(limited.ties, vec![51; 52]);
    let none = plan(&positions, &chains, &ligand, &vec![1; 53], 53, Some(&vec![0; 53])).unwrap();
    assert_eq!(none.mass, 0.);
    assert!(none.centers.is_empty());
}

#[test]
fn ligand_permutation_keeps_protein_and_chain_coordinates() {
    let bytes = |values: &[f32]| {
        values
            .iter()
            .flat_map(|value| value.to_ne_bytes())
            .collect::<Vec<_>>()
    };
    let left = bytes(&[0., 0., 0., 1., 0., 0., 2., 0., 0.]);
    let right = bytes(&[0., 0., 0., 2., 0., 0., 1., 0., 0.]);
    let chains: Vec<u8> = [0_i64, 1, 1]
        .into_iter()
        .flat_map(|value| value.to_ne_bytes())
        .collect();
    assert_eq!(
        permutation(&left, &right, &chains, &[0, 1, 1]).unwrap(),
        vec![0, 2, 1]
    );
    assert!(permutation(&left, &right, &chains, &[0, 0, 1]).is_err());
    let different = bytes(&[0., 0., 0., 3., 0., 0., 1., 0., 0.]);
    assert!(permutation(&left, &different, &chains, &[0, 1, 1]).is_err());
}

/// Deterministic xorshift values in [0, 1).
fn uniform(state: &mut u64) -> f32 {
    *state ^= *state << 13;
    *state ^= *state >> 7;
    *state ^= *state << 17;
    (*state >> 40) as f32 / (1_u64 << 24) as f32
}

fn layouts(seed: u64, count: usize) -> Vec<Vec<[f32; 3]>> {
    let mut state = seed | 1;
    let mut clustered = Vec::new();
    while clustered.len() < count {
        let centre = [
            uniform(&mut state) * 200.,
            uniform(&mut state) * 200.,
            uniform(&mut state) * 40.,
        ];
        let spread = 1. + uniform(&mut state) * 12.;
        for _ in 0..(1 + (uniform(&mut state) * 300.) as usize).min(count - clustered.len()) {
            clustered.push(std::array::from_fn(|axis| {
                centre[axis] + (uniform(&mut state) - 0.5) * spread
            }));
        }
    }
    // Integer lattices and repeated points produce many equal distances.
    let lattice = (0..count)
        .map(|index| {
            [
                (index % 13) as f32,
                (index / 13 % 11) as f32,
                (index / 143) as f32,
            ]
        })
        .collect();
    let repeated = (0..count)
        .map(|index| {
            let value = (index % 37) as f32 * 0.5;
            [value, -value, 3.]
        })
        .collect();
    let planar = (0..count)
        .map(|_| [uniform(&mut state) * 90., uniform(&mut state) * 90., 7.25])
        .collect();
    let mut outliers: Vec<[f32; 3]> = (0..count)
        .map(|_| std::array::from_fn(|_| uniform(&mut state) * 20.))
        .collect();
    outliers[3] = [5000., -4000., 120.];
    outliers[count / 2] = [-3e4, 2e4, 9e3];
    let tiny = (0..count)
        .map(|_| std::array::from_fn(|_| uniform(&mut state) * 1e-4))
        .collect();
    vec![clustered, lattice, repeated, planar, outliers, tiny]
}

/// Roughly normal values with mean zero and unit variance.
fn normal(state: &mut u64) -> f32 {
    ((0..12).map(|_| uniform(state)).sum::<f32>() - 6.) * 1.
}

/// Tokens and chains for the sweep tests, with an optional chain policy.
struct Layout {
    positions: Vec<[f32; 3]>,
    chains: Vec<i64>,
    ligand: Vec<u8>,
    usable: Vec<u8>,
    eligible: Option<Vec<u8>>,
    /// Minimum count and protein flag of each chain.
    policy: Option<(Vec<i64>, Vec<u8>)>,
}

/// Plans from a scan over every token for each centre, with `reference_filter`, and the
/// number of centres whose crop passes the ligand limits.
fn reference_plan(view: &View) -> (CropPlans, usize) {
    let total = view.center_eligible.iter().filter(|&&value| value).count();
    let ligand_count = (0..view.positions.len())
        .filter(|&index| view.ligand[index] && view.center_eligible[index])
        .count();
    let mut result = CropPlans {
        centers: Vec::new(),
        radii: Vec::new(),
        ties: Vec::new(),
        bits: Vec::new(),
        probabilities: Vec::new(),
        words: view.chain_sizes.len().div_ceil(64).max(1),
        mass: 0.,
    };
    let (mut nearest, mut retry, mut over_count) = (Vec::new(), 0, 0);
    for center in 0..view.positions.len() {
        if !view.center_eligible[center] {
            continue;
        }
        view.nearest_all(center, &mut nearest);
        let (bits, over) = reference_filter(view, &nearest);
        retry += usize::from(over && view.ligand[center]);
        over_count += usize::from(over);
        if let Some(bits) = bits {
            let (radius, tie) = *nearest.last().unwrap();
            result.centers.push(center as i64);
            result.radii.push(radius);
            result.ties.push(tie as i64);
            result.bits.extend(bits);
        }
    }
    (finish(result, view, retry, total, ligand_count), over_count)
}

/// Check the sweep plans of `layout` against `reference_plan`, bit for bit, and return the
/// accepted and over-limit centre counts.
fn check_sweep(layout: &Layout, size: usize) -> (usize, usize) {
    let positions: Vec<u8> = layout.positions.iter().flatten().flat_map(|value| value.to_ne_bytes()).collect();
    let chains: Vec<u8> = layout.chains.iter().flat_map(|value| value.to_ne_bytes()).collect();
    let eligible = layout.eligible.as_deref();
    let minimum: Option<Vec<u8>> = layout
        .policy
        .as_ref()
        .map(|(minimum, _)| minimum.iter().flat_map(|value| value.to_ne_bytes()).collect());
    let protein = layout.policy.as_ref().map(|(_, protein)| protein.as_slice());
    let result = plan_with_chain_policy(
        &positions, &chains, &layout.ligand, &layout.usable, size, eligible, minimum.as_deref(), protein,
    )
    .unwrap();
    let mut view = View::read(&positions, &chains, &layout.ligand, &layout.usable, size, eligible).unwrap();
    if let Some((minimum, protein)) = &layout.policy {
        view.chain_minimum = minimum.iter().map(|&value| value as usize).collect();
        view.chain_protein = protein.iter().map(|&value| value != 0).collect();
    }
    // These are the conditions under which the plan uses the sweep.
    assert!(layout.positions.len() >= SWEEP_MIN_TOKENS.max(4 * size));
    assert!(Cells::build(&view).is_some());
    let (expected, over) = reference_plan(&view);
    let floats = |values: &[f32]| values.iter().map(|value| value.to_bits()).collect::<Vec<_>>();
    let doubles = |values: &[f64]| values.iter().map(|value| value.to_bits()).collect::<Vec<_>>();
    assert_eq!(result.centers, expected.centers, "centres at size {size}");
    assert_eq!(floats(&result.radii), floats(&expected.radii), "radii at size {size}");
    assert_eq!(result.ties, expected.ties, "ties at size {size}");
    assert_eq!(result.bits, expected.bits, "chain bits at size {size}");
    assert_eq!(doubles(&result.probabilities), doubles(&expected.probabilities));
    assert_eq!(result.mass.to_bits(), expected.mass.to_bits());
    (expected.centers.len(), over)
}

/// Chain IDs in blocks of `tokens`, with each listed token alone in a chain of its own.
fn chains_with_singles(blocks: impl Iterator<Item = i64>, single: &[bool]) -> Vec<i64> {
    let blocks: Vec<i64> = blocks.collect();
    let first = blocks.iter().max().map_or(0, |&value| value + 1);
    blocks
        .iter()
        .zip(single)
        .enumerate()
        .map(|(index, (&block, &alone))| if alone { first + index as i64 } else { block })
        .collect()
}

#[test]
fn sweep_matches_the_full_scan_on_ties_and_equal_coordinates() {
    let mut state = 11;
    let count = 4100;
    let mut all = layouts(5, count);
    // Signed zeros at one site that holds more tokens than any crop, then a sphere around it,
    // so that many tokens share each distance from the origin.
    all.push(
        (0..count)
            .map(|index| match index {
                0..600 => [0.; 3],
                600..1200 => [-0., 0., -0.],
                _ => {
                    let point: [f32; 3] = std::array::from_fn(|_| normal(&mut state));
                    let norm = point.iter().map(|value| value * value).sum::<f32>().sqrt();
                    point.map(|value| 10. * value / norm)
                }
            })
            .collect(),
    );
    let chains: Vec<i64> = (0..count).map(|index| (index / 60) as i64).collect();
    let ligand: Vec<u8> = chains.iter().map(|&chain| u8::from(chain % 7 == 3)).collect();
    let mut accepted = 0;
    for (index, positions) in all.into_iter().enumerate() {
        let layout = Layout {
            positions,
            chains: chains.clone(),
            ligand: ligand.clone(),
            usable: (0..count).map(|_| u8::from(uniform(&mut state) > 0.05)).collect(),
            eligible: (index % 2 == 1)
                .then(|| (0..count).map(|_| u8::from(uniform(&mut state) < 0.5)).collect()),
            policy: None,
        };
        for size in [50, 384, 512] {
            accepted += check_sweep(&layout, size).0;
        }
    }
    assert!(accepted > 5_000, "only {accepted} accepted centres");
}

#[test]
fn sweep_skips_protein_free_blocks_exactly() {
    // Protein tokens share one site, as when their positions are zeroed, far from nucleic
    // chains and single-token ligands. Most blocks then hold no protein.
    let mut state = 13;
    let count = 4100;
    let protein_rows = count / 3;
    let positions: Vec<[f32; 3]> = (0..count)
        .map(|index| {
            if index < protein_rows {
                [0.; 3]
            } else {
                std::array::from_fn(|_| 300. + 25. * normal(&mut state))
            }
        })
        .collect();
    let mut rejected = 0;
    for share in [0.2, 0.01, 0.002] {
        let single: Vec<bool> =
            (0..count).map(|index| index >= protein_rows && uniform(&mut state) < share).collect();
        let chains = chains_with_singles((0..count).map(|index| (index / 100) as i64), &single);
        let chain_count = chains.iter().max().map_or(0, |&value| value + 1) as usize;
        let mut protein = vec![0_u8; chain_count];
        for &chain in &chains[..protein_rows] {
            protein[chain as usize] = 1;
        }
        let minimum = protein.iter().map(|&value| if value == 1 { 50 } else { 1 }).collect();
        let layout = Layout {
            positions: positions.clone(),
            chains,
            ligand: single.iter().map(|&value| u8::from(value)).collect(),
            usable: (0..count)
                .map(|index| u8::from(index < protein_rows || single[index] || uniform(&mut state) < 0.9))
                .collect(),
            eligible: None,
            policy: Some((minimum, protein)),
        };
        for size in [50, 384, 512] {
            let (accepted, over) = check_sweep(&layout, size);
            assert!(accepted >= protein_rows, "{accepted} accepted at size {size}");
            rejected += count - accepted - over;
        }
    }
    assert!(rejected > 15_000, "only {rejected} centres without protein");
}

#[test]
fn sweep_proves_the_ligand_chain_limit_exactly() {
    // Dense single-token ligands, so that most crops keep more ligand chains than the limit.
    let mut state = 17;
    let count = 4100;
    let mut over_count = 0;
    for (spread, share) in [(25., 0.15), (12., 0.05), (40., 0.3)] {
        let positions: Vec<[f32; 3]> =
            (0..count).map(|_| std::array::from_fn(|_| spread * normal(&mut state))).collect();
        let single: Vec<bool> = (0..count).map(|_| uniform(&mut state) < share).collect();
        let layout = Layout {
            positions,
            chains: chains_with_singles((0..count).map(|index| (index / 150) as i64), &single),
            ligand: single.iter().map(|&value| u8::from(value)).collect(),
            usable: single.iter().map(|&value| u8::from(value || uniform(&mut state) < 0.8)).collect(),
            eligible: None,
            policy: None,
        };
        for size in [50, 384, 512] {
            let (accepted, over) = check_sweep(&layout, size);
            assert!(accepted + over > 0);
            over_count += over;
        }
    }
    assert!(over_count > 20_000, "only {over_count} centres over the limit");
}

/// A filter that recounts every neighbour for each centre; `View::filter` must agree with it.
fn reference_selected(view: &View, initial: &[(f32, usize)]) -> (Vec<usize>, Vec<usize>) {
    let mut counts = vec![0; view.chain_sizes.len()];
    for &(_, index) in initial {
        counts[view.chains[index]] += 1;
    }
    let selected: Vec<usize> = initial
        .iter()
        .map(|&(_, index)| index)
        .filter(|&index| {
            let chain = view.chains[index];
            if view.chain_ligand[chain] {
                counts[chain] == view.chain_sizes[chain]
            } else {
                counts[chain] >= view.chain_minimum[chain]
            }
        })
        .collect();
    (selected, counts)
}

fn reference_filter(view: &View, initial: &[(f32, usize)]) -> (Option<Vec<u64>>, bool) {
    let (selected, mut counts) = reference_selected(view, initial);
    let mut ligand_count = 0;
    let mut ligand_chain = vec![false; counts.len()];
    for &index in &selected {
        if view.ligand[index] {
            ligand_count += 1;
            ligand_chain[view.chains[index]] = true;
        }
    }
    let ligand_chains = ligand_chain.into_iter().filter(|&value| value).count();
    if ligand_count > 64 || ligand_chains > 9 {
        return (None, true);
    }
    counts.fill(0);
    let mut protein = vec![false; counts.len()];
    for &index in &selected {
        if view.usable[index] {
            counts[view.chains[index]] += 1;
            protein[view.chains[index]] |= view.chain_protein[view.chains[index]];
        }
    }
    let retained: Vec<bool> = counts
        .iter()
        .enumerate()
        .map(|(chain, &count)| count > 0 && count >= view.chain_minimum[chain])
        .collect();
    if !retained
        .iter()
        .zip(&protein)
        .any(|(&keep, &protein)| keep && protein)
    {
        return (None, false);
    }
    let mut bits = vec![0; view.chain_sizes.len().div_ceil(64)];
    for (chain, keep) in retained.into_iter().enumerate() {
        if keep {
            bits[chain / 64] |= 1_u64 << (chain % 64);
        }
    }
    (Some(bits), false)
}

#[test]
fn chain_tallies_match_the_reference_filter() {
    let mut state = 23;
    let mut outcomes = [0; 3];
    for trial in 0..400 {
        let count = 200 + trial % 300;
        let chain_count = 2 + trial % 40;
        let chain_ids: Vec<i64> = (0..count)
            .map(|_| (uniform(&mut state) * chain_count as f32) as i64)
            .collect();
        let ligand_chain: Vec<bool> = (0..chain_count)
            .map(|_| uniform(&mut state) < 0.6)
            .collect();
        let ligand: Vec<u8> = chain_ids
            .iter()
            .map(|&chain| u8::from(ligand_chain[chain as usize]))
            .collect();
        let usable: Vec<u8> = (0..count)
            .map(|_| u8::from(uniform(&mut state) > 0.1))
            .collect();
        let bytes: Vec<u8> = (0..count)
            .flat_map(|_| [0_f32; 3])
            .flat_map(f32::to_ne_bytes)
            .collect();
        let chains: Vec<u8> = chain_ids
            .iter()
            .flat_map(|value| value.to_ne_bytes())
            .collect();
        let size = 50 + trial % 150;
        let mut view = View::read(&bytes, &chains, &ligand, &usable, size, None).unwrap();
        for chain in 0..view.chain_minimum.len() {
            if !view.chain_ligand[chain] {
                view.chain_minimum[chain] = 1 + (uniform(&mut state) * 12.) as usize;
            }
        }
        let mut tally = Tally::new(view.chain_sizes.len());
        for _ in 0..20 {
            // Order whole chains randomly so that crops often hold complete ligand chains.
            let priority: Vec<f32> = (0..chain_count).map(|_| uniform(&mut state)).collect();
            let mut initial: Vec<(f32, usize)> = (0..count)
                .map(|index| (priority[chain_ids[index] as usize], index))
                .collect();
            initial.sort_by(compare);
            initial.truncate(size.min(count));
            let expected = reference_filter(&view, &initial);
            outcomes[usize::from(expected.1) + 2 * usize::from(expected.0.is_some())] += 1;
            assert_eq!(view.filter(&initial, &mut tally), expected);
        }
    }
    assert!(outcomes.iter().all(|&seen| seen > 100), "{outcomes:?}");
}
