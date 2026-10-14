//! Exact finite spatial-crop plans over caller-owned token arrays.

use std::cmp::Ordering;
use std::collections::{HashMap, VecDeque};

pub(crate) struct CropPlans {
    pub centers: Vec<i64>,
    pub radii: Vec<f32>,
    pub ties: Vec<i64>,
    pub bits: Vec<u64>,
    pub probabilities: Vec<f64>,
    pub words: usize,
    pub mass: f64,
}

struct View {
    positions: Vec<[f32; 3]>,
    chains: Vec<usize>,
    ligand: Vec<bool>,
    usable: Vec<bool>,
    center_eligible: Vec<bool>,
    chain_sizes: Vec<usize>,
    chain_ligand: Vec<bool>,
    size: usize,
}

fn centers(bytes: &[u8]) -> Result<Vec<[f32; 3]>, String> {
    let result: Vec<[f32; 3]> = bytes
        .chunks_exact(12)
        .map(|row| {
            std::array::from_fn(|axis| {
                let offset = axis * 4;
                f32::from_ne_bytes([
                    row[offset],
                    row[offset + 1],
                    row[offset + 2],
                    row[offset + 3],
                ])
            })
        })
        .collect();
    if result.iter().flatten().any(|value| !value.is_finite()) {
        return Err("Spatial-crop centers must be finite".into());
    }
    Ok(result)
}

impl View {
    fn read(
        positions: &[u8],
        chains: &[u8],
        ligand: &[u8],
        usable: &[u8],
        size: usize,
        center_eligible: Option<&[u8]>,
    ) -> Result<Self, String> {
        let count = ligand.len();
        if !(50..=512).contains(&size)
            || positions.len() / 12 != count
            || !positions.len().is_multiple_of(12)
            || chains.len() / 8 != count
            || !chains.len().is_multiple_of(8)
            || usable.len() != count
            || ligand.iter().chain(usable).any(|&value| value > 1)
            || center_eligible.is_some_and(|values| values.len() != count || values.iter().any(|&value| value > 1))
        {
            return Err("Invalid spatial-crop input shapes, masks or size".into());
        }
        let positions = centers(positions)?;
        let chains: Vec<usize> = chains
            .chunks_exact(8)
            .map(|row| i64::from_ne_bytes(std::array::from_fn(|index| row[index])))
            .map(|value| {
                usize::try_from(value).map_err(|_| "Spatial-crop chain indices must be nonnegative")
            })
            .collect::<Result<_, _>>()?;
        if chains.iter().any(|&chain| chain >= count) {
            return Err("Spatial-crop chain indices exceed the token count".into());
        }
        let chain_count = chains.iter().max().map_or(0, |value| value + 1);
        let ligand: Vec<bool> = ligand.iter().map(|&value| value != 0).collect();
        let mut chain_sizes = vec![0; chain_count];
        let mut chain_ligand = vec![false; chain_count];
        for (index, &chain) in chains.iter().enumerate() {
            chain_sizes[chain] += 1;
            chain_ligand[chain] |= ligand[index];
        }
        Ok(Self {
            positions,
            chains,
            ligand,
            usable: usable.iter().map(|&value| value != 0).collect(),
            center_eligible: center_eligible.map_or_else(|| vec![true; count], |values| values.iter().map(|&value| value != 0).collect()),
            chain_sizes,
            chain_ligand,
            size,
        })
    }

    fn nearest(&self, center: usize) -> Vec<(f32, usize)> {
        let point = self.positions[center];
        let mut distances: Vec<_> = self
            .positions
            .iter()
            .enumerate()
            .map(|(index, other)| {
                let delta = std::array::from_fn::<_, 3, _>(|axis| other[axis] - point[axis]);
                (
                    (delta[0] * delta[0] + delta[1] * delta[1]) + delta[2] * delta[2],
                    index,
                )
            })
            .collect();
        let count = distances.len().min(self.size);
        if count < distances.len() {
            distances.select_nth_unstable_by(count - 1, compare);
            distances.truncate(count);
        } else if let Some((boundary, _)) = distances
            .iter()
            .enumerate()
            .max_by(|left, right| compare(left.1, right.1))
        {
            distances.swap(boundary, count - 1);
        }
        // Filtering uses counts, so only the last boundary pair needs ordering.
        distances
    }

    fn selected(&self, initial: &[(f32, usize)]) -> (Vec<usize>, Vec<usize>) {
        let mut counts = vec![0; self.chain_sizes.len()];
        for &(_, index) in initial {
            counts[self.chains[index]] += 1;
        }
        let selected: Vec<usize> = initial
            .iter()
            .map(|&(_, index)| index)
            .filter(|&index| {
                let chain = self.chains[index];
                if self.chain_ligand[chain] {
                    counts[chain] == self.chain_sizes[chain]
                } else {
                    counts[chain] >= 50
                }
            })
            .collect();
        (selected, counts)
    }

    fn filter(&self, initial: &[(f32, usize)]) -> (Option<Vec<u64>>, bool) {
        let (selected, mut counts) = self.selected(initial);
        let mut ligand_count = 0;
        let mut ligand_chain = vec![false; counts.len()];
        for &index in &selected {
            if self.ligand[index] {
                ligand_count += 1;
                ligand_chain[self.chains[index]] = true;
            }
        }
        let ligand_chains = ligand_chain.into_iter().filter(|&value| value).count();
        if ligand_count > 64 || ligand_chains > 9 {
            return (None, true);
        }
        counts.fill(0);
        let mut protein = vec![false; counts.len()];
        for &index in &selected {
            if self.usable[index] {
                counts[self.chains[index]] += 1;
                protein[self.chains[index]] |= !self.ligand[index];
            }
        }
        let retained: Vec<bool> = counts
            .iter()
            .zip(&protein)
            .map(|(&count, &protein)| count > 0 && (!protein || count >= 50))
            .collect();
        if !retained
            .iter()
            .zip(&protein)
            .any(|(&keep, &protein)| keep && protein)
        {
            return (None, false);
        }
        let mut bits = vec![0; self.chain_sizes.len().div_ceil(64)];
        for (chain, keep) in retained.into_iter().enumerate() {
            if keep {
                bits[chain / 64] |= 1_u64 << (chain % 64);
            }
        }
        (Some(bits), false)
    }
}

fn compare(left: &(f32, usize), right: &(f32, usize)) -> Ordering {
    left.0
        .total_cmp(&right.0)
        .then_with(|| left.1.cmp(&right.1))
}

pub(crate) fn plan(
    positions: &[u8],
    chains: &[u8],
    ligand: &[u8],
    usable: &[u8],
    size: usize,
    center_eligible: Option<&[u8]>,
) -> Result<CropPlans, String> {
    let view = View::read(positions, chains, ligand, usable, size, center_eligible)?;
    let total = view.center_eligible.iter().filter(|&&value| value).count();
    let ligand_count = view.ligand.iter().zip(&view.center_eligible).filter(|&(ligand, eligible)| *ligand && *eligible).count();
    let mut result = CropPlans {
        centers: Vec::new(),
        radii: Vec::new(),
        ties: Vec::new(),
        bits: Vec::new(),
        probabilities: Vec::new(),
        words: view.chain_sizes.len().div_ceil(64).max(1),
        mass: 0.,
    };
    let mut retry = 0;
    if total == 0 || view.ligand.iter().all(|&value| value) {
        return Ok(result);
    }
    for center in 0..view.positions.len() {
        if !view.center_eligible[center] { continue; }
        let nearest = view.nearest(center);
        let (bits, over) = view.filter(&nearest);
        retry += usize::from(over && view.ligand[center]);
        if let Some(bits) = bits {
            let &(radius, tie) = nearest
                .last()
                .ok_or("An accepted spatial crop has no neighbors")?;
            result.centers.push(center as i64);
            result.radii.push(radius);
            result.ties.push(tie as i64);
            result.bits.extend(bits);
        }
    }
    let base = if ligand_count > 0 {
        (0.01 + 0.99 * retry as f64 / ligand_count as f64) / total as f64
    } else {
        1. / total as f64
    };
    result.probabilities = result
        .centers
        .iter()
        .map(|&center| {
            base + if view.ligand[center as usize] {
                0.99 / ligand_count as f64
            } else {
                0.
            }
        })
        .collect();
    result.mass = result.probabilities.iter().sum();
    if result.mass > 0. {
        result
            .probabilities
            .iter_mut()
            .for_each(|value| *value /= result.mass);
    }
    Ok(result)
}

pub(crate) fn permutation(
    before: &[u8],
    after: &[u8],
    chains: &[u8],
    ligand: &[u8],
) -> Result<Vec<i64>, String> {
    let usable = vec![1; ligand.len()];
    let left = View::read(before, chains, ligand, &usable, 50, None)?;
    let right = View::read(after, chains, ligand, &usable, 50, None)?;
    let bits =
        |position: [f32; 3]| position.map(|value| if value == 0. { 0 } else { value.to_bits() });
    let mut positions: HashMap<_, VecDeque<i64>> = HashMap::new();
    for (index, position) in left.positions.iter().enumerate() {
        if left.ligand[index] {
            positions
                .entry((left.chains[index], bits(*position)))
                .or_default()
                .push_back(index as i64);
        }
    }
    right
        .positions
        .iter()
        .enumerate()
        .map(|(index, position)| {
            if !left.ligand[index] {
                if *position == left.positions[index] {
                    return Ok(index as i64);
                }
                return Err("Protonation moved a protein crop center".into());
            }
            positions
                .get_mut(&(left.chains[index], bits(*position)))
                .and_then(VecDeque::pop_front)
                .ok_or_else(|| "Protonation changed a ligand chain's coordinate set".into())
        })
        .collect()
}

#[cfg(test)]
#[path = "crop_tests.rs"]
mod tests;
