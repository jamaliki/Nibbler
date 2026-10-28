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
    /// Ligand tokens and ligand chains that each accepted crop keeps.
    pub ligand_tokens: Vec<u32>,
    pub ligand_chains: Vec<u32>,
    /// The same counts for the crop of every eligible ligand centre, in token order, whether
    /// or not the crop keeps a protein chain. Only uncapped plans fill these.
    pub centre_ligand_tokens: Vec<u32>,
    pub centre_ligand_chains: Vec<u32>,
}

impl CropPlans {
    fn empty(words: usize) -> Self {
        Self {
            centers: Vec::new(),
            radii: Vec::new(),
            ties: Vec::new(),
            bits: Vec::new(),
            probabilities: Vec::new(),
            words,
            mass: 0.,
            ligand_tokens: Vec::new(),
            ligand_chains: Vec::new(),
            centre_ligand_tokens: Vec::new(),
            centre_ligand_chains: Vec::new(),
        }
    }
}

/// Default ligand limits: a crop that keeps more ligand tokens or ligand chains is rejected.
pub(crate) const LIGAND_LIMITS: (u32, u32) = (64, 9);

struct View {
    positions: Vec<[f32; 3]>,
    chains: Vec<usize>,
    ligand: Vec<bool>,
    usable: Vec<bool>,
    center_eligible: Vec<bool>,
    chain_sizes: Vec<usize>,
    chain_ligand: Vec<bool>,
    chain_minimum: Vec<usize>,
    chain_protein: Vec<bool>,
    size: usize,
    /// Most ligand tokens and ligand chains a crop may keep; `u32::MAX` when uncapped.
    token_limit: u32,
    chain_limit: u32,
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
            || center_eligible.is_some_and(|values| {
                values.len() != count || values.iter().any(|&value| value > 1)
            })
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
        // Complete-residue filtering preserves chain IDs and can leave gaps.
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
            center_eligible: center_eligible.map_or_else(
                || vec![true; count],
                |values| values.iter().map(|&value| value != 0).collect(),
            ),
            chain_sizes,
            chain_minimum: chain_ligand
                .iter()
                .map(|&ligand| if ligand { 1 } else { 50 })
                .collect(),
            chain_protein: chain_ligand.iter().map(|&ligand| !ligand).collect(),
            chain_ligand,
            size,
            token_limit: LIGAND_LIMITS.0,
            chain_limit: LIGAND_LIMITS.1,
        })
    }

    fn capped(&self) -> bool {
        self.token_limit != u32::MAX || self.chain_limit != u32::MAX
    }

    fn over(&self, tokens: u32, chains: u32) -> bool {
        tokens > self.token_limit || chains > self.chain_limit
    }

    /// Write the crop size nearest tokens under `compare` into `distances`, with the farthest
    /// last.
    ///
    /// Other entries stay unordered because filtering only counts them.
    fn nearest_all(&self, center: usize, distances: &mut Vec<(f32, usize)>) {
        let point = self.positions[center];
        distances.clear();
        distances.extend(
            self.positions
                .iter()
                .enumerate()
                .map(|(index, other)| (distance(*other, point), index)),
        );
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
    }

    /// Return the retained chain bits, or None, whether the crop has too many ligands, and the
    /// ligand tokens and ligand chains it keeps.
    ///
    /// A neighbour stays when its chain does: a ligand chain must be complete and any other
    /// chain needs its minimum count. `tally` must be zero on entry and is zero on return.
    fn filter(&self, initial: &[(f32, usize)], tally: &mut Tally) -> Filtered {
        for &(_, index) in initial {
            let chain = self.chains[index];
            if tally.tokens[chain] == 0 {
                tally.chains.push(chain);
            }
            tally.tokens[chain] += 1;
            tally.ligands[chain] += usize::from(self.ligand[index]);
            tally.usable[chain] += usize::from(self.usable[index]);
        }
        let (mut ligand_count, mut ligand_chains, mut protein) = (0_u32, 0_u32, false);
        let mut bits = vec![0; self.chain_sizes.len().div_ceil(64)];
        for &chain in &tally.chains {
            let count = tally.tokens[chain];
            let kept = if self.chain_ligand[chain] {
                count == self.chain_sizes[chain]
            } else {
                count >= self.chain_minimum[chain]
            };
            if kept {
                ligand_count += tally.ligands[chain] as u32;
                ligand_chains += u32::from(tally.ligands[chain] > 0);
                let usable = tally.usable[chain];
                if usable > 0 && usable >= self.chain_minimum[chain] {
                    bits[chain / 64] |= 1_u64 << (chain % 64);
                    protein |= self.chain_protein[chain];
                }
            }
            tally.tokens[chain] = 0;
            tally.ligands[chain] = 0;
            tally.usable[chain] = 0;
        }
        tally.chains.clear();
        let ligands = (ligand_count, ligand_chains);
        if self.over(ligand_count, ligand_chains) {
            return (None, true, ligands);
        }
        (protein.then_some(bits), false, ligands)
    }
}

/// Retained chain bits or None, the ligand-limit flag, and the kept ligand tokens and chains.
type Filtered = (Option<Vec<u64>>, bool, (u32, u32));

/// Per-chain counts reused across crop centres.
struct Tally {
    tokens: Vec<usize>,
    ligands: Vec<usize>,
    usable: Vec<usize>,
    chains: Vec<usize>,
}

impl Tally {
    fn new(chains: usize) -> Self {
        Self {
            tokens: vec![0; chains],
            ligands: vec![0; chains],
            usable: vec![0; chains],
            chains: Vec::new(),
        }
    }
}

/// Filter result and the farthest neighbour of one centre.
type Crop = (Filtered, Option<(f32, usize)>);

/// Coordinates as a hash key; signed zeros give equal distances, so they share a key.
fn site_key(point: [f32; 3]) -> [u32; 3] {
    point.map(|value| if value == 0. { 0 } else { value.to_bits() })
}

/// Squared distance with the operation order that every crop radius uses.
fn distance(other: [f32; 3], point: [f32; 3]) -> f32 {
    let delta = std::array::from_fn::<_, 3, _>(|axis| other[axis] - point[axis]);
    (delta[0] * delta[0] + delta[1] * delta[1]) + delta[2] * delta[2]
}

/// Fewest tokens for which the sweep replaces a scan over every token per centre. The
/// sweep also needs at least four times the crop size.
const SWEEP_MIN_TOKENS: usize = 2048;
/// Centres whose crop radius sets the sweep cell size.
const SWEEP_SAMPLES: usize = 16;
/// Cells per crop radius when sizing sweep cells.
const SWEEP_CELL_DIVISOR: f64 = 4.;
/// Sweep blocks span this many cells along each axis.
const BLOCK_CELLS: usize = 2;
/// Largest selection band left to `select_nth_unstable`.
const BAND_TOKENS: usize = 48;
/// Counting passes before the band is taken as it stands.
const BAND_PASSES: usize = 8;
/// Share of the previous squared crop radius used as the first counting threshold.
const GUESS_SHARE: f32 = 0.9;
/// Surveys of one block, each with a longer reach, before the block is searched in full.
const SKIP_TRIES: usize = 3;

/// Neighbour tokens sorted by grid cell with z fastest.
///
/// Each site keeps only its `size` lowest token indices. At every centre those tokens sort
/// before the rest of the site, so no crop can reach the others.
struct Cells {
    origin: [f64; 3],
    cell: f64,
    dims: [usize; 3],
    /// Start of each cell's tokens, then the token count.
    starts: Vec<u32>,
    /// Coordinates of each token, then its index as f32 bits.
    points: Vec<[f32; 4]>,
    /// Counts over the tokens before each entry of `points`, then over all of them: tokens
    /// that are usable and in a protein chain, ligand tokens, and entries of `singles`.
    prefix: Vec<[u32; 3]>,
    /// Coordinates of the tokens that form a ligand chain on their own, in `points` order.
    singles: Vec<[f32; 3]>,
    /// Median crop radius over the sampled centres.
    typical: f64,
}

impl Cells {
    fn build(view: &View) -> Option<Self> {
        let (positions, size) = (&view.positions, view.size);
        let count = positions.len();
        if count <= size || u32::try_from(count).is_err() {
            return None;
        }
        let mut origin = [f64::INFINITY; 3];
        let mut extent = [f64::NEG_INFINITY; 3];
        for point in positions {
            for axis in 0..3 {
                origin[axis] = origin[axis].min(f64::from(point[axis]));
                extent[axis] = extent[axis].max(f64::from(point[axis]));
            }
        }
        let extent: [f64; 3] = std::array::from_fn(|axis| extent[axis] - origin[axis]);
        let largest = extent.iter().copied().fold(0., f64::max);
        if !largest.is_finite() || largest <= 0. {
            return None;
        }
        let mut radii: Vec<f64> = (0..SWEEP_SAMPLES)
            .map(|sample| {
                let point = positions[sample * count / SWEEP_SAMPLES];
                let mut distances: Vec<f32> = positions
                    .iter()
                    .map(|&other| distance(other, point))
                    .collect();
                distances.select_nth_unstable_by(size - 1, f32::total_cmp);
                f64::from(distances[size - 1]).sqrt()
            })
            .filter(|&radius| radius > 0.)
            .collect();
        radii.sort_by(f64::total_cmp);
        let typical = radii.get(radii.len() / 2).copied().unwrap_or(largest);
        let mut cell = (typical / SWEEP_CELL_DIVISOR).max(largest / 256.);
        let dims = loop {
            let dims: [usize; 3] = std::array::from_fn(|axis| (extent[axis] / cell) as usize + 1);
            if dims.iter().product::<usize>() <= 8 * count {
                break dims;
            }
            cell *= 1.25;
        };
        let mut cells = Self {
            origin,
            cell,
            dims,
            starts: Vec::new(),
            points: Vec::with_capacity(count),
            prefix: Vec::with_capacity(count + 1),
            singles: Vec::new(),
            typical,
        };
        let mut order: Vec<(usize, [u32; 3], usize)> = positions
            .iter()
            .enumerate()
            .map(|(token, &point)| (cells.index(cells.home(point)), site_key(point), token))
            .collect();
        order.sort_unstable();
        let mut starts = vec![0_u32; dims.iter().product::<usize>() + 1];
        let (mut counts, mut repeat) = ([0; 3], 0);
        for (row, &(cell, site, token)) in order.iter().enumerate() {
            repeat = if row > 0 && order[row - 1].0 == cell && order[row - 1].1 == site {
                repeat + 1
            } else {
                0
            };
            if repeat < size {
                let [x, y, z] = positions[token];
                let chain = view.chains[token];
                cells.points.push([x, y, z, f32::from_bits(token as u32)]);
                cells.prefix.push(counts);
                starts[cell + 1] += 1;
                let single = view.ligand[token] && view.chain_sizes[chain] == 1;
                if single {
                    cells.singles.push([x, y, z]);
                }
                counts[0] += u32::from(view.usable[token] && view.chain_protein[chain]);
                counts[1] += u32::from(view.ligand[token]);
                counts[2] += u32::from(single);
            }
        }
        cells.prefix.push(counts);
        for cell in 1..starts.len() {
            starts[cell] += starts[cell - 1];
        }
        cells.starts = starts;
        Some(cells)
    }

    fn home(&self, point: [f32; 3]) -> [usize; 3] {
        std::array::from_fn(|axis| {
            (((f64::from(point[axis]) - self.origin[axis]) / self.cell) as i64)
                .clamp(0, self.dims[axis] as i64 - 1) as usize
        })
    }

    fn index(&self, home: [usize; 3]) -> usize {
        (home[0] * self.dims[1] + home[1]) * self.dims[2] + home[2]
    }
}

/// What one centre's crop gives the plan.
#[derive(Clone, Copy)]
struct Outcome {
    over: bool,
    /// Start of the chain words in the shared store, or `u32::MAX` when the crop is rejected.
    bits: u32,
    radius: f32,
    tie: u32,
    /// Ligand tokens and ligand chains the crop keeps, when it was measured.
    ligands: (u32, u32),
}

/// Contribution of one chain to the filter of the current crop.
#[derive(Clone, Copy, Default, PartialEq)]
struct Share {
    ligands: u32,
    retained: bool,
}

/// Chain tallies of the current crop, kept up to date one token at a time.
///
/// The sums match `View::filter`: a chain stays if a ligand chain is complete or another
/// chain reaches its minimum, and a kept chain is retained if its usable count reaches the
/// minimum too.
struct Tallies {
    chains: Vec<u32>,
    /// Per token: one, the ligand flag and the usable flag, in 16-bit fields.
    weights: Vec<u64>,
    /// Per chain sums of `weights` over the crop.
    counts: Vec<u64>,
    shares: Vec<Share>,
    dirty: Vec<u32>,
    is_dirty: Vec<bool>,
    touched: Vec<u32>,
    is_touched: Vec<bool>,
    ligand_count: u32,
    ligand_chains: u32,
    protein: u32,
    words: Vec<u64>,
}

impl Tallies {
    fn new(view: &View) -> Self {
        let chain_count = view.chain_sizes.len();
        Self {
            chains: view.chains.iter().map(|&chain| chain as u32).collect(),
            weights: view
                .ligand
                .iter()
                .zip(&view.usable)
                .map(|(&ligand, &usable)| 1 | u64::from(ligand) << 16 | u64::from(usable) << 32)
                .collect(),
            counts: vec![0; chain_count],
            shares: vec![Share::default(); chain_count],
            dirty: Vec::new(),
            is_dirty: vec![false; chain_count],
            touched: Vec::new(),
            is_touched: vec![false; chain_count],
            ligand_count: 0,
            ligand_chains: 0,
            protein: 0,
            words: vec![0; chain_count.div_ceil(64)],
        }
    }

    fn reset(&mut self) {
        for &chain in &self.touched {
            let chain = chain as usize;
            self.counts[chain] = 0;
            self.shares[chain] = Share::default();
            self.is_touched[chain] = false;
        }
        self.touched.clear();
        self.ligand_count = 0;
        self.ligand_chains = 0;
        self.protein = 0;
        self.words.fill(0);
    }

    /// Add a token to the crop or remove it.
    fn change(&mut self, token: usize, inside: bool) {
        let chain = self.chains[token] as usize;
        let weight = self.weights[token];
        self.counts[chain] = if inside {
            self.counts[chain] + weight
        } else {
            self.counts[chain] - weight
        };
        if !self.is_dirty[chain] {
            self.is_dirty[chain] = true;
            self.dirty.push(chain as u32);
        }
    }

    /// Refresh the shares of changed chains; return the ligand-limit flag and acceptance.
    fn settle(&mut self, view: &View) -> (bool, bool) {
        for slot in 0..self.dirty.len() {
            let chain = self.dirty[slot] as usize;
            self.is_dirty[chain] = false;
            if !self.is_touched[chain] {
                self.is_touched[chain] = true;
                self.touched.push(chain as u32);
            }
            let count = self.counts[chain];
            let (tokens, ligands, usable) = (
                (count & 0xffff) as usize,
                ((count >> 16) & 0xffff) as u32,
                (count >> 32) as usize,
            );
            let kept = if view.chain_ligand[chain] {
                tokens == view.chain_sizes[chain]
            } else {
                tokens >= view.chain_minimum[chain]
            };
            let share = if kept {
                Share {
                    ligands,
                    retained: usable > 0 && usable >= view.chain_minimum[chain],
                }
            } else {
                Share::default()
            };
            let old = self.shares[chain];
            if share != old {
                self.shares[chain] = share;
                self.ligand_count = self.ligand_count + share.ligands - old.ligands;
                self.ligand_chains =
                    self.ligand_chains + u32::from(share.ligands > 0) - u32::from(old.ligands > 0);
                if share.retained != old.retained {
                    self.words[chain / 64] ^= 1_u64 << (chain % 64);
                    if view.chain_protein[chain] {
                        self.protein =
                            self.protein + u32::from(share.retained) - u32::from(old.retained);
                    }
                }
            }
        }
        self.dirty.clear();
        let over = view.over(self.ligand_count, self.ligand_chains);
        (over, !over && self.protein > 0)
    }
}

/// Call `visit` with the index of every nonzero flag, in increasing order.
///
/// Flags are read eight at a time as one word, so that runs of zeros cost one test.
fn each_flag(flags: &[u8], mut visit: impl FnMut(usize)) {
    let chunks = flags.chunks_exact(8);
    let tail = chunks.remainder();
    for (chunk, bytes) in chunks.enumerate() {
        let Ok(bytes) = <[u8; 8]>::try_from(bytes) else {
            continue;
        };
        let mut word = u64::from_le_bytes(bytes);
        while word != 0 {
            let offset = word.trailing_zeros() as usize / 8;
            word &= word - 1;
            visit(8 * chunk + offset);
        }
    }
    let base = flags.len() - tail.len();
    for (offset, &flag) in tail.iter().enumerate() {
        if flag != 0 {
            visit(base + offset);
        }
    }
}

/// Result of one boundary search.
enum Search {
    /// The gathered tokens may not hold the crop.
    Miss,
    /// The crop keeps more ligand chains than the limit; the value estimates its radius.
    Over(f32),
    /// Key of the crop boundary.
    Found(u64),
}

/// Threshold at which the count should reach `size`, between two counted thresholds.
fn interpolate(size: usize, (low, below): (u32, usize), (high, within): (u32, usize)) -> f64 {
    // Counts grow with distance to the power 1.5 when density is even.
    let (start, end) = (
        f64::from(f32::from_bits(low)).powf(1.5),
        f64::from(f32::from_bits(high)).powf(1.5),
    );
    let share = (size - below) as f64 / (within - below) as f64;
    (start + share * (end - start)).powf(2. / 3.)
}

/// Crop search that visits centres block by block.
///
/// All centres in a block share one gathered list of nearby tokens. For each centre the
/// sweep computes every distance in that list, finds the crop boundary by counting tokens
/// under trial thresholds, and updates the chain tallies only for tokens whose membership
/// changed since the previous centre it measured.
///
/// Two tests skip that work when they prove a centre is rejected. A block whose gathered
/// tokens surely hold every crop of its centres, but hold too few protein-usable tokens for
/// any protein chain, rejects all of its centres. A crop whose nearest tokens already
/// include more single-token ligand chains than the limit is rejected as over the limit.
struct Sweep<'a> {
    view: &'a View,
    cells: Cells,
    /// Gathered tokens of the current block.
    xs: Vec<f32>,
    ys: Vec<f32>,
    zs: Vec<f32>,
    tokens: Vec<u32>,
    /// Ranges of `cells.points` and `cells.singles` that the current block gathers.
    ranges: Vec<(usize, usize)>,
    single_ranges: Vec<(usize, usize)>,
    /// Gathered single-token ligand chains.
    singles: Vec<[f32; 3]>,
    /// Squared distances from the current centre, as f32 bits.
    distances: Vec<u32>,
    /// Largest of `distances`.
    farthest: u32,
    /// Crop membership of the gathered tokens, as 0 or 1.
    member: Vec<u8>,
    flags: Vec<u8>,
    band: Vec<u64>,
    /// Largest f32 squared distance that the gathered tokens cover, or None if none.
    cover: Option<u32>,
    /// Smallest minimum of a protein chain; a crop with fewer protein-usable tokens is
    /// rejected.
    protein_floor: u32,
    tallies: Tallies,
}

impl<'a> Sweep<'a> {
    fn new(view: &'a View, cells: Cells) -> Self {
        let protein_floor = view
            .chain_minimum
            .iter()
            .zip(&view.chain_protein)
            .filter(|&(_, &protein)| protein)
            .map(|(&minimum, _)| u32::try_from(minimum).unwrap_or(u32::MAX))
            .min()
            .unwrap_or(u32::MAX);
        Self {
            view,
            cells,
            xs: Vec::new(),
            ys: Vec::new(),
            zs: Vec::new(),
            tokens: Vec::new(),
            ranges: Vec::new(),
            single_ranges: Vec::new(),
            singles: Vec::new(),
            distances: Vec::new(),
            farthest: 0,
            member: Vec::new(),
            flags: Vec::new(),
            band: Vec::new(),
            cover: None,
            protein_floor,
            tallies: Tallies::new(view),
        }
    }

    /// Find the cells within `reach` of the block; return their protein-usable and ligand
    /// token counts.
    fn survey(&mut self, low: [usize; 3], high: [usize; 3], reach: f64) -> (u32, u32) {
        let cells = &self.cells;
        // Cells this many steps past the block are at least (steps - 1) cells away; the 1e-3
        // allowance covers round-off when binning coordinates.
        let gap = |step: usize| ((step as f64 - 1.) * cells.cell - 1e-3).max(0.);
        let span = |budget: f64| ((budget.max(0.).sqrt() + 1e-3) / cells.cell) as usize + 1;
        let budget = reach * reach;
        let range = |axis: usize, steps: usize| {
            (
                low[axis].saturating_sub(steps),
                (high[axis] + steps).min(cells.dims[axis] - 1),
            )
        };
        let (x0, x1) = range(0, span(budget));
        let (y0, y1) = range(1, span(budget));
        // Each column adds at most one range. Writing every column and keeping the nonempty
        // ones avoids a branch that empty columns would often mispredict.
        self.ranges.clear();
        self.ranges.resize((x1 + 1 - x0) * (y1 + 1 - y0), (0, 0));
        let mut used = 0;
        for x in x0..=x1 {
            let gx = gap(low[0].saturating_sub(x).max(x.saturating_sub(high[0])));
            let left = budget - gx * gx;
            if left < 0. {
                continue;
            }
            let (y0, y1) = range(1, span(left));
            for y in y0..=y1 {
                let gy = gap(low[1].saturating_sub(y).max(y.saturating_sub(high[1])));
                let rest = left - gy * gy;
                if rest < 0. {
                    continue;
                }
                let (z0, z1) = range(2, span(rest));
                let (from, to) = (cells.index([x, y, z0]), cells.index([x, y, z1]) + 1);
                let (first, last) = (cells.starts[from] as usize, cells.starts[to] as usize);
                self.ranges[used] = (first, last);
                used += usize::from(first < last);
            }
        }
        self.ranges.truncate(used);
        let (mut total, mut protein, mut ligands) = (0, 0, 0);
        self.single_ranges.clear();
        for &(first, last) in &self.ranges {
            let (before, after) = (cells.prefix[first], cells.prefix[last]);
            total += last - first;
            protein += after[0] - before[0];
            ligands += after[1] - before[1];
            self.single_ranges
                .push((before[2] as usize, after[2] as usize));
        }
        self.cover = if total == cells.points.len() {
            Some(u32::MAX)
        } else {
            // A computed f32 squared distance is at least (1 - 3e-7) times the true one, less
            // any underflow, so a token under this bound lies within `reach` of the centre.
            let bound = budget * (1. - 1e-5) - 1e-30;
            let value = bound as f32;
            let value = if f64::from(value) > bound {
                f32::from_bits(value.to_bits().saturating_sub(1))
            } else {
                value
            };
            (bound >= 0.).then_some(value.to_bits())
        };
        (protein, ligands)
    }

    /// Count the tokens that lie within the cover of every point in the block.
    ///
    /// If the count reaches the crop size, no crop of a centre in the block reaches past
    /// the cover, so the gathered tokens hold it.
    fn inner(&self, low: [usize; 3], high: [usize; 3]) -> usize {
        let cells = &self.cells;
        let cover = match self.cover {
            None => return 0,
            Some(u32::MAX) => return cells.points.len(),
            Some(cover) => f64::from(f32::from_bits(cover)),
        };
        // A token this far away has a computed squared distance under the cover, by the
        // round-off bound in `survey`.
        let budget = (cover - 1e-29) / (1. + 1e-4);
        // Along one axis, a token in cell `at` lies at most this far from any point in the
        // block, with the binning allowance.
        let far = |axis: usize, at: usize| {
            (at + 1)
                .saturating_sub(low[axis])
                .max((high[axis] + 1).saturating_sub(at)) as f64
                * cells.cell
                + 1e-3
        };
        // Cells along one axis whose farthest point stays within sqrt(budget).
        let range = |axis: usize, budget: f64| {
            let steps = ((budget.max(0.).sqrt() - 2e-3) / cells.cell).floor();
            let width = (high[axis] + 1 - low[axis]) as f64;
            (steps >= width).then(|| {
                let steps = steps as usize;
                (
                    (high[axis] + 1).saturating_sub(steps),
                    (low[axis] + steps - 1).min(cells.dims[axis] - 1),
                )
            })
        };
        let mut count = 0;
        let Some((x0, x1)) = range(0, budget) else {
            return 0;
        };
        for x in x0..=x1 {
            let left = budget - far(0, x).powi(2);
            let Some((y0, y1)) = range(1, left) else {
                continue;
            };
            for y in y0..=y1 {
                let Some((z0, z1)) = range(2, left - far(1, y).powi(2)) else {
                    continue;
                };
                count += (cells.starts[cells.index([x, y, z1]) + 1]
                    - cells.starts[cells.index([x, y, z0])]) as usize;
            }
        }
        count
    }

    /// Copy the surveyed tokens and reset the crop.
    fn fill(&mut self) {
        let cells = &self.cells;
        self.xs.clear();
        self.ys.clear();
        self.zs.clear();
        self.tokens.clear();
        for &(first, last) in &self.ranges {
            for point in &cells.points[first..last] {
                self.xs.push(point[0]);
                self.ys.push(point[1]);
                self.zs.push(point[2]);
                self.tokens.push(point[3].to_bits());
            }
        }
        let total = self.tokens.len();
        self.singles.clear();
        for &(first, last) in &self.single_ranges {
            self.singles.extend(&cells.singles[first..last]);
        }
        self.member.clear();
        self.member.resize(total, 0);
        self.tallies.reset();
    }

    /// Compute the distance of every gathered token from `point`; return how many lie at or
    /// under each threshold.
    fn measure(&mut self, point: [f32; 3], first: u32, second: u32) -> (usize, usize) {
        self.distances.resize(self.tokens.len(), 0);
        let (mut farthest, mut low, mut high) = (0_u32, 0_u32, 0_u32);
        let others = self.xs.iter().zip(&self.ys).zip(&self.zs);
        for (value, ((&x, &y), &z)) in self.distances.iter_mut().zip(others) {
            let bits = distance([x, y, z], point).to_bits();
            *value = bits;
            farthest = farthest.max(bits);
            low += u32::from(bits <= first);
            high += u32::from(bits <= second);
        }
        self.farthest = farthest;
        (low as usize, high as usize)
    }

    /// Count distances at or under each threshold.
    fn count(&self, first: u32, second: u32) -> (usize, usize) {
        let (mut low, mut high) = (0_u32, 0_u32);
        for &value in &self.distances {
            low += u32::from(value <= first);
            high += u32::from(value <= second);
        }
        (low as usize, high as usize)
    }

    /// Whether more single-token ligand chains than the limit lie at or under `limit`.
    ///
    /// The caller ensures that fewer tokens than the crop size lie at or under `limit`.
    /// Those tokens are then all in the crop, and each of these chains is complete in it.
    fn over(&self, point: [f32; 3], limit: u32) -> bool {
        self.view.chain_limit != u32::MAX
            && self
                .singles
                .iter()
                .filter(|&&other| distance(other, point).to_bits() <= limit)
                .count()
                > self.view.chain_limit as usize
    }

    /// Find the crop boundary of the centre at `point`.
    ///
    /// A key holds a distance's bits above its token index; distances are sums of squares,
    /// never negative or NaN, so integer order matches `compare`.
    fn search(&mut self, point: [f32; 3], guess: u32) -> Search {
        let size = self.view.size;
        let Some(cover) = self.cover else {
            return Search::Miss;
        };
        let guess = guess.min(cover);
        let (at_guess, at_cover) = self.measure(point, guess, cover);
        // No distance exceeds `farthest`, so these thresholds keep their counts.
        let cover = cover.min(self.farthest);
        let guess = guess.min(cover);
        if at_cover < size {
            return Search::Miss;
        }
        if at_guess < size && self.over(point, guess) {
            let estimate = interpolate(size, (guess, at_guess), (cover, at_cover));
            return Search::Over(estimate as f32);
        }
        // The boundary distance lies in [low, high]: `below` distances are under `low` and
        // `within` are at or under `high`.
        let (mut low, mut below, mut high, mut within) = if at_guess >= size {
            (0, 0, guess, at_guess)
        } else {
            (guess + 1, at_guess, cover, at_cover)
        };
        for _ in 0..BAND_PASSES {
            if within - below <= BAND_TOKENS || low >= high {
                break;
            }
            let target = interpolate(size, (low, below), (high, within));
            let spread = 0.03 * target;
            let pick = |value: f64| ((value.max(0.) as f32).to_bits()).clamp(low, high - 1);
            let first = pick(target - spread);
            let second = pick(target + spread).max(first);
            let (at_first, at_second) = self.count(first, second);
            if at_first >= size {
                (high, within) = (first, at_first);
            } else if at_second >= size {
                (low, below, high, within) = (first + 1, at_first, second, at_second);
            } else {
                (low, below) = (second + 1, at_second);
            }
        }
        // Fewer than `size` distances lie under `low`.
        if low > 0 && (at_guess >= size || low > guess + 1) && self.over(point, low - 1) {
            return Search::Over(f32::from_bits(low));
        }
        self.flags.clear();
        self.flags.extend(
            self.distances
                .iter()
                .map(|&value| u8::from((value >= low) & (value <= high))),
        );
        self.band.clear();
        let (distances, tokens, band) = (&self.distances, &self.tokens, &mut self.band);
        each_flag(&self.flags, |index| {
            band.push(u64::from(distances[index]) << 32 | u64::from(tokens[index]));
        });
        let rank = size - below - 1;
        self.band.select_nth_unstable(rank);
        Search::Found(self.band[rank])
    }

    /// Move the crop to the tokens at or under `boundary` and return its filter result.
    fn update(&mut self, boundary: u64) -> (bool, bool) {
        let (distance, tie) = ((boundary >> 32) as u32, boundary as u32);
        self.flags.clear();
        self.flags.extend(
            self.distances
                .iter()
                .zip(&self.tokens)
                .zip(&self.member)
                .map(|((&value, &token), &member)| {
                    u8::from((value < distance) | ((value == distance) & (token <= tie))) ^ member
                }),
        );
        let (member, tokens, tallies) = (&mut self.member, &self.tokens, &mut self.tallies);
        each_flag(&self.flags, |index| {
            member[index] ^= 1;
            tallies.change(tokens[index] as usize, member[index] != 0);
        });
        self.tallies.settle(self.view)
    }

    fn run(mut self) -> (Vec<Outcome>, Vec<u64>) {
        let view = self.view;
        let size = view.size;
        let count = view.positions.len();
        let rejected = Outcome {
            over: false,
            bits: u32::MAX,
            radius: 0.,
            tie: 0,
            ligands: (0, 0),
        };
        // A block's gathered ligand tokens bound the ligand tokens and chains of its crops, so
        // few enough of them prove that no crop is over a limit. Uncapped plans measure every
        // ligand centre to report its counts.
        let ligand_bound = if view.capped() {
            view.token_limit.min(view.chain_limit)
        } else {
            0
        };
        let mut outcomes = vec![rejected; count];
        let mut store = Vec::new();
        let blocks: [usize; 3] =
            std::array::from_fn(|axis| self.cells.dims[axis].div_ceil(BLOCK_CELLS));
        // Visit blocks in a serpentine order so that consecutive blocks touch.
        let rank = |home: [usize; 3]| -> usize {
            let block: [usize; 3] = std::array::from_fn(|axis| home[axis] / BLOCK_CELLS);
            let y = if block[0].is_multiple_of(2) {
                block[1]
            } else {
                blocks[1] - 1 - block[1]
            };
            let row = block[0] * blocks[1] + y;
            let z = if row.is_multiple_of(2) {
                block[2]
            } else {
                blocks[2] - 1 - block[2]
            };
            row * blocks[2] + z
        };
        let mut centres: Vec<(usize, [u32; 3], usize)> = (0..count)
            .filter(|&center| view.center_eligible[center])
            .map(|center| {
                let point = view.positions[center];
                (rank(self.cells.home(point)), site_key(point), center)
            })
            .collect();
        centres.sort_unstable();
        let slack = 0.3 * self.cells.cell * BLOCK_CELLS as f64;
        let mut reach = self.cells.typical + slack;
        let mut guess = (self.cells.typical * self.cells.typical) as f32;
        let mut start = 0;
        while start < centres.len() {
            let block = centres[start].0;
            let end = start + centres[start..].partition_point(|centre| centre.0 == block);
            let home = self.cells.home(view.positions[centres[start].2]);
            let low: [usize; 3] =
                std::array::from_fn(|axis| home[axis] / BLOCK_CELLS * BLOCK_CELLS);
            let high: [usize; 3] = std::array::from_fn(|axis| {
                (low[axis] + BLOCK_CELLS - 1).min(self.cells.dims[axis] - 1)
            });
            // Whether every crop in the block is rejected: for any centre when it keeps no
            // protein chain, and also for capped ligand centres when it cannot be over a limit.
            let (mut skip, mut skip_ligand) = (false, false);
            for _ in 0..SKIP_TRIES {
                let (protein, ligands) = self.survey(low, high, reach);
                if protein >= self.protein_floor {
                    break;
                }
                if self.inner(low, high) >= size {
                    (skip, skip_ligand) = (true, view.capped() && ligands <= ligand_bound);
                    break;
                }
                reach *= 1.5;
            }
            let skipped = |center: usize| {
                if view.ligand[center] {
                    skip_ligand
                } else {
                    skip
                }
            };
            if (start..end).all(|row| skipped(centres[row].2)) {
                start = end;
                continue;
            }
            self.fill();
            let mut farthest = 0_f32;
            // Whether the outcome of the previous row holds for every centre at its site.
            let mut shared = false;
            for row in start..end {
                let (_, site, center) = centres[row];
                if shared && centres[row - 1].1 == site {
                    outcomes[center] = outcomes[centres[row - 1].2];
                    continue;
                }
                if skipped(center) {
                    // A skipped ligand centre has over == false; a skipped other centre
                    // leaves it unknown.
                    shared = skip_ligand;
                    continue;
                }
                shared = true;
                let point = view.positions[center];
                let found = loop {
                    match self.search(point, (guess * GUESS_SHARE).to_bits()) {
                        Search::Miss => {
                            reach *= 1.2;
                            self.survey(low, high, reach);
                            self.fill();
                        }
                        Search::Over(estimate) => {
                            guess = estimate;
                            break None;
                        }
                        Search::Found(boundary) => break Some(boundary),
                    }
                };
                let Some(boundary) = found else {
                    outcomes[center].over = true;
                    farthest = farthest.max(guess);
                    continue;
                };
                let (over, accepted) = self.update(boundary);
                let radius = f32::from_bits((boundary >> 32) as u32);
                let bits = if accepted {
                    let at = store.len() as u32;
                    store.extend_from_slice(&self.tallies.words);
                    at
                } else {
                    u32::MAX
                };
                outcomes[center] = Outcome {
                    over,
                    bits,
                    radius,
                    tie: boundary as u32,
                    ligands: (self.tallies.ligand_count, self.tallies.ligand_chains),
                };
                guess = radius;
                farthest = farthest.max(radius);
            }
            reach = f64::from(farthest).sqrt() + slack;
            start = end;
        }
        (outcomes, store)
    }
}

fn compare(left: &(f32, usize), right: &(f32, usize)) -> Ordering {
    left.0
        .total_cmp(&right.0)
        .then_with(|| left.1.cmp(&right.1))
}

#[cfg(test)]
pub(crate) fn plan(
    positions: &[u8],
    chains: &[u8],
    ligand: &[u8],
    usable: &[u8],
    size: usize,
    center_eligible: Option<&[u8]>,
) -> Result<CropPlans, String> {
    plan_with_chain_policy(
        positions,
        chains,
        ligand,
        usable,
        size,
        center_eligible,
        None,
        None,
        Some(LIGAND_LIMITS),
    )
}

/// Plan every crop under the ligand limits, or under none with `ligand_limits` None.
///
/// Uncapped plans keep each crop that retains a protein chain and report the ligand counts of
/// every eligible ligand centre's crop. Rejecting the plans over chosen limits then gives
/// those limits' plans, and those counts give their retry weight.
#[allow(clippy::too_many_arguments)]
pub(crate) fn plan_with_chain_policy(
    positions: &[u8],
    chains: &[u8],
    ligand: &[u8],
    usable: &[u8],
    size: usize,
    center_eligible: Option<&[u8]>,
    chain_minimum: Option<&[u8]>,
    chain_protein: Option<&[u8]>,
    ligand_limits: Option<(u32, u32)>,
) -> Result<CropPlans, String> {
    let mut view = View::read(positions, chains, ligand, usable, size, center_eligible)?;
    (view.token_limit, view.chain_limit) = match ligand_limits {
        Some((tokens, chains)) if tokens < u32::MAX && chains < u32::MAX => (tokens, chains),
        Some(_) => return Err("Ligand crop limits must be below 2^32 - 1".into()),
        None => (u32::MAX, u32::MAX),
    };
    if let (Some(minimum), Some(protein)) = (chain_minimum, chain_protein) {
        if minimum.len() != view.chain_sizes.len() * 8
            || protein.len() != view.chain_sizes.len()
            || protein.iter().any(|&value| value > 1)
        {
            return Err("Crop chain policy must cover every chain ID".into());
        }
        view.chain_minimum = minimum
            .chunks_exact(8)
            .map(|row| {
                usize::try_from(i64::from_ne_bytes(std::array::from_fn(|i| row[i])))
                    .map_err(|_| "Crop chain minimum must be positive")
            })
            .collect::<Result<Vec<_>, _>>()?;
        view.chain_protein = protein.iter().map(|&value| value != 0).collect();
        if view
            .chain_minimum
            .iter()
            .enumerate()
            .any(|(chain, &count)| {
                count == 0 || view.chain_protein[chain] && (count < 50 || view.chain_ligand[chain])
            })
        {
            return Err(
                "Protein crop chains require at least 50 residues and cannot be ligands".into(),
            );
        }
    } else if chain_minimum.is_some() || chain_protein.is_some() {
        return Err("Crop chain minima and protein flags must be supplied together".into());
    }
    let total = view.center_eligible.iter().filter(|&&value| value).count();
    let ligand_count = view
        .ligand
        .iter()
        .zip(&view.center_eligible)
        .filter(|&(ligand, eligible)| *ligand && *eligible)
        .count();
    let mut result = CropPlans::empty(view.chain_sizes.len().div_ceil(64).max(1));
    let mut retry = 0;
    // A view of ligands alone keeps no protein chain; uncapped plans still report its counts.
    if total == 0 || view.capped() && view.ligand.iter().all(|&value| value) {
        return Ok(result);
    }
    let cells = (view.positions.len() >= SWEEP_MIN_TOKENS.max(4 * view.size))
        .then(|| Cells::build(&view))
        .flatten();
    if let Some(cells) = cells {
        let (outcomes, store) = Sweep::new(&view, cells).run();
        let words = result.words;
        for (center, outcome) in outcomes.iter().enumerate() {
            if !view.center_eligible[center] {
                continue;
            }
            retry += usize::from(outcome.over && view.ligand[center]);
            if view.ligand[center] && !view.capped() {
                result.centre_ligand_tokens.push(outcome.ligands.0);
                result.centre_ligand_chains.push(outcome.ligands.1);
            }
            if outcome.bits != u32::MAX {
                let at = outcome.bits as usize;
                result.centers.push(center as i64);
                result.radii.push(outcome.radius);
                result.ties.push(i64::from(outcome.tie));
                result.bits.extend_from_slice(&store[at..at + words]);
                result.ligand_tokens.push(outcome.ligands.0);
                result.ligand_chains.push(outcome.ligands.1);
            }
        }
        return Ok(finish(result, &view, retry, total, ligand_count));
    }
    let mut nearest = Vec::new();
    let mut tally = Tally::new(view.chain_sizes.len());
    // Centres at equal coordinates see equal distances, so they share one search and filter.
    let mut repeats: HashMap<[u32; 3], usize> = HashMap::new();
    for center in 0..view.positions.len() {
        if view.center_eligible[center] {
            *repeats.entry(site_key(view.positions[center])).or_default() += 1;
        }
    }
    repeats.retain(|_, count| *count > 1);
    let mut shared: HashMap<[u32; 3], Crop> = HashMap::new();
    for center in 0..view.positions.len() {
        if !view.center_eligible[center] {
            continue;
        }
        let key = site_key(view.positions[center]);
        let ((bits, over, ligands), boundary) = if let Some(crop) = shared.get(&key) {
            crop.clone()
        } else {
            view.nearest_all(center, &mut nearest);
            let crop = (view.filter(&nearest, &mut tally), nearest.last().copied());
            if repeats.contains_key(&key) {
                shared.insert(key, crop.clone());
            }
            crop
        };
        retry += usize::from(over && view.ligand[center]);
        if view.ligand[center] && !view.capped() {
            result.centre_ligand_tokens.push(ligands.0);
            result.centre_ligand_chains.push(ligands.1);
        }
        if let Some(bits) = bits {
            let (radius, tie) = boundary.ok_or("An accepted spatial crop has no neighbors")?;
            result.centers.push(center as i64);
            result.radii.push(radius);
            result.ties.push(tie as i64);
            result.bits.extend(bits);
            result.ligand_tokens.push(ligands.0);
            result.ligand_chains.push(ligands.1);
        }
    }
    Ok(finish(result, &view, retry, total, ligand_count))
}

/// Weight the accepted centres and normalise the weights.
fn finish(
    mut result: CropPlans,
    view: &View,
    retry: usize,
    total: usize,
    ligand_count: usize,
) -> CropPlans {
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
    result
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
