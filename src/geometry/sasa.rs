//! Shrake-Rupley areas for selected atoms, using a cell-sorted neighbour grid.

use std::collections::HashMap;
use std::f64::consts::PI;
use std::sync::OnceLock;

/// Default number of Fibonacci sphere points per selected atom.
pub const SASA_POINTS: usize = 256;

fn sphere(count: usize) -> Vec<[f64; 3]> {
    let angle = PI * (3.0 - 5.0_f64.sqrt());
    let mut points: Vec<_> = (0..count)
        .map(|i| {
            let z = 1.0 - 2.0 * (i as f64 + 0.5) / count as f64;
            let r = (1.0 - z * z).sqrt();
            let phi = angle * i as f64;
            [r * phi.cos(), r * phi.sin(), z]
        })
        .collect();
    if count > 4096 || points.is_empty() {
        for band in points.chunks_mut(16) {
            band.sort_by(|a, b| a[1].atan2(a[0]).total_cmp(&b[1].atan2(b[0])));
        }
        return points;
    }
    // A nearest-neighbour tour keeps consecutive points under the same occluder.
    let mut tour = Vec::with_capacity(count);
    let mut current = points.swap_remove(0);
    tour.push(current);
    while !points.is_empty() {
        let distance = |p: &[f64; 3]| (0..3).map(|k| (p[k] - current[k]).powi(2)).sum::<f64>();
        let next = (0..points.len())
            .min_by(|&a, &b| distance(&points[a]).total_cmp(&distance(&points[b])))
            .unwrap_or(0);
        current = points.swap_remove(next);
        tour.push(current);
    }
    tour
}

pub(crate) fn validate_sasa_options(probe_radius: f64, point_count: usize) -> Result<(), String> {
    if !probe_radius.is_finite() || !(0.0..=100.0).contains(&probe_radius) {
        return Err("probe_radius must be finite and between 0 and 100 A".into());
    }
    if !(1..=1_000_000).contains(&point_count) {
        return Err("point_count must be between 1 and 1000000".into());
    }
    Ok(())
}

fn validate_atoms(
    coordinates: &[[f64; 3]],
    radii: &[f64],
    selected: &[usize],
) -> Result<(), String> {
    if coordinates.len() != radii.len()
        || coordinates
            .iter()
            .flatten()
            .any(|x| !x.is_finite() || x.abs() > 1e9)
        || radii
            .iter()
            .any(|r| !r.is_finite() || *r <= 0.0 || *r > 100.0)
    {
        return Err(
            "Coordinates and radii must have equal lengths and finite, valid values".into(),
        );
    }
    let mut seen = vec![false; radii.len()];
    for &i in selected {
        if i >= radii.len() || std::mem::replace(&mut seen[i], true) {
            return Err("atom_indices must contain distinct, in-range atom indices".into());
        }
    }
    Ok(())
}

const SPARSE: u32 = u32::MAX;

/// Atoms sorted by cell, with coordinates and expanded radii stored by column
/// so the neighbour scan reads contiguous memory.
struct Grid {
    width: f64,
    lower: [i64; 3],
    upper: [i64; 3],
    shape: [usize; 3],
    /// Dense cell start offsets; empty when the grid is sparse.
    starts: Vec<u32>,
    sparse: HashMap<[i64; 3], (u32, u32)>,
    x: Vec<f64>,
    y: Vec<f64>,
    z: Vec<f64>,
    expanded: Vec<f64>,
    /// Sorted position of each input atom, or `SPARSE` outside the grid.
    position: Vec<u32>,
}

/// Neighbour caps of one atom, widest first, as unit directions and cosine
/// thresholds stored by column and padded to whole blocks.
#[derive(Default)]
struct Caps {
    near: Vec<u32>,
    /// Directions, then squared distances and radii, later thresholds.
    columns: [Vec<f64>; 5],
    keys: Vec<u64>,
    count: usize,
    blocks: Vec<Block>,
}

/// Four caps by column, so one block test compiles to vector arithmetic.
#[derive(Clone, Copy)]
struct Block {
    x: [f64; LANES],
    y: [f64; LANES],
    z: [f64; LANES],
    threshold: [f64; LANES],
}

impl Block {
    /// Never blocks a point.
    const EMPTY: Self = Self { x: [0.0; LANES], y: [0.0; LANES], z: [0.0; LANES], threshold: [f64::INFINITY; LANES] };

    /// Bit `lane` is set when that cap blocks point `p`.
    #[inline(always)]
    fn hits(&self, p: [f64; 3]) -> u32 {
        let mut mask = 0;
        for lane in 0..LANES {
            mask |= u32::from(p[0] * self.x[lane] + p[1] * self.y[lane] + p[2] * self.z[lane] > self.threshold[lane]) << lane;
        }
        mask
    }
}

impl Grid {
    fn cell(&self, xyz: [f64; 3]) -> [i64; 3] {
        xyz.map(|x| (x / self.width).floor() as i64)
    }

    fn inside(&self, c: [i64; 3]) -> bool {
        (0..3).all(|k| c[k] >= self.lower[k] && c[k] <= self.upper[k])
    }

    fn slot(&self, c: [i64; 3]) -> usize {
        let d: [usize; 3] = std::array::from_fn(|k| (c[k] - self.lower[k]) as usize);
        (d[0] * self.shape[1] + d[1]) * self.shape[2] + d[2]
    }

    fn new(coordinates: &[[f64; 3]], radii: &[f64], selected: &[usize], probe: f64) -> Self {
        let expanded: Vec<_> = radii.iter().map(|r| r + probe).collect();
        // Half-width cells shrink the scanned volume. The margin keeps every
        // contact within two cells despite rounding in the cell division.
        let reach = 2.0 * expanded.iter().copied().fold(0.0_f64, f64::max);
        // The lower bound keeps cell indices within i64 for very small radii.
        let width = (0.5 * reach * (1.0 + 1e-9)).max(1e-6);
        let cell = |i: usize, k: usize| (coordinates[i][k] / width).floor() as i64;
        let lower: [i64; 3] =
            std::array::from_fn(|k| selected.iter().map(|&i| cell(i, k)).min().unwrap_or(0) - 2);
        let upper: [i64; 3] =
            std::array::from_fn(|k| selected.iter().map(|&i| cell(i, k)).max().unwrap_or(0) + 2);
        let shape: [usize; 3] = std::array::from_fn(|k| (upper[k] - lower[k] + 1) as usize);
        let count = shape.into_iter().try_fold(1_usize, usize::checked_mul);
        let mut grid = Self {
            width,
            lower,
            upper,
            shape,
            starts: Vec::new(),
            sparse: HashMap::new(),
            x: Vec::new(),
            y: Vec::new(),
            z: Vec::new(),
            expanded: Vec::new(),
            position: vec![SPARSE; coordinates.len()],
        };
        let cells: Vec<[i64; 3]> = coordinates.iter().map(|&xyz| grid.cell(xyz)).collect();
        // Sparse coordinates never force a large dense allocation.
        let order = match count.filter(|&n| n <= 4_000_000) {
            Some(n) => {
                let mut starts = vec![0_u32; n + 1];
                let mut keys = Vec::with_capacity(coordinates.len());
                for (i, &c) in cells.iter().enumerate() {
                    if grid.inside(c) {
                        let slot = grid.slot(c);
                        starts[slot + 1] += 1;
                        keys.push((slot, i));
                    }
                }
                for slot in 0..n {
                    starts[slot + 1] += starts[slot];
                }
                let mut fill = starts.clone();
                let mut order = vec![0_usize; keys.len()];
                for &(slot, i) in &keys {
                    order[fill[slot] as usize] = i;
                    fill[slot] += 1;
                }
                grid.starts = starts;
                order
            }
            None => {
                let mut order: Vec<usize> =
                    (0..coordinates.len()).filter(|&i| grid.inside(cells[i])).collect();
                order.sort_by_key(|&i| cells[i]);
                let mut begin = 0;
                for (index, &i) in order.iter().enumerate() {
                    if order.get(index + 1).is_none_or(|&next| cells[next] != cells[i]) {
                        grid.sparse.insert(cells[i], (begin as u32, index as u32 + 1));
                        begin = index + 1;
                    }
                }
                order
            }
        };
        for (position, &i) in order.iter().enumerate() {
            grid.position[i] = position as u32;
        }
        grid.x = order.iter().map(|&i| coordinates[i][0]).collect();
        grid.y = order.iter().map(|&i| coordinates[i][1]).collect();
        grid.z = order.iter().map(|&i| coordinates[i][2]).collect();
        grid.expanded = order.iter().map(|&i| expanded[i]).collect();
        grid
    }

    /// Visit the sorted ranges of the 5 x 5 x 5 cells around one cell.
    fn ranges(&self, c: [i64; 3], mut visit: impl FnMut(usize, usize) -> bool) {
        for dx in -2..=2 {
            for dy in -2..=2 {
                if self.sparse.is_empty() {
                    // Cells along z share one sorted run.
                    let first = self.slot([c[0] + dx, c[1] + dy, c[2] - 2]);
                    if !visit(self.starts[first] as usize, self.starts[first + 5] as usize) {
                        return;
                    }
                } else {
                    for dz in -2..=2 {
                        if let Some(&(begin, end)) = self.sparse.get(&[c[0] + dx, c[1] + dy, c[2] + dz])
                            && !visit(begin as usize, end as usize)
                        {
                            return;
                        }
                    }
                }
            }
        }
    }

    /// Collect the caps that occlude selected atom `i`, or return false when buried.
    fn caps(&self, i: usize, caps: &mut Caps) -> bool {
        let own = self.position[i] as usize;
        let centre = [self.x[own], self.y[own], self.z[own]];
        let r = self.expanded[own];
        // Keep contacts without branching, then classify the few that remain.
        let mut count = 0;
        self.ranges(self.cell(centre), |begin, end| {
            if caps.near.len() < count + end - begin {
                caps.near.resize(count + end - begin, 0);
            }
            for j in begin..end {
                let delta = [self.x[j] - centre[0], self.y[j] - centre[1], self.z[j] - centre[2]];
                let d2 = delta[0] * delta[0] + delta[1] * delta[1] + delta[2] * delta[2];
                caps.near[count] = j as u32;
                count += usize::from(d2 < (r + self.expanded[j]).powi(2));
            }
            true
        });
        caps.near.truncate(count);
        let columns = &mut caps.columns;
        for column in columns.iter_mut() {
            column.clear();
        }
        for &j in &caps.near {
            let j = j as usize;
            if j == own {
                continue;
            }
            let delta = [self.x[j] - centre[0], self.y[j] - centre[1], self.z[j] - centre[2]];
            let d2 = delta[0] * delta[0] + delta[1] * delta[1] + delta[2] * delta[2];
            let other = self.expanded[j];
            if other >= r && d2 <= (other - r).powi(2) {
                return false;
            }
            if r >= other && d2 <= (r - other).powi(2) {
                continue;
            }
            for (column, value) in columns.iter_mut().zip([delta[0], delta[1], delta[2], d2, other]) {
                column.push(value);
            }
        }
        // Unit directions and cosine thresholds, computed by column.
        let [x, y, z, d2, other] = columns;
        let count = d2.len();
        for k in 0..count {
            let distance = d2[k].sqrt();
            let threshold = (d2[k] + r * r - other[k] * other[k]) / (2.0 * r * distance);
            x[k] /= distance;
            y[k] /= distance;
            z[k] /= distance;
            d2[k] = threshold;
        }
        let threshold = &*d2;
        // Wide caps first: they block the most points. The low key bits hold the
        // cap index, which only perturbs the order of near-equal thresholds.
        let packed = count <= 0x10000;
        caps.keys.clear();
        caps.keys.extend(threshold.iter().enumerate().map(|(k, &value)| {
            if packed { order_key(value) & !0xFFFF | k as u64 } else { k as u64 }
        }));
        if packed {
            caps.keys.sort_unstable();
        }
        caps.count = count;
        caps.blocks.clear();
        caps.blocks.resize(count.div_ceil(LANES), Block::EMPTY);
        for (index, &key) in caps.keys.iter().enumerate() {
            let k = if packed { key & 0xFFFF } else { key } as usize;
            let block = &mut caps.blocks[index / LANES];
            let lane = index % LANES;
            block.x[lane] = x[k];
            block.y[lane] = y[k];
            block.z[lane] = z[k];
            block.threshold[lane] = threshold[k];
        }
        true
    }

    fn area(&self, i: usize, points: &[[f64; 3]], caps: &mut Caps) -> f64 {
        let r = self.expanded[self.position[i] as usize];
        if !self.caps(i, caps) {
            return 0.0;
        }
        let exposed = exposed_points(points, caps);
        4.0 * PI * r * r * exposed as f64 / points.len() as f64
    }
}

/// Map a float to an integer with the same order as `f64::total_cmp`.
fn order_key(value: f64) -> u64 {
    let bits = value.to_bits();
    if bits >> 63 == 1 { !bits } else { bits | 1 << 63 }
}

const LANES: usize = 4;

fn exposed_points(points: &[[f64; 3]], caps: &Caps) -> usize {
    if caps.count == 0 {
        return points.len();
    }
    let blocks = &caps.blocks;
    // Consecutive tour points often share an occluder.
    let mut last = 0;
    let mut exposed = 0;
    'points: for &p in points {
        if blocks[last].hits(p) != 0 {
            continue;
        }
        for (index, block) in blocks.iter().enumerate() {
            if block.hits(p) != 0 {
                last = index;
                continue 'points;
            }
        }
        exposed += 1;
    }
    exposed
}

/// Estimate selected atoms' solvent-accessible areas in square angstroms.
///
/// Coordinates and positive van der Waals radii describe all occluding atoms.
/// Results follow `atom_indices` order. Other atoms receive no surface samples.
/// Coincident spheres occlude one another, including equal-radius duplicates.
/// This function does not infer elements, select conformers, or expand symmetry.
///
/// # Errors
///
/// Rejects mismatched arrays, repeated or invalid indices, non-finite coordinates,
/// coordinates beyond +/-1e9 A, radii outside (0, 100] A, and invalid options.
pub fn sasa(
    coordinates: &[[f64; 3]],
    radii: &[f64],
    atom_indices: &[usize],
    probe_radius: f64,
    point_count: usize,
) -> Result<Vec<f64>, String> {
    validate_sasa_options(probe_radius, point_count)?;
    validate_atoms(coordinates, radii, atom_indices)?;
    if atom_indices.is_empty() {
        return Ok(Vec::new());
    }
    let grid = Grid::new(coordinates, radii, atom_indices, probe_radius);
    let mut caps = Caps::default();
    static DEFAULT_SPHERE: OnceLock<Vec<[f64; 3]>> = OnceLock::new();
    let custom;
    let points = if point_count == SASA_POINTS {
        DEFAULT_SPHERE
            .get_or_init(|| sphere(SASA_POINTS))
            .as_slice()
    } else {
        custom = sphere(point_count);
        &custom
    };
    Ok(atom_indices.iter().map(|&i| grid.area(i, points, &mut caps)).collect())
}
