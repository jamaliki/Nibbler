//! Shrake-Rupley areas for selected atoms, using a sparse neighbour grid.

use std::collections::{HashMap, HashSet};
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
    for band in points.chunks_mut(16) {
        band.sort_by(|a, b| a[1].atan2(a[0]).total_cmp(&b[1].atan2(b[0])));
    }
    points
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
    let unique: HashSet<_> = selected.iter().copied().collect();
    if unique.len() != selected.len() || unique.iter().any(|&i| i >= radii.len()) {
        return Err("atom_indices must contain distinct, in-range atom indices".into());
    }
    Ok(())
}

fn adjacent(c: [i64; 3]) -> impl Iterator<Item = [i64; 3]> {
    (-1..=1).flat_map(move |x| {
        (-1..=1).flat_map(move |y| (-1..=1).map(move |z| [c[0] + x, c[1] + y, c[2] + z]))
    })
}

const EMPTY: usize = usize::MAX;

enum Cells {
    Dense(Vec<usize>, [usize; 3]),
    Sparse(HashMap<[i64; 3], usize>),
}

struct Grid<'a> {
    coordinates: &'a [[f64; 3]],
    expanded: Vec<f64>,
    width: f64,
    lower: [i64; 3],
    upper: [i64; 3],
    cells: Cells,
    next: Vec<usize>,
}

impl<'a> Grid<'a> {
    fn cell(&self, xyz: [f64; 3]) -> [i64; 3] {
        xyz.map(|x| (x / self.width).floor() as i64)
    }

    fn slot(&self, c: [i64; 3], shape: [usize; 3]) -> usize {
        let d: [usize; 3] = std::array::from_fn(|k| (c[k] - self.lower[k]) as usize);
        (d[0] * shape[1] + d[1]) * shape[2] + d[2]
    }

    fn head(&self, c: [i64; 3]) -> usize {
        if (0..3).any(|k| c[k] < self.lower[k] || c[k] > self.upper[k]) {
            return EMPTY;
        }
        match &self.cells {
            Cells::Dense(heads, shape) => heads[self.slot(c, *shape)],
            Cells::Sparse(heads) => heads.get(&c).copied().unwrap_or(EMPTY),
        }
    }

    fn insert(&mut self, c: [i64; 3], i: usize) {
        if (0..3).any(|k| c[k] < self.lower[k] || c[k] > self.upper[k]) {
            return;
        }
        let slot = match &self.cells {
            Cells::Dense(_, shape) => self.slot(c, *shape),
            Cells::Sparse(_) => 0,
        };
        self.next[i] = match &mut self.cells {
            Cells::Dense(heads, _) => std::mem::replace(&mut heads[slot], i),
            Cells::Sparse(heads) => heads.insert(c, i).unwrap_or(EMPTY),
        };
    }

    fn new(coordinates: &'a [[f64; 3]], radii: &[f64], selected: &[usize], probe: f64) -> Self {
        let expanded: Vec<_> = radii.iter().map(|r| r + probe).collect();
        // The lower bound keeps cell indices within i64 for very small radii.
        let width = (2.0 * expanded.iter().copied().fold(0.0_f64, f64::max)).max(1e-6);
        let cell = |i: usize, k: usize| (coordinates[i][k] / width).floor() as i64;
        let lower: [i64; 3] =
            std::array::from_fn(|k| selected.iter().map(|&i| cell(i, k)).min().unwrap_or(0) - 1);
        let upper: [i64; 3] =
            std::array::from_fn(|k| selected.iter().map(|&i| cell(i, k)).max().unwrap_or(0) + 1);
        let shape: [usize; 3] = std::array::from_fn(|k| (upper[k] - lower[k] + 1) as usize);
        let count = shape.into_iter().try_fold(1_usize, usize::checked_mul);
        // Sparse coordinates never force a large dense allocation.
        let cells = match count {
            Some(n) if n <= 1_000_000 => Cells::Dense(vec![EMPTY; n], shape),
            _ => Cells::Sparse(HashMap::new()),
        };
        let mut grid = Self {
            coordinates,
            expanded,
            width,
            lower,
            upper,
            cells,
            next: vec![EMPTY; coordinates.len()],
        };
        for (i, &xyz) in coordinates.iter().enumerate() {
            grid.insert(grid.cell(xyz), i);
        }
        grid
    }

    fn area(&self, i: usize, points: &[[f64; 3]]) -> f64 {
        let r = self.expanded[i];
        let mut caps = Vec::new();
        for c in adjacent(self.cell(self.coordinates[i])) {
            let mut j = self.head(c);
            while j != EMPTY {
                let current = j;
                j = self.next[j];
                if i == current {
                    continue;
                }
                let delta = std::array::from_fn::<_, 3, _>(|k| {
                    self.coordinates[current][k] - self.coordinates[i][k]
                });
                let d2: f64 = delta.iter().map(|x| x * x).sum();
                let other = self.expanded[current];
                if d2 >= (r + other).powi(2) {
                    continue;
                }
                if other >= r && d2 <= (other - r).powi(2) {
                    return 0.0;
                }
                if r >= other && d2 <= (r - other).powi(2) {
                    continue;
                }
                let distance = d2.sqrt();
                caps.push((
                    delta.map(|x| x / distance),
                    (d2 + r * r - other * other) / (2.0 * r * distance),
                ));
            }
        }
        caps.sort_unstable_by(|a, b| a.1.total_cmp(&b.1));
        let exposed = exposed_points(points, &caps);
        4.0 * PI * r * r * exposed as f64 / points.len() as f64
    }
}

fn exposed_points(points: &[[f64; 3]], caps: &[([f64; 3], f64)]) -> usize {
    if caps.is_empty() {
        return points.len();
    }
    // Adjacent Fibonacci points often share an occluder.
    let mut last = 0;
    let mut exposed = 0;
    for p in points {
        let blocked = |index: usize| {
            let (d, threshold) = caps[index];
            p[0] * d[0] + p[1] * d[1] + p[2] * d[2] > threshold
        };
        if blocked(last) {
            continue;
        }
        if let Some(index) = (0..caps.len()).find(|&index| blocked(index)) {
            last = index;
        } else {
            exposed += 1;
        }
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
    Ok(atom_indices.iter().map(|&i| grid.area(i, points)).collect())
}
