//! Evaluate explicit chemical restraints without rebuilding their topology.
//!
//! The caller owns chemical interpretation and symmetry expansion. This graph
//! does not infer bonds, chemical states, or Probe contact rules.

mod evaluate;
mod math;

use std::collections::HashSet;

/// A bond's atom indices, ideal distance, weight, slack and reporting flag.
pub type Bond = (usize, usize, f64, f64, f64, bool);
/// An angle's atom indices, ideal angle in degrees and weight.
pub type Angle = (usize, usize, usize, f64, f64);
/// A chiral centre's indices, ideal volume, weight and sign/pseudochiral flags.
pub type Chirality = (usize, usize, usize, usize, f64, f64, bool, bool);
/// A plane's atom indices and corresponding weights.
pub type Plane = (Vec<usize>, Vec<f64>);

/// An immutable, validated sparse restraint graph.
#[derive(Clone)]
pub struct RestraintGraph {
    pub(crate) atom_ids: Vec<String>,
    pub(crate) hydrogens: Vec<bool>,
    pub(crate) bonds: Vec<Bond>,
    pub(crate) angles: Vec<Angle>,
    pub(crate) chiralities: Vec<Chirality>,
    pub(crate) planes: Vec<Plane>,
    pub(crate) parents: Vec<Vec<usize>>,
}

/// A restraint finding with the source proxy index retained.
#[derive(Debug)]
pub struct Finding {
    /// Check name.
    pub check: &'static str,
    /// Index in the corresponding source proxy array.
    pub proxy_index: usize,
    /// Atom indices, in source proxy order.
    pub indices: Vec<usize>,
    /// Observed distance, angle, volume, or plane RMS distance.
    pub model: f64,
    /// Target distance, angle, volume, or zero for planes.
    pub target: f64,
    /// Target minus observed value; planes use maximum absolute distance.
    pub delta: f64,
    /// Restraint standard deviation; planes have individual atom weights.
    pub sigma: Option<f64>,
    /// Absolute standardized deviation, or maximum atom deviation for planes.
    pub score: f64,
    /// Weighted sum of squared deviations.
    pub residual: f64,
}

/// Results from all supplied restraints and explicit hydrogen parents.
#[derive(Debug)]
pub struct Evaluation {
    /// A reportable score lies within floating-point comparison tolerance.
    /// Callers that require reference decisions should use reference evaluation.
    pub near_cutoff: bool,
    /// Check names and evaluated proxy counts.
    pub counts: Vec<(&'static str, usize)>,
    /// Geometric outliers, excluding caller-marked pseudochiral centres.
    pub findings: Vec<Finding>,
    /// Hydrogen atom indices and their bonded neighbours when invalid.
    pub hydrogen_failures: Vec<(usize, Vec<usize>)>,
}

fn validate_indices(indices: &[usize], count: usize, minimum: usize) -> Result<(), String> {
    if indices.len() < minimum || indices.iter().any(|&index| index >= count) {
        return Err("Restraint atom indices are out of range or missing".into());
    }
    if indices.iter().copied().collect::<HashSet<_>>().len() != indices.len() {
        return Err("Restraint atom indices must be distinct".into());
    }
    Ok(())
}

fn validate_weight(weight: f64) -> Result<(), String> {
    if !weight.is_finite() || weight <= 0.0 {
        return Err("Restraint weights must be finite and positive".into());
    }
    Ok(())
}

impl RestraintGraph {
    /// Validate and cache explicit restraints and hydrogen connectivity.
    pub fn new(
        atom_ids: Vec<String>,
        hydrogens: Vec<bool>,
        bonds: Vec<Bond>,
        angles: Vec<Angle>,
        chiralities: Vec<Chirality>,
        planes: Vec<Plane>,
    ) -> Result<Self, String> {
        let count = atom_ids.len();
        if count == 0 || hydrogens.len() != count || atom_ids.iter().any(String::is_empty) {
            return Err("Atom IDs and hydrogen flags must cover every atom".into());
        }
        if atom_ids.iter().collect::<HashSet<_>>().len() != count {
            return Err("Atom IDs must be unique".into());
        }
        let mut parents = vec![Vec::new(); count];
        for &(i, j, target, weight, slack, _) in &bonds {
            validate_indices(&[i, j], count, 2)?;
            validate_weight(weight)?;
            if !target.is_finite() || target <= 0.0 || !slack.is_finite() || slack < 0.0 {
                return Err("Bond targets and slack must be finite and valid".into());
            }
            if hydrogens[i] {
                parents[i].push(j);
            }
            if hydrogens[j] {
                parents[j].push(i);
            }
        }
        for &(i, j, k, target, weight) in &angles {
            validate_indices(&[i, j, k], count, 3)?;
            validate_weight(weight)?;
            if !target.is_finite() || !(0.0..=180.0).contains(&target) {
                return Err("Angle targets must lie between zero and 180 degrees".into());
            }
        }
        for &(i, j, k, l, target, weight, _, _) in &chiralities {
            validate_indices(&[i, j, k, l], count, 4)?;
            validate_weight(weight)?;
            if !target.is_finite() {
                return Err("Chiral targets must be finite".into());
            }
        }
        for (indices, weights) in &planes {
            validate_indices(indices, count, 3)?;
            if weights.len() != indices.len() {
                return Err("Plane weights must cover every plane atom".into());
            }
            for &weight in weights {
                validate_weight(weight)?;
            }
        }
        Ok(Self {
            atom_ids,
            hydrogens,
            bonds,
            angles,
            chiralities,
            planes,
            parents,
        })
    }
}
