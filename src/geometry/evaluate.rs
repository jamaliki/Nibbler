use super::math::{cross, dot, norm, plane_deltas, sub};
use super::{Evaluation, Finding, RestraintGraph};

impl RestraintGraph {
    /// Evaluate supplied coordinates, failing closed on undefined geometry.
    pub fn evaluate(&self, sites: &[[f64; 3]], cutoff: f64) -> Result<Evaluation, String> {
        if sites.len() != self.atom_ids.len() || !sites.iter().flatten().all(|x| x.is_finite()) {
            return Err("Coordinates must cover every atom and contain finite values".into());
        }
        if !cutoff.is_finite() || cutoff < 0.0 {
            return Err("Sigma cutoff must be finite and nonnegative".into());
        }
        let mut result = Evaluation {
            near_cutoff: false,
            counts: vec![
                ("bonds", self.bonds.len()),
                ("angles", self.angles.len()),
                ("chiralities", self.chiralities.len()),
                ("planarities", self.planes.len()),
                ("hydrogens", self.hydrogens.iter().filter(|&&h| h).count()),
            ],
            findings: Vec::new(),
            hydrogen_failures: Vec::new(),
        };
        for (proxy_index, &(i, j, target, weight, slack, report)) in self.bonds.iter().enumerate() {
            let model = norm(sub(sites[i], sites[j]));
            let difference = target - model;
            let delta_slack = difference.signum() * (difference.abs() - slack).max(0.0);
            record(
                &mut result,
                cutoff,
                report,
                "bonds",
                proxy_index,
                &[i, j],
                model,
                target,
                difference,
                weight,
            )?;
            if report
                && let Some(finding) = result.findings.last_mut()
                && finding.check == "bonds"
                && finding.proxy_index == proxy_index
            {
                finding.residual = weight * delta_slack * delta_slack;
            }
        }
        for (proxy_index, &(i, j, k, target, weight)) in self.angles.iter().enumerate() {
            let (a, b) = (sub(sites[i], sites[j]), sub(sites[k], sites[j]));
            let denominator = norm(a) * norm(b);
            if denominator == 0.0 || !denominator.is_finite() {
                return Err("Angle coordinates contain a zero or overflowed arm".into());
            }
            let model = (dot(a, b) / denominator)
                .clamp(-1.0, 1.0)
                .acos()
                .to_degrees();
            record(
                &mut result,
                cutoff,
                true,
                "angles",
                proxy_index,
                &[i, j, k],
                model,
                target,
                target - model,
                weight,
            )?;
        }
        for (proxy_index, &(i, j, k, l, target, weight, both_signs, pseudo)) in
            self.chiralities.iter().enumerate()
        {
            let model = dot(
                sub(sites[j], sites[i]),
                cross(sub(sites[k], sites[i]), sub(sites[l], sites[i])),
            );
            let delta = target - if both_signs { model.abs() } else { model };
            record(
                &mut result,
                cutoff,
                !pseudo,
                "chiralities",
                proxy_index,
                &[i, j, k, l],
                model,
                target,
                delta,
                weight,
            )?;
        }
        for (proxy_index, (indices, weights)) in self.planes.iter().enumerate() {
            let deltas = plane_deltas(sites, indices, weights)?;
            let score = deltas
                .iter()
                .zip(weights)
                .map(|(d, w)| d.abs() * w.sqrt())
                .fold(0.0_f64, f64::max);
            let residual = deltas
                .iter()
                .zip(weights)
                .map(|(d, w)| d * d * w)
                .sum::<f64>();
            let model = (deltas.iter().map(|d| d * d).sum::<f64>() / deltas.len() as f64).sqrt();
            if !score.is_finite() || !residual.is_finite() || !model.is_finite() {
                return Err("Plane evaluation overflowed".into());
            }
            result.near_cutoff |= (score - cutoff).abs() <= 1e-8 * cutoff.max(1.0);
            if score > cutoff {
                result.findings.push(Finding {
                    check: "planarities",
                    proxy_index,
                    indices: indices.clone(),
                    model,
                    target: 0.0,
                    delta: deltas.iter().map(|d| d.abs()).fold(0.0, f64::max),
                    sigma: None,
                    score,
                    residual,
                });
            }
        }
        for (i, &hydrogen) in self.hydrogens.iter().enumerate() {
            if hydrogen && (self.parents[i].len() != 1 || self.hydrogens[self.parents[i][0]]) {
                result.hydrogen_failures.push((i, self.parents[i].clone()));
            }
        }
        Ok(result)
    }
}

#[allow(clippy::too_many_arguments)]
fn record(
    result: &mut Evaluation,
    cutoff: f64,
    report: bool,
    check: &'static str,
    proxy_index: usize,
    indices: &[usize],
    model: f64,
    target: f64,
    delta: f64,
    weight: f64,
) -> Result<(), String> {
    let score = delta.abs() * weight.sqrt();
    let residual = weight * delta * delta;
    if !model.is_finite() || !score.is_finite() || !residual.is_finite() {
        return Err("Restraint evaluation overflowed".into());
    }
    result.near_cutoff |= report && (score - cutoff).abs() <= 1e-8 * cutoff.max(1.0);
    if report && score > cutoff {
        result.findings.push(Finding {
            check,
            proxy_index,
            indices: indices.to_vec(),
            model,
            target,
            delta,
            sigma: Some(1.0 / weight.sqrt()),
            score,
            residual,
        });
    }
    Ok(())
}
