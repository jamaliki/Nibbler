// The cyclic Jacobi implementation adapts
// cctbx_project/scitbx/matrix/eigensystem.h (CCTBX 2026.9).
// Copyright (c) 2006 - 2026, The Regents of the University of California,
// through Lawrence Berkeley National Laboratory. See
// licenses/CCTBX-LICENSE.txt and licenses/CCTBX-COPYRIGHT.txt.
// The upstream implementation cites BTL3's matrix algorithms, Ian Tickle
// (1992), David Moss (1993), and Mark Williams (1998).

pub(super) fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

pub(super) fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

pub(super) fn norm(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
}

pub(super) fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

// Match scitbx's cyclic Jacobi thresholds and rotations. Its stopping tolerance
// affects outlier decisions for poorly conditioned weighted planes.
fn smallest_eigenvector(mut matrix: [[f64; 3]; 3]) -> Result<[f64; 3], String> {
    let mut vectors = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    let norm = (2.0
        * (matrix[1][0] * matrix[1][0]
            + matrix[2][0] * matrix[2][0]
            + matrix[2][1] * matrix[2][1]))
        .sqrt();
    let limit = 1e-10 * norm / 3.0;
    let mut threshold = norm;
    let mut iterations = 0;
    while threshold > limit {
        threshold /= 3.0;
        loop {
            let mut changed = false;
            for (p, q) in [(0, 1), (0, 2), (1, 2)] {
                if matrix[p][q] * matrix[p][q] <= threshold * threshold {
                    continue;
                }
                changed = true;
                let (pp, qq, pq) = (matrix[p][p], matrix[q][q], matrix[p][q]);
                let x = 0.5 * (pp - qq);
                let mut y = -pq / (pq * pq + x * x).sqrt();
                if x < 0.0 {
                    y = -y;
                }
                let s = y / (2.0 * (1.0 + (1.0 - y * y).sqrt())).sqrt();
                let c = (1.0 - s * s).sqrt();
                for k in 0..3 {
                    if k != p && k != q {
                        let (kp, kq) = (matrix[k][p], matrix[k][q]);
                        matrix[k][p] = kp * c - kq * s;
                        matrix[p][k] = matrix[k][p];
                        matrix[k][q] = kp * s + kq * c;
                        matrix[q][k] = matrix[k][q];
                    }
                    let (vp, vq) = (vectors[k][p], vectors[k][q]);
                    vectors[k][p] = vp * c - vq * s;
                    vectors[k][q] = vp * s + vq * c;
                }
                matrix[p][p] = pp * c * c + qq * s * s - 2.0 * pq * s * c;
                matrix[q][q] = pp * s * s + qq * c * c + 2.0 * pq * s * c;
                matrix[p][q] = (pp - qq) * s * c + pq * (c * c - s * s);
                matrix[q][p] = matrix[p][q];
            }
            iterations += 1;
            if iterations > 1000 || !matrix.iter().flatten().all(|v| v.is_finite()) {
                return Err("Plane eigensolver did not converge".into());
            }
            if !changed {
                break;
            }
        }
    }
    let mut order = [0, 1, 2];
    order.sort_by(|&i, &j| matrix[i][i].total_cmp(&matrix[j][j]));
    let scale = matrix.iter().flatten().fold(0.0_f64, |a, b| a.max(b.abs()));
    if matrix[order[1]][order[1]] <= scale * 32.0 * f64::EPSILON
        || matrix[order[1]][order[1]] - matrix[order[0]][order[0]] <= scale * 1e-10
    {
        return Err("Plane coordinates do not define a unique plane".into());
    }
    let i = order[0];
    Ok([vectors[0][i], vectors[1][i], vectors[2][i]])
}

pub(super) fn plane_deltas(
    coordinates: &[[f64; 3]],
    indices: &[usize],
    weights: &[f64],
) -> Result<Vec<f64>, String> {
    let sum: f64 = weights.iter().sum();
    let mut centre = [0.0; 3];
    for (&i, &weight) in indices.iter().zip(weights) {
        for (axis, value) in centre.iter_mut().enumerate() {
            *value += coordinates[i][axis] * weight;
        }
    }
    for value in &mut centre {
        *value /= sum;
    }
    let mut tensor = [[0.0; 3]; 3];
    for (&i, &weight) in indices.iter().zip(weights) {
        let x = sub(coordinates[i], centre);
        for a in 0..3 {
            for b in 0..3 {
                tensor[a][b] += weight * x[a] * x[b];
            }
        }
    }
    if !tensor.iter().flatten().all(|value| value.is_finite()) {
        return Err("Plane tensor overflowed".into());
    }
    let normal = smallest_eigenvector(tensor)?;
    Ok(indices
        .iter()
        .map(|&i| dot(normal, sub(coordinates[i], centre)))
        .collect())
}
