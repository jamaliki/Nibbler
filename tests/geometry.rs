//! Tests for explicit chemical-restraint evaluation.

#![allow(clippy::expect_used)]

use _core::geometry::RestraintGraph;

fn graph(
    bonds: Vec<_core::geometry::Bond>,
    angles: Vec<_core::geometry::Angle>,
    chirals: Vec<_core::geometry::Chirality>,
    planes: Vec<_core::geometry::Plane>,
) -> RestraintGraph {
    match RestraintGraph::new(
        (0..5).map(|i| format!("atom{i}")).collect(),
        vec![false, false, false, false, true],
        bonds,
        angles,
        chirals,
        planes,
    ) {
        Ok(graph) => graph,
        Err(error) => unreachable!("Valid test graph: {error}"),
    }
}

fn sites() -> [[f64; 3]; 5] {
    [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
        [1.0, 1.0, 0.1],
    ]
}

#[test]
fn cutoff_is_strict_and_bond_slack_does_not_change_validation_score() {
    let g = graph(vec![(0, 1, 1.5, 100.0, 0.4, true)], vec![], vec![], vec![]);
    let at = g.evaluate(&sites(), 5.0).expect("Valid sites");
    assert!(at.findings.is_empty());
    let below = g.evaluate(&sites(), 4.999).expect("Valid sites");
    assert_eq!(below.findings[0].score, 5.0);
    assert!((below.findings[0].residual - 1.0).abs() < 1e-12);
    assert_eq!(below.hydrogen_failures, vec![(4, vec![])]);
}

#[test]
fn angle_and_chirality_use_reference_units_and_signs() {
    let g = graph(
        vec![],
        vec![(1, 0, 2, 100.0, 1.0)],
        vec![
            (0, 1, 2, 3, -1.0, 100.0, false, false),
            (0, 1, 2, 3, -1.0, 100.0, false, true),
        ],
        vec![],
    );
    let result = g.evaluate(&sites(), 5.0).expect("Valid sites");
    assert_eq!(result.findings.len(), 2);
    assert_eq!(result.findings[0].model, 90.0);
    assert_eq!(result.findings[1].model, 1.0);
    assert_eq!(result.findings[1].score, 20.0);
    assert!(result.counts.contains(&("chiralities", 2)));
    let g = graph(
        vec![],
        vec![],
        vec![(0, 1, 2, 3, 1.0, 100.0, true, false)],
        vec![],
    );
    let mut xyz = sites();
    xyz[3][2] = -1.0;
    assert!(
        g.evaluate(&xyz, 5.0)
            .expect("Valid sites")
            .findings
            .is_empty()
    );
}

#[test]
fn weighted_plane_is_translation_and_rotation_invariant() {
    let g = graph(
        vec![],
        vec![],
        vec![],
        vec![(vec![0, 1, 2, 4], vec![10000.0, 100.0, 400.0, 900.0])],
    );
    let result = g.evaluate(&sites(), 0.0).expect("Valid sites");
    let f = &result.findings[0];
    assert!(f.score > 0.0 && f.residual > 0.0);
    let transformed = sites().map(|p| [p[2] + 10.0, p[0] - 23.0, p[1] + 2.0]);
    let transformed_result = g.evaluate(&transformed, 0.0).expect("Valid sites");
    assert!((f.score - transformed_result.findings[0].score).abs() < 1e-11);
    assert!((f.residual - transformed_result.findings[0].residual).abs() < 1e-11);
}

#[test]
fn hydrogen_parent_count_uses_nonreporting_bonds_and_excludes_hydrogen_parents() {
    let g = graph(vec![(0, 4, 1.0, 100.0, 0.0, false)], vec![], vec![], vec![]);
    let result = g.evaluate(&sites(), 0.0).expect("Valid sites");
    assert!(result.findings.is_empty());
    assert!(result.hydrogen_failures.is_empty());
    let g = graph(
        vec![
            (0, 4, 1.0, 100.0, 0.0, false),
            (1, 4, 1.0, 100.0, 0.0, true),
        ],
        vec![],
        vec![],
        vec![],
    );
    assert_eq!(
        g.evaluate(&sites(), 5.0)
            .expect("Valid sites")
            .hydrogen_failures,
        vec![(4, vec![0, 1])]
    );
}

#[test]
fn malformed_and_undefined_geometry_fail_closed() {
    assert!(
        RestraintGraph::new(
            vec!["a".into()],
            vec![false],
            vec![(0, 0, 1.0, 1.0, 0.0, true)],
            vec![],
            vec![],
            vec![]
        )
        .is_err()
    );
    let g = graph(vec![], vec![(1, 0, 2, 90.0, 1.0)], vec![], vec![]);
    let mut xyz = sites();
    xyz[1] = xyz[0];
    assert!(g.evaluate(&xyz, 5.0).is_err());
    xyz[1][0] = f64::NAN;
    assert!(g.evaluate(&xyz, 5.0).is_err());
    let g = graph(vec![], vec![], vec![], vec![(vec![0, 1, 2], vec![1.0; 3])]);
    xyz = sites();
    xyz[2] = [2.0, 0.0, 0.0];
    assert!(g.evaluate(&xyz, 5.0).is_err());
    let g = graph(
        vec![],
        vec![],
        vec![],
        vec![(vec![0, 1, 2, 3], vec![1.0; 4])],
    );
    xyz[..4].copy_from_slice(&[
        [1.0, 1.0, 1.0],
        [-1.0, -1.0, 1.0],
        [1.0, -1.0, -1.0],
        [-1.0, 1.0, -1.0],
    ]);
    assert!(g.evaluate(&xyz, 5.0).is_err());
}
