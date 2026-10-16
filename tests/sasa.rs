//! Analytic and boundary checks for selected-atom surface areas.

use _core::geometry::sasa;
use std::f64::consts::PI;

#[test]
fn isolated_contained_tangent_and_coincident_spheres() -> Result<(), String> {
    let full = 4.0 * PI * 3.1_f64.powi(2);
    let run = |xyz: &[[f64; 3]], radii: &[f64]| sasa(xyz, radii, &[0], 1.4, 256);
    assert!((run(&[[0.; 3]], &[1.7])?[0] - full).abs() < 1e-10);
    assert_eq!(run(&[[0.; 3], [0.; 3]], &[1.7, 1.7])?, [0.]);
    assert_eq!(run(&[[0.; 3], [0.1, 0., 0.]], &[1.7, 2.])?, [0.]);
    assert!((run(&[[0.; 3], [6.2, 0., 0.]], &[1.7, 1.7])?[0] - full).abs() < 1e-10);
    assert!((run(&[[0.; 3], [0.1, 0., 0.]], &[1.7, 1.])?[0] - full).abs() < 1e-10);
    Ok(())
}

#[test]
fn partial_equal_spheres_match_analytic_area() -> Result<(), String> {
    let result = sasa(&[[0.; 3], [3.1, 0., 0.]], &[1.7; 2], &[1, 0], 1.4, 100_000)?;
    let expected = 3.0 * PI * 3.1_f64.powi(2);
    assert!(result.iter().all(|a| (a - expected).abs() < 0.03));
    Ok(())
}

#[test]
fn tiny_radii_and_large_negative_coordinates_do_not_overflow_cells() -> Result<(), String> {
    let result = sasa(&[[-1e9; 3], [1e9; 3]], &[1e-100; 2], &[1, 0], 0., 1)?;
    assert!(result.iter().all(|a| a.is_finite() && *a > 0.));
    Ok(())
}

#[test]
fn rejects_invalid_arrays_and_options_even_for_empty_selection() -> Result<(), String> {
    assert!(sasa(&[[0.; 3]], &[], &[], 1.4, 256).is_err());
    assert!(sasa(&[[f64::NAN; 3]], &[1.7], &[], 1.4, 256).is_err());
    assert!(sasa(&[[0.; 3]], &[0.], &[0], 1.4, 256).is_err());
    assert!(sasa(&[[0.; 3]], &[1.7], &[1], 1.4, 256).is_err());
    assert!(sasa(&[[0.; 3]], &[1.7], &[0, 0], 1.4, 256).is_err());
    assert!(sasa(&[], &[], &[], -1., 256).is_err());
    assert!(sasa(&[], &[], &[], 1.4, 0).is_err());
    assert!(sasa(&[], &[], &[], 1.4, 1_000_001).is_err());
    assert!(sasa(&[], &[], &[], 1.4, 256)?.is_empty());
    Ok(())
}
