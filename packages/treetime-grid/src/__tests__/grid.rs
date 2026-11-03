#[cfg(test)]
mod tests {
  use crate::{Grid, GridFn, MaxGridPoints, coordinates_coincide};
  use eyre::Report;
  use ndarray::{Array1, array};
  use pretty_assertions::assert_eq;
  use proptest::prelude::*;
  use rstest::rstest;
  use treetime_utils::assert_error;

  const ULP_2000: f64 = 1024.0 * f64::EPSILON;

  const ULP_2048: f64 = 2048.0 * f64::EPSILON;

  #[rustfmt::skip]
  #[rstest]
  #[case::at_x_min(0.0, 1.0, 5, 0.0, 0)]
  #[case::at_x_max(0.0, 1.0, 5, 4.0, 3)]
  #[case::at_midpoint(0.0, 1.0, 5, 2.0, 2)]
  #[case::between_0_1(0.0, 1.0, 5, 0.5, 0)]
  #[case::between_1_2(0.0, 1.0, 5, 1.5, 1)]
  #[case::below_min_small(0.0, 1.0, 5, -0.5, 0)]
  #[case::below_min_large(0.0, 1.0, 5, -100.0, 0)]
  #[case::below_positive_min(5.0, 1.0, 5, 4.0, 0)]
  #[case::above_max_small(0.0, 1.0, 5, 4.5, 3)]
  #[case::above_max_large(0.0, 1.0, 5, 100.0, 3)]
  #[trace]
  fn test_find_interval_index(
    #[case] x_min: f64,
    #[case] dx: f64,
    #[case] n_points: usize,
    #[case] query: f64,
    #[case] expected: usize,
  ) -> Result<(), Report> {
    let grid = Grid::from_start_dx(x_min, dx, n_points)?;
    let actual = grid.find_interval_index(query);
    assert_eq!(expected, actual);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::equal(                          ( 2.0,           2.0,            1.0),           true)]
  #[case::within_rounding_of_the_scale(   ( 0.0,           8e-16,          1.0),           true)]
  #[case::beyond_rounding_of_the_scale(   ( 0.0,           1e-15,          1.0),           false)]
  #[case::opposite_signs_within_rounding( (-4e-16,         4e-16,          1.0),           true)]
  #[case::tiny_values_at_a_tiny_scale(    ( 0.0,           1e-17,          1e-6),          false)]
  #[case::date_scale(                     ( 2010.0,        2010.0 + 1e-12, 2010.0),        true)]
  #[case::equal_infinities(               ( f64::INFINITY, f64::INFINITY,  1.0),           true)]
  #[case::infinity_and_finite(            ( f64::INFINITY, 1.0,            f64::INFINITY), false)]
  #[case::nan(                            ( f64::NAN,      f64::NAN,       1.0),           false)]
  #[trace]
  fn test_coordinates_coincide(#[case] (a, b, scale): (f64, f64, f64), #[case] expected: bool) {
    let actual = coordinates_coincide(a, b, scale);
    assert_eq!(expected, actual);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::one_ulp_cells(          (2000.0, 2000.0 + 4.0 * ULP_2000),  5, array![2000.0, 2000.0 + ULP_2000,       2000.0 + 2.0 * ULP_2000, 2000.0 + 3.0 * ULP_2000,  2000.0 + 4.0 * ULP_2000])]
  #[case::cells_rounded_unevenly( (2048.0, 2048.0 + 17.0 * ULP_2048), 5, array![2048.0, 2048.0 + 4.0 * ULP_2048, 2048.0 + 8.0 * ULP_2048, 2048.0 + 13.0 * ULP_2048, 2048.0 + 17.0 * ULP_2048])]
  #[trace]
  fn test_grid_from_range_n_points_accepts_strictly_increasing_coordinates(
    #[case] (x_min, x_max): (f64, f64),
    #[case] n_points: usize,
    #[case] expected: Array1<f64>,
  ) -> Result<(), Report> {
    let actual = Grid::from_range_n_points(x_min, x_max, n_points)?.to_array();
    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  fn test_grid_slice_of_an_unevenly_rounded_grid_is_accepted() -> Result<(), Report> {
    let grid = Grid::from_range_n_points(2048.0, 2048.0 + 17.0 * ULP_2048, 5)?;
    let expected = array![2048.0, 2048.0 + 4.0 * ULP_2048];
    let actual = grid.slice(0..=1)?.to_array();
    assert_eq!(expected, actual);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::one_point(              (0.0,                   1.0),                    1, "Grid must have at least 2 points, got 1")]
  #[case::nan_x_min(              (f64::NAN,              1.0),                    3, "Grid x_min (NaN) must be finite")]
  #[case::nan_x_max(              (0.0,                   f64::NAN),               3, "Grid x_max (NaN) must be finite")]
  #[case::negative_infinite_x_min((f64::NEG_INFINITY,     1.0),                    3, "Grid x_min (-inf) must be finite")]
  #[case::infinite_x_max(         (0.0,                   f64::INFINITY),          3, "Grid x_max (inf) must be finite")]
  #[case::reversed(               (1.0,                   0.0),                    3, "x_max (0.0) must be greater than x_min (1.0)")]
  #[case::empty(                  (1.0,                   1.0),                    3, "x_max (1.0) must be greater than x_min (1.0)")]
  #[case::span_overflows(         (-1e308,                1e308),                  3, "Grid of 3 points over [-1e308, 1e308] has spacing inf, which is not a finite positive number")]
  #[case::spacing_underflows(     (0.0,                   5e-324),                 3, "Grid of 3 points over [0.0, 5e-324] has spacing 0.0, which is not a finite positive number")]
  #[case::repeated_first_point(   (2000.0,                2000.0 + ULP_2000),      3, "Grid of 3 points over [2000.0, 2000.0000000000002] with spacing 1.1368683772161603e-13 does not increase strictly: x[1] = 2000.0 is not above x[0] = 2000.0")]
  #[case::repeated_stored_end(    (2000.0 + ULP_2000,     2000.0 + 2.0 * ULP_2000), 3, "Grid of 3 points over [2000.0000000000002, 2000.0000000000005] with spacing 1.1368683772161603e-13 does not increase strictly: x[2] = 2000.0000000000005 is not above x[1] = 2000.0000000000005")]
  #[case::inexact_point_count(    (0.0,                   1.0),                    (1 << 53) + 2, "Grid of 9007199254740994 points over [0.0, 1.0] has more points than 9007199254740991, the largest count with exact coordinates")]
  #[case::first_inexact_count(    (0.0,                   1.0),                    1 << 53,       "Grid of 9007199254740992 points over [0.0, 1.0] has more points than 9007199254740991, the largest count with exact coordinates")]
  #[trace]
  fn test_grid_from_range_n_points_rejects_invalid_ranges(#[case] (x_min, x_max): (f64, f64), #[case] n_points: usize, #[case] expected: &str) {
    assert_error!(Grid::from_range_n_points(x_min, x_max, n_points), expected);
  }

  #[test]
  fn test_grid_from_start_dx_rejects_inexact_point_count() {
    assert_error!(
      Grid::from_start_dx(0.0, 1.0, (1 << 53) + 2),
      "Grid of 9007199254740994 points from 0.0 with spacing 1.0 has more points than 9007199254740991, the largest count with exact coordinates"
    );
  }

  #[test]
  #[should_panic(expected = "index 5 is outside the grid of 5 points")]
  fn test_grid_x_at_panics_past_the_last_point() {
    Grid::from_range_n_points(0.0, 1.0, 5).unwrap().x_at(5);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::nan_x_min(     (f64::NAN,      1.0),               "Grid x_min (NaN) must be finite")]
  #[case::infinite_x_min((f64::INFINITY, 1.0),               "Grid x_min (inf) must be finite")]
  #[case::nan_dx(        (0.0,           f64::NAN),          "Grid spacing dx (NaN) must be a finite positive number")]
  #[case::infinite_dx(   (0.0,           f64::INFINITY),     "Grid spacing dx (inf) must be a finite positive number")]
  #[case::zero_dx(       (0.0,           0.0),               "Grid spacing dx (0.0) must be a finite positive number")]
  #[case::negative_dx(   (0.0,           -1.0),              "Grid spacing dx (-1.0) must be a finite positive number")]
  #[case::end_overflows( (1e308,         1e308),             "Grid x_max (inf) must be finite")]
  #[trace]
  fn test_gridfn_from_start_dx_values_rejects_invalid_inputs(#[case] (x_min, dx): (f64, f64), #[case] expected: &str) {
    assert_error!(GridFn::from_start_dx_values(x_min, dx, array![0.0, 1.0, 2.0]), expected);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::nan_x_min(     (f64::NAN, 1.0,           0.1),      "Grid x_min (NaN) must be finite")]
  #[case::infinite_x_max((0.0,      f64::INFINITY, 0.1),      "Grid x_max (inf) must be finite")]
  #[case::nan_x_max(     (0.0,      f64::NAN,      0.1),      "Grid x_max (NaN) must be finite")]
  #[case::nan_dx(        (0.0,      1.0,           f64::NAN), "Grid spacing dx (NaN) must be a finite positive number")]
  #[case::zero_dx(       (0.0,      1.0,           0.0),      "Grid spacing dx (0.0) must be a finite positive number")]
  #[case::reversed(      (1.0,      0.0,           0.1),      "x_max (0.0) must be greater than x_min (1.0)")]
  #[trace]
  fn test_grid_from_range_dx_rejects_invalid_inputs(#[case] (x_min, x_max, dx): (f64, f64, f64), #[case] expected: &str) {
    assert_error!(Grid::from_range_dx(x_min, x_max, dx, MaxGridPoints::default()), expected);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::nan_interior(                 array![0.0, f64::NAN, 2.0],                                                                                  "x array is not a grid: x[1] is NaN, which is not finite")]
  #[case::infinite_end(                 array![0.0, 1.0, f64::INFINITY],                                                                             "x array is not a grid: x[2] is inf, which is not finite")]
  #[case::repeated_coordinates(         array![2000.0, 2000.0 + 12.0 * ULP_2000, 2000.0 + 12.0 * ULP_2000, 2000.0 + 24.0 * ULP_2000],                "x array is not a grid: x[2] = 2000.0000000000027 is not above x[1] = 2000.0000000000027")]
  #[case::reversed(                     array![2.0, 1.0, 0.0],                                                                                       "x array is not a grid: x[1] = 1.0 is not above x[0] = 2.0")]
  #[case::nonuniform(                   array![0.0, 1.0, 3.0],                                                                                       "x array is not a uniform grid: x[1] is 1.0, but the uniform grid from x[0] = 0.0 to x[2] = 3.0 has 1.5 there")]
  #[case::nonuniform_below_end_rounding(array![2000.0, 2000.0 + 11.0 * ULP_2000, 2000.0 + 13.0 * ULP_2000, 2000.0 + 24.0 * ULP_2000],                "x array is not a uniform grid: x[1] is 2000.0000000000025, but the uniform grid from x[0] = 2000.0 to x[3] = 2000.0000000000055 has 2000.0000000000018 there")]
  #[case::nonuniform_at_ulp_spacing(    array![2000.0, 2000.0 + 2.0 * ULP_2000, 2000.0 + 3.0 * ULP_2000, 2000.0 + 4.0 * ULP_2000],                   "x array is not a uniform grid: x[1] is 2000.0000000000005, but the uniform grid from x[0] = 2000.0 to x[3] = 2000.000000000001 has 2000.0000000000002 there")]
  #[trace]
  fn test_grid_from_array_rejects_coordinates_it_cannot_represent(#[case] x: Array1<f64>, #[case] expected: &str) {
    assert_error!(Grid::from_array(&x), expected);
  }

  #[test]
  fn test_grid_from_array_accepts_the_rounding_of_a_generator() -> Result<(), Report> {
    let x = array![
      2010.3,
      2010.301,
      2010.302,
      2010.3029999999999,
      2010.3039999999999,
      2010.305,
      2010.306
    ];
    let grid = Grid::from_array(&x)?;
    let expected = (2010.3, 2010.306, 7);
    let actual = (grid.x_min(), grid.x_max(), grid.n_points());
    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  fn test_grid_shift_keeps_coordinates_two_ulps_apart() -> Result<(), Report> {
    let grid_fn = GridFn::from_range_values((0.0, 1e-12), array![0.0, 1.0, 2.0])?;
    let expected = array![2000.0, 2000.0 + 2.0 * ULP_2000, 2000.0 + 4.0 * ULP_2000];
    let actual = grid_fn.shift_x(2000.0)?.x();
    assert_eq!(expected, actual);
    Ok(())
  }

  proptest! {
    #[test]
    fn test_prop_grid_from_range_n_points_accepts_exactly_strictly_increasing_coordinates(
      exponent in -40_i32..40,
      mantissa in 1.0_f64..2.0,
      negative in any::<bool>(),
      span_ulps in 1_u32..400,
      n_points in 2_usize..64,
    ) {
      let binade = 2.0_f64.powi(exponent);
      let x_min = if negative { -mantissa * binade } else { mantissa * binade };
      let x_max = x_min + f64::from(span_ulps) * f64::EPSILON * binade;
      prop_assume!(x_max > x_min);
      let expected = helpers::coordinates_increase_strictly(x_min, x_max, n_points);
      let actual = Grid::from_range_n_points(x_min, x_max, n_points).is_ok();
      prop_assert_eq!(expected, actual);
    }
  }

  mod helpers {
    use itertools::Itertools;
    use treetime_utils::cast::ToF64Exact;

    pub(super) fn coordinates_increase_strictly(x_min: f64, x_max: f64, n_points: usize) -> bool {
      let dx = (x_max - x_min) / (n_points - 1).to_f64_exact();
      (0..n_points - 1)
        .map(|i| x_min + dx * i.to_f64_exact())
        .chain([x_max])
        .tuple_windows()
        .all(|(left, right)| left < right)
    }
  }
}
