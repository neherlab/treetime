#[cfg(test)]
mod tests {
  use crate::{Grid, GridFn, coordinates_coincide};
  use eyre::Report;
  use ndarray::{Array1, array};
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime_utils::assert_error;

  const ULP_2000: f64 = 1024.0 * f64::EPSILON;

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
  #[case::spacing_below_one_ulp(  array![2000.0, 2000.0 + 2.0 * ULP_2000, 2000.0 + 3.0 * ULP_2000, 2000.0 + 4.0 * ULP_2000], "Grid of 4 points over [2000.0, 2000.000000000001] has spacing 3.0316490059097606e-13, not above the rounding resolution 1.7763568394002513e-12 of its coordinates")]
  #[case::nonuniform(             array![0.0, 1.0, 3.0],                                                                        "x array is not a uniform grid: x[1] is 1.0, but the uniform grid from x[0] = 0.0 to x[2] = 3.0 has 1.5 there")]
  #[trace]
  fn test_grid_from_array_rejects_coordinates_it_cannot_represent(#[case] x: Array1<f64>, #[case] expected: &str) {
    assert_error!(Grid::from_array(&x), expected);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::range_below_resolution((2000.0, 2000.0 + 4.0 * ULP_2000), 5, "Grid of 5 points over [2000.0, 2000.000000000001] has spacing 2.2737367544323206e-13, not above the rounding resolution 1.7763568394002513e-12 of its coordinates")]
  #[trace]
  fn test_grid_from_range_n_points_rejects_spacing_below_rounding_resolution(#[case] (x_min, x_max): (f64, f64), #[case] n_points: usize, #[case] expected: &str) {
    assert_error!(Grid::from_range_n_points(x_min, x_max, n_points), expected);
  }

  #[test]
  fn test_grid_shift_that_collapses_the_spacing_is_error() -> Result<(), Report> {
    let grid_fn = GridFn::from_range_values((0.0, 1e-12), array![0.0, 1.0, 2.0])?;
    assert_error!(
      grid_fn.shift_x(2000.0),
      "Grid of 3 points over [2000.0, 2000.000000000001] has spacing 4.547473508864641e-13, not above the rounding resolution 1.7763568394002513e-12 of its coordinates"
    );
    Ok(())
  }
}
