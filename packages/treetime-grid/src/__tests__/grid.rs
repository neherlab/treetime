#[cfg(test)]
mod tests {
  use crate::{Grid, coordinates_coincide};
  use eyre::Report;
  use pretty_assertions::assert_eq;
  use rstest::rstest;

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
}
