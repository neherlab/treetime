#[cfg(test)]
mod tests {
  use crate::grid::Grid;
  use crate::grid_fn::GridFn;
  use crate::hard_approach_law::Side;
  use crate::max_grid_points::{GridPointLimitExceeded, MaxGridPoints};
  use eyre::Report;
  use ndarray::array;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use schemars::schema_for;
  use serde_json::json;
  use std::str::FromStr;
  use treetime_utils::assert_error;
  use treetime_utils::io::json::{JsonPretty, json_read_str, json_write_str};

  #[test]
  fn test_max_grid_points_default_is_one_million() {
    assert_eq!(1_000_000, MaxGridPoints::default().get());
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::minimum(       "1000",    1_000)]
  #[case::default_value( "1000000", 1_000_000)]
  #[case::padded(        " 2000 ",  2_000)]
  #[trace]
  fn test_max_grid_points_from_str_accepts(#[case] text: &str, #[case] expected: usize) {
    assert_eq!(expected, MaxGridPoints::from_str(text).unwrap().get());
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::below_minimum("999",  "the grid point limit must be at least 1000 points, got 999")]
  #[case::zero(         "0",    "the grid point limit must be at least 1000 points, got 0")]
  #[case::negative(     "-5",   "'-5' is not a whole number: invalid digit found in string")]
  #[case::fraction(     "1e6",  "'1e6' is not a whole number: invalid digit found in string")]
  #[trace]
  fn test_max_grid_points_from_str_rejects(#[case] text: &str, #[case] expected: &str) {
    assert_eq!(expected, MaxGridPoints::from_str(text).unwrap_err().to_string());
  }

  #[test]
  fn test_max_grid_points_serde_roundtrip_is_a_plain_number() -> Result<(), Report> {
    let limit = MaxGridPoints::new(1_000)?;
    let text = json_write_str(&limit, JsonPretty(false))?;
    assert_eq!(("1000", limit), (text.as_str(), json_read_str::<MaxGridPoints>(&text)?));
    Ok(())
  }

  #[test]
  fn test_max_grid_points_serde_rejects_below_minimum() {
    assert_error!(
      json_read_str::<MaxGridPoints>("999"),
      "When parsing JSON: Unexpected: invalid value: the grid point limit must be at least 1000 points, got 999 at line 1 column 1"
    );
  }

  #[test]
  fn test_max_grid_points_schema_states_the_minimum() {
    let expected = json!({
      "$schema": "https://json-schema.org/draft/2020-12/schema",
      "title": "MaxGridPoints",
      "type": "integer",
      "format": "uint",
      "minimum": 1000,
    });
    assert_eq!(expected, schema_for!(MaxGridPoints).to_value());
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::below_the_limit(   999.0,  999)]
  #[case::exactly_the_limit( 1000.0, 1_000)]
  #[trace]
  fn test_max_grid_points_point_count_accepts(#[case] required: f64, #[case] expected: usize) -> Result<(), Report> {
    let actual = MaxGridPoints::new(1_000)?.point_count(required, (0.0, 1.0), 0.5)?;
    assert_eq!(expected, actual);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::one_above_the_limit(1001.0,        "A grid over [0, 1] with spacing 0.5 needs 1001 points, more than the limit of 1000")]
  #[case::beyond_usize(       1e300,         "A grid over [0, 1] with spacing 0.5 needs more points than the limit of 1000")]
  #[case::infinite(           f64::INFINITY, "A grid over [0, 1] with spacing 0.5 needs more points than the limit of 1000")]
  #[case::not_a_number(       f64::NAN,      "Cannot build a grid over [0, 1] with spacing 0.5: the number of points is not a number")]
  #[case::negative(           -2.0,          "Cannot build a grid over [0, 1] with spacing 0.5: the number of points is -2")]
  #[trace]
  fn test_max_grid_points_point_count_rejects(#[case] required: f64, #[case] expected: &str) {
    assert_error!(
      MaxGridPoints::new(1_000).unwrap().point_count(required, (0.0, 1.0), 0.5),
      expected
    );
  }

  #[test]
  fn test_max_grid_points_exceeded_error_is_downcastable() {
    let report = MaxGridPoints::new(1_000)
      .unwrap()
      .point_count(1001.0, (0.0, 1.0), 0.5)
      .unwrap_err()
      .wrap_err("outer context");
    let expected = GridPointLimitExceeded {
      required: 1001.0,
      limit: MaxGridPoints::new(1_000).unwrap(),
      range: (0.0, 1.0),
      dx: 0.5,
    };
    let actual = report
      .chain()
      .find_map(|cause| cause.downcast_ref::<GridPointLimitExceeded>())
      .copied();
    assert_eq!(Some(expected), actual);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::exactly_the_limit((0.0, 999.0), 1.0, 1_000)]
  #[case::rounds_to_limit(  (0.0, 999.4), 1.0, 1_000)]
  #[trace]
  fn test_max_grid_points_grid_from_range_dx_accepts(
    #[case] (x_min, x_max): (f64, f64),
    #[case] dx: f64,
    #[case] expected: usize,
  ) -> Result<(), Report> {
    let grid = Grid::from_range_dx(x_min, x_max, dx, MaxGridPoints::new(1_000)?)?;
    assert_eq!(expected, grid.n_points());
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::one_above_the_limit((0.0, 1000.0),   1.0,    "A grid over [0, 1000] with spacing 1 needs 1001 points, more than the limit of 1000")]
  #[case::ratio_overflows(    (0.0, 1.0),      1e-320, "A grid over [0, 1] with spacing 1.0e-320 needs more points than the limit of 1000")]
  #[case::infinite_range(     (0.0, f64::MAX), 1e-10,  "A grid over [0, 1.79769e308] with spacing 1.0e-10 needs more points than the limit of 1000")]
  #[case::not_a_number(       (0.0, f64::NAN), 1.0,    "Cannot build a grid over [0, NaN] with spacing 1: the number of points is not a number")]
  #[trace]
  fn test_max_grid_points_grid_from_range_dx_rejects(
    #[case] (x_min, x_max): (f64, f64),
    #[case] dx: f64,
    #[case] expected: &str,
  ) {
    assert_error!(
      Grid::from_range_dx(x_min, x_max, dx, MaxGridPoints::new(1_000).unwrap()),
      expected
    );
  }

  #[test]
  fn test_max_grid_points_resample_range_dx_clamped_at_the_limit() -> Result<(), Report> {
    let grid_fn = GridFn::from_range_values((0.0, 1.0), array![0.0, 1.0])?;
    let resampled = grid_fn.resample_range_dx_clamped((0.0, 999.0), 1.0, Side::Left, MaxGridPoints::new(1_000)?)?;
    assert_eq!(1_000, resampled.len());
    Ok(())
  }

  #[test]
  fn test_max_grid_points_resample_range_dx_clamped_above_the_limit() {
    let grid_fn = GridFn::from_range_values((0.0, 1.0), array![0.0, 1.0]).unwrap();
    assert_error!(
      grid_fn.resample_range_dx_clamped((0.0, 1000.0), 1.0, Side::Left, MaxGridPoints::new(1_000).unwrap()),
      "A grid over [0, 1000] with spacing 1 needs 1001 points, more than the limit of 1000"
    );
  }
}
