#[cfg(test)]
mod tests {
  use crate::Grid;
  use eyre::Report;
  use pretty_assertions::assert_eq;
  use rstest::rstest;

  #[rustfmt::skip]
  #[rstest]
  // Interior points
  #[case::at_x_min(0.0, 1.0, 5, 0.0, 0)]
  #[case::at_x_max(0.0, 1.0, 5, 4.0, 3)]
  #[case::at_midpoint(0.0, 1.0, 5, 2.0, 2)]
  #[case::between_0_1(0.0, 1.0, 5, 0.5, 0)]
  #[case::between_1_2(0.0, 1.0, 5, 1.5, 1)]
  // Below grid minimum - should return 0
  #[case::below_min_small(0.0, 1.0, 5, -0.5, 0)]
  #[case::below_min_large(0.0, 1.0, 5, -100.0, 0)]
  #[case::below_positive_min(5.0, 1.0, 5, 4.0, 0)]
  // Above grid maximum - should return last valid interval
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
}
