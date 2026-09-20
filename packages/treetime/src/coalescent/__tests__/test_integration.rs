#[cfg(test)]
mod tests {
  use crate::coalescent::integration::compute_integral_merger_rate;
  use crate::coalescent::integration::compute_merger_rate_per_lineage_scalar;
  use crate::coalescent::integration::compute_merger_rate_total_scalar;
  use approx::assert_abs_diff_eq;
  use eyre::Report;
  use ndarray::Array1;
  use ndarray::array;
  use rstest::rstest;
  use treetime_distribution::{Distribution, DistributionFunction};
  use treetime_grid::piecewise_constant_fn::PiecewiseConstantFn;
  use treetime_utils::pretty_assert_ulps_eq;

  #[rustfmt::skip]
  #[rstest]
  #[case::below_clamp(1.0, 1.0, 0.25, 0.375)]
  #[case::at_clamp(   1.5, 2.0, 0.125, 0.1875)]
  #[case::above_clamp(3.0, 2.0, 0.5, 1.5)]
  #[trace]
  fn test_integration_compute_merger_rates_scalar_representable(
    #[case] k: f64,
    #[case] tc: f64,
    #[case] expected_per_lineage: f64,
    #[case] expected_total: f64,
  ) {
    pretty_assert_ulps_eq!(expected_per_lineage, compute_merger_rate_per_lineage_scalar(k, tc), max_ulps = 0);
    pretty_assert_ulps_eq!(expected_total, compute_merger_rate_total_scalar(k, tc), max_ulps = 0);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::below_clamp(1.0, 3.0, 1.0 / 12.0, 1.0 / 8.0)]
  #[case::at_clamp(   1.5, 3.0, 1.0 / 12.0, 1.0 / 8.0)]
  #[case::above_clamp(2.0, 3.0, 1.0 / 6.0,  1.0 / 3.0)]
  #[trace]
  fn test_integration_compute_merger_rates_scalar_ulps(
    #[case] k: f64,
    #[case] tc: f64,
    #[case] expected_per_lineage: f64,
    #[case] expected_total: f64,
  ) {
    pretty_assert_ulps_eq!(expected_per_lineage, compute_merger_rate_per_lineage_scalar(k, tc));
    pretty_assert_ulps_eq!(expected_total, compute_merger_rate_total_scalar(k, tc));
  }

  #[test]
  fn test_integration_compute_merger_rates_scalar_preserves_v0_extreme_ordering() {
    pretty_assert_ulps_eq!(
      2.781342323134e-309,
      compute_merger_rate_per_lineage_scalar(2.0, f64::MAX),
      max_ulps = 0
    );
    pretty_assert_ulps_eq!(
      5.562684646268003e-309,
      compute_merger_rate_total_scalar(2.0, f64::MAX),
      max_ulps = 0
    );

    pretty_assert_ulps_eq!(
      1.1250000000000002e308,
      compute_merger_rate_total_scalar(1.5e154, 1.0),
      max_ulps = 0
    );
  }

  #[test]
  fn test_integration_compute_merger_rates_scalar_propagates_nan_lineage_count() {
    assert!(compute_merger_rate_per_lineage_scalar(f64::NAN, 1.0).is_nan());
    assert!(compute_merger_rate_total_scalar(f64::NAN, 1.0).is_nan());
  }

  #[test]
  fn test_integration_constant_tc_accumulates_from_present_to_past() -> Result<(), Report> {
    let lineage_counts = PiecewiseConstantFn::new(array![2000.0, 2010.0], array![1.0, 2.0, 0.0]);

    let actual = compute_integral_merger_rate(&Distribution::constant(0.01), &lineage_counts)?;

    pretty_assert_ulps_eq!(actual.values()[0], 500.0, max_ulps = 1000);
    pretty_assert_ulps_eq!(actual.values()[1], 0.0);
    Ok(())
  }

  #[test]
  fn test_integration_multiple_segments() -> Result<(), Report> {
    let lineage_counts = PiecewiseConstantFn::new(array![2000.0, 2005.0, 2010.0], array![1.0, 1.0, 5.0, 0.0]);

    let actual = compute_integral_merger_rate(&Distribution::constant(0.01), &lineage_counts)?;

    pretty_assert_ulps_eq!(actual.values()[0], 1125.0, max_ulps = 1000);
    pretty_assert_ulps_eq!(actual.values()[1], 1000.0, max_ulps = 1000);
    pretty_assert_ulps_eq!(actual.values()[2], 0.0);
    Ok(())
  }

  #[test]
  fn test_integration_varying_tc_uses_calendar_midpoint() -> Result<(), Report> {
    let lineage_counts = PiecewiseConstantFn::new(array![2000.0, 2010.0], array![1.0, 3.0, 0.0]);
    let tc = Distribution::Function(DistributionFunction::from_range_values(
      (2000.0, 2010.0),
      Array1::linspace(0.01, 0.05, 100),
    )?);

    let actual = compute_integral_merger_rate(&tc, &lineage_counts)?;

    pretty_assert_ulps_eq!(actual.values()[0], 10.0 / 0.03, max_ulps = 1000);
    pretty_assert_ulps_eq!(actual.values()[1], 0.0);
    Ok(())
  }

  #[test]
  #[ignore = "varying-Tc midpoint quadrature error (kb/issues/N-coalescent-skyline-quadrature-contract-undecided.md)"]
  fn test_integration_varying_tc_converges_with_refined_lineage_grid() -> Result<(), Report> {
    let n_segments = 1000;
    let breakpoints = Array1::linspace(2000.0, 2010.0, n_segments + 1);
    let values = Array1::from_iter(
      std::iter::once(1.0)
        .chain(std::iter::repeat_n(3.0, n_segments))
        .chain(std::iter::once(0.0)),
    );
    let lineage_counts = PiecewiseConstantFn::new(breakpoints, values);
    let tc = Distribution::function(array![2000.0, 2010.0], array![0.01, 0.05])?;

    let actual = compute_integral_merger_rate(&tc, &lineage_counts)?;

    let expected = 250.0 * 5.0_f64.ln();
    assert_abs_diff_eq!(expected, actual.values()[0], epsilon = 1e-6);
    Ok(())
  }

  #[test]
  fn test_integration_rejects_insufficient_breakpoints() {
    let lineage_counts = PiecewiseConstantFn::new(array![2000.0], array![1.0, 0.0]);
    let error = compute_integral_merger_rate(&Distribution::constant(1.0), &lineage_counts).unwrap_err();
    assert!(error.to_string().contains("at least 2 breakpoints"));
  }

  #[test]
  fn test_integration_rejects_nonpositive_tc() {
    let lineage_counts = PiecewiseConstantFn::new(array![2000.0, 2010.0], array![1.0, 2.0, 0.0]);

    let error = compute_integral_merger_rate(&Distribution::constant(0.0), &lineage_counts).unwrap_err();

    assert!(error.to_string().contains("finite and positive"));
  }
}
