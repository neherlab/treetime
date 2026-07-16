#[cfg(test)]
mod tests {
  use crate::coalescent::integration::compute_integral_merger_rate;
  use crate::coalescent::integration::compute_merger_rates;
  use crate::coalescent::integration::compute_merger_rates_scalar;
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
    // Oracle: n=max(0.5, k-1), κ=n/(2Tc), and λ=n(n+1)/(2Tc),
    // matching TreeTime v0 merger_models.py:194-213.
    let actual = compute_merger_rates_scalar(k, tc);

    pretty_assert_ulps_eq!(expected_per_lineage, actual.per_lineage, max_ulps = 0);
    pretty_assert_ulps_eq!(expected_total, actual.total, max_ulps = 0);
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
    // Oracle: the analytical v0 formulas cited in the representable scalar cases.
    let actual = compute_merger_rates_scalar(k, tc);

    pretty_assert_ulps_eq!(expected_per_lineage, actual.per_lineage);
    pretty_assert_ulps_eq!(expected_total, actual.total);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::empty(        array![],                          array![],                          array![],                                     array![])]
  #[case::equal_lengths(array![1.0, 1.5, 2.0, 4.0, 50.0], array![1.0, 2.0, 3.0, 4.0, 10.0], array![0.25, 0.125, 1.0 / 6.0, 0.375, 2.45], array![0.375, 0.1875, 1.0 / 3.0, 1.5, 122.5])]
  #[case::broadcast_tc( array![1.0, 1.5, 2.0],             array![2.0],                       array![0.125, 0.125, 0.25],                    array![0.1875, 0.1875, 0.5])]
  #[case::broadcast_k(  array![2.0],                       array![1.0, 2.0, 4.0],             array![0.5, 0.25, 0.125],                    array![1.0, 0.5, 0.25])]
  #[trace]
  fn test_integration_compute_merger_rates_array_cases(
    #[case] k: Array1<f64>,
    #[case] tc: Array1<f64>,
    #[case] expected_per_lineage: Array1<f64>,
    #[case] expected_total: Array1<f64>,
  ) {
    let actual = compute_merger_rates(&k, &tc);

    // Oracle: element-wise evaluation of the v0 formulas cited above.
    pretty_assert_ulps_eq!(expected_per_lineage, actual.per_lineage);
    pretty_assert_ulps_eq!(expected_total, actual.total);
  }

  #[test]
  #[should_panic(expected = "ndarray: could not broadcast array")]
  fn test_integration_compute_merger_rates_rejects_incompatible_shapes() {
    let k = array![1.0, 2.0];
    let tc = array![1.0, 2.0, 3.0];

    compute_merger_rates(&k, &tc);
  }

  #[test]
  fn test_integration_compute_merger_rates_scalar_preserves_v0_extreme_ordering() {
    // Oracle: TreeTime v0 merger_models.py:194-213 evaluates 0.5 in the
    // numerator before division, avoiding denominator overflow for large Tc.
    let large_tc = compute_merger_rates_scalar(2.0, f64::MAX);
    pretty_assert_ulps_eq!(2.781342323134e-309, large_tc.per_lineage, max_ulps = 0);
    pretty_assert_ulps_eq!(5.562684646268003e-309, large_tc.total, max_ulps = 0);

    // The same ordering halves n before multiplying n(n+1), keeping this
    // representable case finite where multiplying n(n+1) first overflows.
    let large_k = compute_merger_rates_scalar(1.5e154, 1.0);
    pretty_assert_ulps_eq!(1.1250000000000002e308, large_k.total, max_ulps = 0);
  }

  #[test]
  fn test_integration_compute_merger_rates_scalar_propagates_nan_lineage_count() {
    // Oracle: numpy.maximum in TreeTime v0 merger_models.py:194-213
    // propagates a NaN lineage count into both merger rates.
    let actual = compute_merger_rates_scalar(f64::NAN, 1.0);

    assert!(actual.per_lineage.is_nan());
    assert!(actual.total.is_nan());
  }

  #[test]
  fn test_integration_constant_tc_accumulates_from_present_to_past() -> Result<(), Report> {
    let lineage_counts = PiecewiseConstantFn::new(array![2000.0, 2010.0], array![1.0, 2.0, 0.0]);

    let actual = compute_integral_merger_rate(&Distribution::constant(0.01), &lineage_counts)?;

    // κ=50/year over ten years, so H(2000)=500 and H(2010)=0.
    pretty_assert_ulps_eq!(actual.values()[0], 500.0, max_ulps = 1000);
    pretty_assert_ulps_eq!(actual.values()[1], 0.0);
    Ok(())
  }

  #[test]
  fn test_integration_multiple_segments() -> Result<(), Report> {
    let lineage_counts = PiecewiseConstantFn::new(array![2000.0, 2005.0, 2010.0], array![1.0, 1.0, 5.0, 0.0]);

    let actual = compute_integral_merger_rate(&Distribution::constant(0.01), &lineage_counts)?;

    // κ=25/year for five years plus κ=200/year for five years.
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

    // Current skyline contract: one midpoint evaluation, Tc(2005)=0.03.
    pretty_assert_ulps_eq!(actual.values()[0], 10.0 / 0.03, max_ulps = 1000);
    pretty_assert_ulps_eq!(actual.values()[1], 0.0);
    Ok(())
  }

  // Production midpoint quadrature is inaccurate when Tc varies inside a lineage interval.
  // See kb/issues/N-coalescent-skyline-quadrature-contract-undecided.md.
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

    // Analytical oracle: ∫ 1/(0.01+0.004t) dt after shifting t to [0,10].
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
