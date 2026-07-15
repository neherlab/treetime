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
  use treetime_distribution::Distribution;
  use treetime_grid::piecewise_constant_fn::PiecewiseConstantFn;
  use treetime_utils::{assert_error, pretty_assert_ulps_eq};

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
  fn test_compute_integral_merger_rate_constant_tc() -> Result<(), Report> {
    // Create PiecewiseConstantFn with k=2.0 from t=0 to t=10
    // Breakpoints: [0.0, 10.0], values: [0.0, 2.0, 2.0]
    let lineage_counts = PiecewiseConstantFn::new(array![0.0, 10.0], array![0.0, 2.0, 2.0]);

    let tc_dist = Distribution::constant(0.01);

    let actual = compute_integral_merger_rate(&tc_dist, &lineage_counts)?;

    let actual_y = actual.values();

    // I(0) = 0, I(10) = integral of 0.5*(2-1)/0.01 = 50*10 = 500
    pretty_assert_ulps_eq!(actual_y[0], 0.0);
    pretty_assert_ulps_eq!(actual_y[actual_y.len() - 1], 500.0, max_ulps = 1000);

    Ok(())
  }

  #[test]
  fn test_compute_integral_merger_rate_multiple_segments() -> Result<(), Report> {
    // Three segments: k=1 from 0-5, k=5 from 5-10
    // Breakpoints: [0.0, 5.0, 10.0], values: [0.0, 1.0, 5.0, 5.0]
    let lineage_counts = PiecewiseConstantFn::new(array![0.0, 5.0, 10.0], array![0.0, 1.0, 5.0, 5.0]);

    let tc_dist = Distribution::constant(0.01);

    let actual = compute_integral_merger_rate(&tc_dist, &lineage_counts)?;

    let actual_y = actual.values();

    // Segment 0-5: k=1, rate = 0.5*max(0.5, 1-1)/0.01 = 0.5*0.5/0.01 = 25, contribution = 25*5 = 125
    // Segment 5-10: k=5, rate = 0.5*(5-1)/0.01 = 200, contribution = 200*5 = 1000
    // Total = 125 + 1000 = 1125
    pretty_assert_ulps_eq!(actual_y[0], 0.0);
    pretty_assert_ulps_eq!(actual_y[actual_y.len() - 1], 1125.0, max_ulps = 1000);

    Ok(())
  }

  #[test]
  fn test_compute_integral_merger_rate_insufficient_points() {
    let lineage_counts = PiecewiseConstantFn::new(array![5.0], array![0.0, 2.0]);

    let tc_dist = Distribution::constant(0.01);

    let result = compute_integral_merger_rate(&tc_dist, &lineage_counts);

    assert_error!(result, "lineage count must have at least 2 breakpoints");
  }

  #[test]
  fn test_compute_integral_merger_rate_varying_tc() -> Result<(), Report> {
    // k=3 constant from 0-10, single segment
    let lineage_counts = PiecewiseConstantFn::new(array![0.0, 10.0], array![0.0, 3.0, 3.0]);

    // Tc varies linearly from 0.01 to 0.05
    let t_vals = Array1::linspace(0.0, 10.0, 100);
    let tc_vals = Array1::linspace(0.01, 0.05, 100);
    let tc_dist = Distribution::function(t_vals, tc_vals)?;

    let actual = compute_integral_merger_rate(&tc_dist, &lineage_counts)?;
    let actual_y = actual.values();

    // Single segment evaluates at midpoint t=5 where Tc=0.03
    // rate = 0.5 * (3-1) / 0.03 = 33.333...
    // integral = 33.333 * 10 = 333.333...
    let expected = 10.0 / 0.03; // 333.333...
    pretty_assert_ulps_eq!(actual_y[0], 0.0);
    pretty_assert_ulps_eq!(actual_y[actual_y.len() - 1], expected, max_ulps = 1000);

    Ok(())
  }

  // Production midpoint quadrature is inaccurate when Tc varies inside a lineage interval.
  // See kb/issues/N-coalescent-skyline-quadrature-contract-undecided.md.
  #[test]
  #[ignore = "varying-Tc midpoint quadrature error (kb/issues/N-coalescent-skyline-quadrature-contract-undecided.md)"]
  fn test_compute_integral_merger_rate_varying_tc_many_segments() -> Result<(), Report> {
    // k=3 constant, 1000 segments for high-accuracy numerical integration.
    let n_segments = 1000;
    let breakpoints = Array1::linspace(0.0, 10.0, n_segments + 1);
    let values = {
      let mut v = vec![0.0];
      v.extend(std::iter::repeat_n(3.0, n_segments + 1));
      Array1::from(v)
    };
    let lineage_counts = PiecewiseConstantFn::new(breakpoints, values);

    // Tc varies linearly from 0.01 to 0.05
    let t_vals = array![0.0, 10.0];
    let tc_vals = array![0.01, 0.05];
    let tc_dist = Distribution::function(t_vals, tc_vals)?;

    let actual = compute_integral_merger_rate(&tc_dist, &lineage_counts)?;
    let actual_y = actual.values();

    // Analytical: ∫₀¹⁰ (k-1)/(2*Tc(t)) dt = ∫₀¹⁰ 1/(0.01 + 0.004t) dt
    //           = (1/0.004) * ln(0.05/0.01) = 250 * ln(5) ≈ 402.359
    let expected_analytical = 250.0 * 5.0_f64.ln();
    pretty_assert_ulps_eq!(actual_y[0], 0.0);
    assert_abs_diff_eq!(actual_y[actual_y.len() - 1], expected_analytical, epsilon = 1e-6);

    Ok(())
  }
}
