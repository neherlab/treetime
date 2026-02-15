#[cfg(test)]
mod tests {
  use crate::commands::timetree::coalescent::integration::compute_integral_merger_rate;
  use crate::commands::timetree::coalescent::integration::compute_merger_rates;
  use crate::commands::timetree::coalescent::piecewise_constant_fn::PiecewiseConstantFn;
  use approx::assert_abs_diff_eq;
  use eyre::Report;
  use ndarray::Array1;
  use ndarray::array;
  use treetime_distribution::Distribution;
  use treetime_utils::pretty_assert_ulps_eq;

  #[test]
  fn test_merger_rates() -> Result<(), Report> {
    let k = array![2.0, 3.0, 4.0];
    let tc = array![0.001, 0.002, 0.003];

    let (actual_branch_rate, actual_total_rate) = compute_merger_rates(&k, &tc);

    let expected_branch_rate = array![500.0, 500.0, 500.0];
    let expected_total_rate = array![1000.0, 1500.0, 2000.0];

    pretty_assert_ulps_eq!(actual_branch_rate, expected_branch_rate);
    pretty_assert_ulps_eq!(actual_total_rate, expected_total_rate);

    Ok(())
  }

  #[test]
  fn test_merger_rates_edge_cases() -> Result<(), Report> {
    let k = array![1.0, 0.5, 2.0];
    let tc = array![0.001, 0.002, 0.003];

    let (actual_branch_rate, actual_total_rate) = compute_merger_rates(&k, &tc);

    let expected_branch_rate = array![250.0, 125.0, 166.66666666666666];
    let expected_total_rate = array![375.0, 187.5, 333.3333333333333];

    pretty_assert_ulps_eq!(actual_branch_rate, expected_branch_rate);
    pretty_assert_ulps_eq!(actual_total_rate, expected_total_rate);

    Ok(())
  }

  #[test]
  fn test_merger_rates_large_k() -> Result<(), Report> {
    let k = array![10.0, 20.0, 50.0];
    let tc = array![0.01, 0.01, 0.01];

    let (actual_branch_rate, actual_total_rate) = compute_merger_rates(&k, &tc);

    let expected_branch_rate = array![450.0, 950.0, 2450.0];
    let expected_total_rate = array![4500.0, 19000.0, 122500.0];

    pretty_assert_ulps_eq!(actual_branch_rate, expected_branch_rate);
    pretty_assert_ulps_eq!(actual_total_rate, expected_total_rate);

    Ok(())
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

    assert!(result.is_err());
    assert!(result.unwrap_err().to_string().contains("at least 2 breakpoints"));
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

  #[test]
  fn test_compute_integral_merger_rate_varying_tc_many_segments() -> Result<(), Report> {
    // k=3 constant, many segments for better numerical accuracy
    // Use 11 breakpoints (10 segments) to approach analytical solution
    let breakpoints = array![0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0];
    let values = array![0.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0];
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
    // With 10 segments, numerical approximation should be within 2% of analytical
    assert_abs_diff_eq!(actual_y[actual_y.len() - 1], expected_analytical, epsilon = 10.0);

    Ok(())
  }
}
