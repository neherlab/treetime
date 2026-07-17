#[cfg(test)]
mod tests {
  use crate::DistributionPlain as Distribution;
  use crate::distribution_ops::log_cost::distribution_apply_neg_log_weight;
  use approx::assert_abs_diff_eq;
  use ndarray::array;

  /// Empty passes through unchanged.
  #[test]
  fn test_log_cost_empty_passthrough() {
    let actual = distribution_apply_neg_log_weight(&Distribution::Empty, |_| Ok(1.0)).unwrap();
    assert_eq!(Distribution::Empty, actual);
  }

  /// A Formula has no grid and is rejected.
  #[test]
  fn test_log_cost_rejects_formula() {
    let formula = Distribution::Formula(crate::DistributionFormula::new(|_| Ok(1.0), 0.0, 10.0));
    let error = distribution_apply_neg_log_weight(&formula, |_| Ok(1.0)).unwrap_err();
    assert!(error.to_string().contains("concrete Point, Range, or Function"));
  }

  /// Log-space accumulation preserves the product peak when the weight spans
  /// more than ~700 across the grid.
  ///
  /// The product peak sits at t=1, where the weight is largest (800) but the
  /// message amplitude is largest too. A naive `amplitude * exp(min_weight -
  /// weight)` shift evaluates `exp(-800)`, which underflows to exactly 0.0
  /// before multiplication, destroying the peak. Shifting by the combined
  /// minimum `weight - ln(amplitude)` keeps the peak at unit magnitude.
  ///
  /// Analytical oracle: the normalized product `amplitude_i * exp(-weight_i)`
  /// peaks at t=1, and `product(0) / product(1) = exp((ln a0 - w0) - (ln a1 -
  /// w1)) = exp((-405 - 0) - (405 - 800)) = exp(-10)`.
  // Exact float equality is intentional: the naive shift underflows to exactly
  // 0.0, and the normalized peak is exactly 1.0 by construction.
  #[allow(clippy::float_cmp)]
  #[test]
  fn test_log_cost_preserves_peak_under_wide_weight_span() {
    // Amplitudes e^-405, e^405, e^-405 are all finite, representable f64 values
    // (e^405 ~ 8e175, e^-405 ~ 1e-176); their product with exp(-weight) is not.
    let amplitudes = array![(-405.0_f64).exp(), (405.0_f64).exp(), (-405.0_f64).exp()];
    let times = array![0.0, 1.0, 2.0];
    let distribution = Distribution::function(times.clone(), amplitudes).unwrap();

    let weight = |time: f64| Ok(if (time - 1.0).abs() < 1e-9 { 800.0 } else { 0.0 });
    let actual = distribution_apply_neg_log_weight(&distribution, weight).unwrap();

    // The naive cost-only shift zeroes the peak: exp(-800) underflows, so
    // e^405 * exp(-800) = e^405 * 0.0 = 0.0.
    let naive_peak = (405.0_f64).exp() * (0.0_f64 - 800.0).exp();
    assert_eq!(0.0, naive_peak, "the naive shift must underflow the peak to zero");

    // The log-space result keeps the grid, peaks at t=1, and stays finite.
    assert_eq!(times, actual.t());
    assert!(matches!(actual, Distribution::Function(_)));
    let peak = actual.eval(1.0).unwrap();
    assert_eq!(1.0, peak, "normalized peak must be 1.0 at t=1");
    assert!(actual.eval(0.0).unwrap() < peak);
    assert!(actual.eval(2.0).unwrap() < peak);

    // Off-peak ratio matches the analytical value exp(-10). The tolerance
    // absorbs the ln/exp round-trip on amplitudes of magnitude e^405.
    assert_abs_diff_eq!((-10.0_f64).exp(), actual.eval(0.0).unwrap(), epsilon = 1e-15);
    assert_abs_diff_eq!((-10.0_f64).exp(), actual.eval(2.0).unwrap(), epsilon = 1e-15);
  }

  /// Zero-amplitude grid points contribute no mass (their negative-log value is
  /// +inf and is excluded from the shift).
  // Exact float equality is intentional: a zero-amplitude point yields exactly 0.0.
  #[allow(clippy::float_cmp)]
  #[test]
  fn test_log_cost_zero_amplitude_stays_zero() {
    let amplitudes = array![0.0, 1.0, 0.5];
    let times = array![0.0, 1.0, 2.0];
    let distribution = Distribution::function(times, amplitudes).unwrap();

    let actual = distribution_apply_neg_log_weight(&distribution, |_| Ok(0.0)).unwrap();

    assert_eq!(0.0, actual.eval(0.0).unwrap());
    assert!(actual.eval(1.0).unwrap() > 0.0);
  }
}
