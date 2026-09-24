#[cfg(test)]
mod tests {
  use crate::Distribution;
  use crate::policy::NegLog;
  use ndarray::Array1;
  use proptest::prelude::*;
  use treetime_utils::{prop_assert_array_abs_diff_eq, prop_assert_array_finite, prop_assert_array_nonneg};

  proptest! {
    #[test]
    #[allow(clippy::float_cmp, reason = "the peak is subtracted from itself, which is exactly zero")]
    fn test_prop_distribution_normalize_peak_is_zero(y in generators::neglog_values()) {
      let normalized = helpers::normalized_values(&y);
      prop_assert_eq!(0.0, normalized.iter().copied().fold(f64::INFINITY, f64::min));
    }

    #[test]
    fn test_prop_distribution_normalize_preserves_finite_nonneg(y in generators::neglog_values()) {
      let normalized = helpers::normalized_values(&y);
      prop_assert_array_finite!(normalized);
      prop_assert_array_nonneg!(normalized);
    }

    #[test]
    fn test_prop_distribution_normalize_preserves_differences(y in generators::neglog_values()) {
      let normalized = helpers::normalized_values(&y);
      let expected = &y - y[0];
      prop_assert_array_abs_diff_eq!(expected, &normalized - normalized[0], epsilon = 1e-12);
    }

    #[test]
    fn test_prop_distribution_normalize_invariant_under_common_offset(
      y in generators::neglog_values(),
      offset in -1e3_f64..1e3,
    ) {
      let expected = helpers::normalized_values(&y);
      let actual = helpers::normalized_values(&(&y + offset));
      prop_assert_array_abs_diff_eq!(expected, actual, epsilon = 1e-12);
    }
  }

  mod generators {
    use super::*;

    pub(super) fn neglog_values() -> impl Strategy<Value = Array1<f64>> {
      prop::collection::vec(-1e3_f64..1e3, 3..50).prop_map(Array1::from)
    }
  }

  mod helpers {
    use super::*;

    pub(super) fn normalized_values(y: &Array1<f64>) -> Array1<f64> {
      let x = Array1::linspace(0.0, 1.0, y.len());
      let distribution = Distribution::<NegLog>::function(x, y.clone()).unwrap();
      let Distribution::Function(f) = distribution.normalize().unwrap() else {
        panic!("a finite function normalizes to a function");
      };
      f.y().clone()
    }
  }
}
