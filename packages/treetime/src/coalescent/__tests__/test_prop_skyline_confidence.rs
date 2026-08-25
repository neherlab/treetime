#[cfg(test)]
mod tests {
  use crate::coalescent::__tests__::helpers::setup_graph;
  use crate::coalescent::skyline::{SkylineParams, optimize_skyline};
  use proptest::prelude::*;
  use treetime_utils::{prop_assert_array_finite, prop_assert_array_nonneg, prop_assert_array_positive};

  proptest! {
    #![proptest_config(ProptestConfig::with_cases(64))]

    #[test]
    fn test_prop_skyline_confidence_band_brackets_estimate(
      n_points in 1_usize..12,
      stiffness in 0.01_f64..10.0,
      n_std in 0.0_f64..4.0,
    ) {
      let graph = setup_graph().unwrap();
      let params = SkylineParams {
        n_points,
        stiffness,
        n_std,
        tolerance: 1e-10,
        max_iter: 1000,
      };

      let result = optimize_skyline(&graph, &params).unwrap();

      prop_assert_array_finite!(result.log_tc_variances);
      prop_assert_array_positive!(result.log_tc_variances);
      prop_assert_array_finite!(result.tc_lower_bounds);
      prop_assert_array_finite!(result.tc_upper_bounds);
      prop_assert_array_nonneg!(&result.tc_values - &result.tc_lower_bounds, epsilon = 1e-14);
      prop_assert_array_nonneg!(&result.tc_upper_bounds - &result.tc_values, epsilon = 1e-14);
    }
  }
}
