#[cfg(test)]
mod tests {
  use crate::coalescent::integration::compute_merger_rates;
  use crate::coalescent::integration::compute_merger_rates_scalar;
  use ndarray::Array1;
  use proptest::prelude::*;
  use treetime_utils::prop_assert_array_ulps_eq;

  proptest! {
    #![proptest_config(ProptestConfig::with_cases(100))]

    #[test]
    fn test_prop_integration_compute_merger_rates_array_matches_scalar(
      inputs in prop::collection::vec((-1e6_f64..1e6_f64, 1e-6_f64..1e6_f64), 0..64),
    ) {
      let k = Array1::from_iter(inputs.iter().map(|(k, _)| *k));
      let tc = Array1::from_iter(inputs.iter().map(|(_, tc)| *tc));
      let expected_per_lineage = Array1::from_iter(
        inputs.iter().map(|(k, tc)| compute_merger_rates_scalar(*k, *tc).per_lineage),
      );
      let expected_total = Array1::from_iter(
        inputs.iter().map(|(k, tc)| compute_merger_rates_scalar(*k, *tc).total),
      );

      let actual = compute_merger_rates(&k, &tc);

      prop_assert_array_ulps_eq!(expected_per_lineage, actual.per_lineage, max_ulps = 0);
      prop_assert_array_ulps_eq!(expected_total, actual.total, max_ulps = 0);
    }
  }
}
