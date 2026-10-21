#[cfg(test)]
mod tests {
  use crate::substitution_model::SubstitutionModel;
  use generators::{gen_pi, gen_w};
  use proptest::{prop_assert_eq, proptest};
  use treetime_utils::sync::random::get_random_number_generator;
  use treetime_utils::{prop_assert_array_ulps_eq, prop_assert_ulps_eq};

  proptest! {
    #[test]
    fn test_prop_substitution_model_random_has_unit_average_rate(seed in 0..u64::MAX, mu in 0.0..10.0_f64) {
      let model = SubstitutionModel::random(mu, &mut get_random_number_generator(seed)).unwrap();
      let pi = model.pi();
      let average_rate = pi.dot(&model.w().dot(pi));
      prop_assert_ulps_eq!(1.0, average_rate, max_ulps = 17);
    }

    #[test]
    fn test_prop_substitution_model_random_is_symmetric_with_zero_diagonal(seed in 0..u64::MAX) {
      let model = SubstitutionModel::random(1.0, &mut get_random_number_generator(seed)).unwrap();
      let w = model.w();
      prop_assert_eq!((w.t(), vec![0.0; 4]), (w.view(), w.diag().to_vec()));
    }

    #[test]
    fn test_prop_substitution_model_random_pi_is_a_distribution(seed in 0..u64::MAX) {
      let model = SubstitutionModel::random(1.0, &mut get_random_number_generator(seed)).unwrap();
      let positive = model.pi().iter().all(|&value| value > 0.0);
      prop_assert_eq!(true, positive);
      prop_assert_ulps_eq!(1.0, model.pi().sum(), max_ulps = 7);
    }

    #[test]
    fn test_prop_substitution_model_common_factor_of_pi_leaves_model_unchanged(
      pi in gen_pi(),
      w in gen_w(),
      exponent in -300.0..308.0_f64,
    ) {
      let expected = SubstitutionModel::new(pi.clone(), w.clone(), 1.0).unwrap();
      let actual = SubstitutionModel::new(pi * 10.0_f64.powf(exponent), w, 1.0).unwrap();
      prop_assert_array_ulps_eq!(expected.pi(), actual.pi(), max_ulps = 14);
      prop_assert_array_ulps_eq!(expected.w(), actual.w(), max_ulps = 46);
    }
  }

  mod generators {
    use ndarray::{Array1, Array2};
    use proptest::collection::vec;
    use proptest::strategy::Strategy;

    pub(super) fn gen_pi() -> impl Strategy<Value = Array1<f64>> {
      vec(0.01..1.0_f64, 4).prop_map(Array1::from_vec)
    }

    pub(super) fn gen_w() -> impl Strategy<Value = Array2<f64>> {
      vec(0.0..1.0_f64, 16).prop_map(|entries| {
        let upper = Array2::from_shape_vec((4, 4), entries).unwrap();
        let mut w = &upper + &upper.t();
        w.diag_mut().fill(0.0);
        w
      })
    }
  }
}
