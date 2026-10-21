#[cfg(test)]
mod tests {
  use crate::substitution_model::SubstitutionModel;
  use proptest::{prop_assert_eq, proptest};
  use treetime_utils::prop_assert_ulps_eq;
  use treetime_utils::sync::random::get_random_number_generator;

  proptest! {
    #[test]
    fn test_prop_substitution_model_random_has_unit_average_rate(seed in 0..u64::MAX, mu in 0.0..10.0_f64) {
      let model = SubstitutionModel::random(mu, &mut get_random_number_generator(seed)).unwrap();
      let pi = model.pi();
      let average_rate = pi.dot(&model.w().dot(pi));
      prop_assert_ulps_eq!(1.0, average_rate, max_ulps = 4);
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
      prop_assert_ulps_eq!(1.0, model.pi().sum(), max_ulps = 2);
    }
  }
}
