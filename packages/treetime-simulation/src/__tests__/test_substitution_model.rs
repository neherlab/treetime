#[cfg(test)]
mod tests {
  use crate::substitution_model::SubstitutionModel;
  use eyre::Report;
  use helpers::{asymmetric, off_diagonal_ones, uniform_pi, with_pair};
  use ndarray::{Array1, Array2, array};
  use rstest::rstest;
  use treetime_utils::sync::random::get_random_number_generator;
  use treetime_utils::{assert_error, pretty_assert_ulps_eq};

  #[test]
  fn test_substitution_model_jc69_scales_exchangeabilities_to_unit_average_rate() {
    let model = SubstitutionModel::jc69(0.5).unwrap();
    let off_diagonal = 4.0 / 3.0;
    let expected_pi = array![0.25, 0.25, 0.25, 0.25];
    let expected_w = Array2::from_shape_fn((4, 4), |(row, column)| if row == column { 0.0 } else { off_diagonal });
    pretty_assert_ulps_eq!(expected_pi, model.pi().clone(), max_ulps = 0);
    pretty_assert_ulps_eq!(expected_w, model.w().clone(), max_ulps = 1);
  }

  #[test]
  fn test_substitution_model_random_pi_sums_to_one_within_rounding_at_a_seed_three_ulps_below_one() {
    let model = SubstitutionModel::random(1.0, &mut get_random_number_generator(15_997_397_837_458_979_626)).unwrap();

    pretty_assert_ulps_eq!(1.0, model.pi().sum(), max_ulps = 7);
  }

  #[test]
  fn test_substitution_model_f81_normalizes_pi_and_scales_by_heterozygosity() {
    let model = SubstitutionModel::f81(array![1.0, 2.0, 3.0, 4.0], 1.0).unwrap();
    let off_diagonal = 1.0 / 0.7;
    let expected_pi = array![0.1, 0.2, 0.3, 0.4];
    let expected_w = Array2::from_shape_fn((4, 4), |(row, column)| if row == column { 0.0 } else { off_diagonal });
    pretty_assert_ulps_eq!(expected_pi, model.pi().clone(), max_ulps = 1);
    pretty_assert_ulps_eq!(expected_w, model.w().clone(), max_ulps = 2);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::pi_wrong_length(     array![0.5, 0.5],                 off_diagonal_ones(),   1.0,           "Equilibrium frequencies (pi) need 4 entries for the states A, C, G, T, got 2")]
  #[case::pi_zero_entry(       array![0.5, 0.0, 0.25, 0.25],     off_diagonal_ones(),   1.0,           "Equilibrium frequencies (pi) must be finite and positive, got [0.5, 0, 0.25, 0.25]")]
  #[case::pi_negative_entry(   array![0.5, -0.1, 0.3, 0.3],      off_diagonal_ones(),   1.0,           "Equilibrium frequencies (pi) must be finite and positive, got [0.5, -0.1, 0.3, 0.3]")]
  #[case::pi_nan_entry(        array![0.5, f64::NAN, 0.3, 0.3],  off_diagonal_ones(),   1.0,           "Equilibrium frequencies (pi) must be finite and positive, got [0.5, NaN, 0.3, 0.3]")]
  #[case::w_wrong_shape(       uniform_pi(),                     Array2::ones((3, 3)),  1.0,           "Exchangeabilities (w) must be a 4x4 matrix, got shape [3, 3]")]
  #[case::w_negative_entry(    uniform_pi(),                     with_pair(-1.0),       1.0,           "Exchangeabilities (w) must be finite and non-negative, got [[0, -1, 1, 1],\n [-1, 0, 1, 1],\n [1, 1, 0, 1],\n [1, 1, 1, 0]]")]
  #[case::w_infinite_entry(    uniform_pi(),                     with_pair(f64::INFINITY), 1.0,        "Exchangeabilities (w) must be finite and non-negative, got [[0, inf, 1, 1],\n [inf, 0, 1, 1],\n [1, 1, 0, 1],\n [1, 1, 1, 0]]")]
  #[case::w_non_zero_diagonal( uniform_pi(),                     Array2::ones((4, 4)),  1.0,           "Exchangeabilities (w) must have a zero diagonal, got [[1, 1, 1, 1],\n [1, 1, 1, 1],\n [1, 1, 1, 1],\n [1, 1, 1, 1]]")]
  #[case::w_asymmetric(        uniform_pi(),                     asymmetric(),          1.0,           "Exchangeabilities (w) must be symmetric, got [[0, 2, 1, 1],\n [1, 0, 1, 1],\n [1, 1, 0, 1],\n [1, 1, 1, 0]]")]
  #[case::w_zero(              uniform_pi(),                     Array2::zeros((4, 4)), 1.0,           "Average rate pi^T w pi of the exchangeabilities (w) must be positive, got 0")]
  #[case::mu_negative(         uniform_pi(),                     off_diagonal_ones(),   -1.0,          "Substitution rate (mu) must be finite and non-negative, got -1")]
  #[case::mu_nan(              uniform_pi(),                     off_diagonal_ones(),   f64::NAN,      "Substitution rate (mu) must be finite and non-negative, got NaN")]
  #[case::mu_infinite(         uniform_pi(),                     off_diagonal_ones(),   f64::INFINITY, "Substitution rate (mu) must be finite and non-negative, got inf")]
  #[trace]
  fn test_substitution_model_rejects_invalid_input(
    #[case] pi: Array1<f64>,
    #[case] w: Array2<f64>,
    #[case] mu: f64,
    #[case] expected: &str,
  ) {
    let result: Result<SubstitutionModel, Report> = SubstitutionModel::new(pi, w, mu);
    assert_error!(result, expected);
  }

  #[test]
  fn test_substitution_model_accepts_zero_mu() {
    let model = SubstitutionModel::jc69(0.0).unwrap();
    let expected = 0.0;
    let actual = model.mu();
    pretty_assert_ulps_eq!(expected, actual, max_ulps = 0);
  }

  mod helpers {
    use ndarray::{Array1, Array2, array};

    pub(super) fn uniform_pi() -> Array1<f64> {
      array![0.25, 0.25, 0.25, 0.25]
    }

    pub(super) fn off_diagonal_ones() -> Array2<f64> {
      with_pair(1.0)
    }

    pub(super) fn with_pair(value: f64) -> Array2<f64> {
      let mut w = Array2::from_shape_fn((4, 4), |(row, column)| if row == column { 0.0 } else { 1.0 });
      w[[0, 1]] = value;
      w[[1, 0]] = value;
      w
    }

    pub(super) fn asymmetric() -> Array2<f64> {
      let mut w = off_diagonal_ones();
      w[[0, 1]] = 2.0;
      w
    }
  }
}
