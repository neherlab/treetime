#[cfg(test)]
mod tests {
  use crate::optimize::branch_length::{is_valid_branch_length, is_valid_branch_length_value};
  use crate::optimize::eval::evaluate_site_contributions;
  use crate::optimize::indel::poisson_indel_log_lh;
  use approx::assert_abs_diff_eq;
  use ndarray::array;
  use pretty_assertions::assert_eq;
  use rstest::rstest;

  #[rustfmt::skip]
  #[rstest]
  #[case::negative(         -0.1,           false)]
  #[case::negative_infinity(f64::NEG_INFINITY, false)]
  #[case::zero(              0.0,            true)]
  #[case::positive(          0.1,            true)]
  #[case::positive_infinity(f64::INFINITY, false)]
  #[case::nan(              f64::NAN,      false)]
  #[trace]
  fn test_branch_length_validation_scalar_domain(#[case] branch_length: f64, #[case] expected: bool) {
    assert_eq!(expected, is_valid_branch_length_value(branch_length));
    assert_eq!(expected, is_valid_branch_length(Some(branch_length)));
  }

  #[test]
  fn test_branch_length_validation_missing_is_invalid() {
    assert!(!is_valid_branch_length(None));
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::negative(         -0.1)]
  #[case::negative_infinity(f64::NEG_INFINITY)]
  #[case::positive_infinity(f64::INFINITY)]
  #[case::nan(              f64::NAN)]
  #[trace]
  fn test_branch_length_validation_substitution_evaluator_rejects_invalid(#[case] branch_length: f64) {
    let coefficients = array![1.0, 0.0];
    let eigvals = array![0.0, -1.0];
    let result = evaluate_site_contributions(
      std::iter::once((1.0, coefficients.view())),
      &eigvals,
      branch_length,
      true,
    );
    let error = result.expect_err("invalid branch length must be rejected");
    assert!(error.to_string().contains("finite and non-negative"));
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::negative_no_indels(  0, -0.1)]
  #[case::negative_with_indels(1, -0.1)]
  #[case::negative_infinity(   1, f64::NEG_INFINITY)]
  #[case::positive_infinity(   1, f64::INFINITY)]
  #[case::nan(                 1, f64::NAN)]
  #[trace]
  fn test_branch_length_validation_poisson_evaluator_rejects_invalid(#[case] k: usize, #[case] branch_length: f64) {
    let result = poisson_indel_log_lh(k, 1.0, branch_length);
    let error = result.expect_err("invalid branch length must be rejected");
    assert!(error.to_string().contains("finite and non-negative"));
  }

  #[test]
  fn test_branch_length_validation_poisson_zero_boundary_depends_on_indel_count() {
    let no_indels = poisson_indel_log_lh(0, 1.0, 0.0).expect("zero is valid without indels");
    assert_abs_diff_eq!(0.0, no_indels.log_lh, epsilon = 1e-15);

    let error = poisson_indel_log_lh(1, 1.0, 0.0).expect_err("zero is invalid with observed indels");
    assert!(error.to_string().contains("positive branch length"));
  }
}
