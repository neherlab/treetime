#[cfg(test)]
mod tests {
  use crate::commands::mugration::run::validate_weights;
  use maplit::btreemap;
  use rstest::rstest;
  use std::path::Path;
  use treetime_utils::assert_error;
  use treetime_utils::o;

  #[test]
  fn test_validate_weights_accepts_positive_finite_weights() {
    let weights = btreemap! { o!("usa") => 0.5, o!("uk") => 1e-9 };
    validate_weights(&weights, Path::new("weights.tsv")).unwrap();
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::zero(     0.0,           "0")]
  #[case::negative( -1.5,          "-1.5")]
  #[case::nan(      f64::NAN,      "NaN")]
  #[case::infinite( f64::INFINITY, "inf")]
  #[trace]
  fn test_validate_weights_names_file_and_state_of_invalid_weight(#[case] weight: f64, #[case] shown: &str) {
    let weights = btreemap! { o!("uk") => 2.0, o!("usa") => weight };
    assert_error!(
      validate_weights(&weights, Path::new("data/weights.tsv")),
      format!("Weights file 'data/weights.tsv': weight of state 'usa' is {shown}, expected a finite positive number")
    );
  }
}
