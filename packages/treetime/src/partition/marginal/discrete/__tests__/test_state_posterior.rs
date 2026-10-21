#[cfg(test)]
mod tests {
  use crate::partition::marginal::discrete::state_posterior::StatePosterior;
  use eyre::Report;
  use ndarray::{Array1, array};
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime_utils::cast::ToF64Exact;
  use treetime_utils::{assert_error, pretty_assert_abs_diff_eq};

  #[rustfmt::skip]
  #[rstest]
  #[case::first_of_two(  array![1.0, 0.0],           0)]
  #[case::last_of_five(  array![0.0, 0.0, 0.0, 0.0, 1.0], 4)]
  #[trace]
  fn test_state_posterior_of_one_hot_profile_has_zero_entropy(
    #[case] profile: Array1<f64>,
    #[case] expected_state: usize,
  ) -> Result<(), Report> {
    let posterior = StatePosterior::new(profile)?;

    assert_eq!(expected_state, posterior.state());
    pretty_assert_abs_diff_eq!(0.0, posterior.entropy(), epsilon = 1e-10);
    Ok(())
  }

  #[rstest]
  #[trace]
  fn test_state_posterior_of_uniform_profile_has_log_n_entropy_and_first_state(
    #[values(2, 3, 5)] n_states: usize,
  ) -> Result<(), Report> {
    let posterior = StatePosterior::new(Array1::from_elem(n_states, 1.0 / n_states.to_f64_exact()))?;

    assert_eq!(0, posterior.state());
    pretty_assert_abs_diff_eq!(n_states.to_f64_exact().ln(), posterior.entropy(), epsilon = 1e-10);
    Ok(())
  }

  #[test]
  fn test_state_posterior_of_tie_takes_first_most_likely_state() -> Result<(), Report> {
    let posterior = StatePosterior::new(array![0.1, 0.45, 0.45])?;

    assert_eq!(1, posterior.state());
    Ok(())
  }

  #[test]
  fn test_state_posterior_of_empty_profile_is_an_error() {
    assert_error!(
      StatePosterior::new(Array1::zeros(0)),
      "A state posterior needs a profile with at least one state. This is an internal error. Please report it to developers."
    );
  }
}
