#[cfg(test)]
mod tests {
  use crate::partition::marginal::posterior_states::most_likely_state;
  use ndarray::{Array1, array};
  use pretty_assertions::assert_eq;
  use proptest::prelude::*;
  use rstest::rstest;

  #[rustfmt::skip]
  #[rstest]
  #[case::clear_maximum(        array![0.1, 0.6, 0.3],                  Some(1))]
  #[case::exact_tie(            array![0.1, 0.45, 0.45],                Some(1))]
  #[case::rounding_tie(         array![0.3916531565081366, 0.39165315650834764, 0.21669368698351576], Some(0))]
  #[case::difference_above(     array![0.5 - 6e-10, 0.5 + 6e-10],      Some(1))]
  #[case::empty(                Array1::<f64>::zeros(0),                None)]
  #[trace]
  fn test_posterior_states_most_likely_state_takes_first_state_within_tie_tolerance(
    #[case] posterior: Array1<f64>,
    #[case] expected: Option<usize>,
  ) {
    assert_eq!(expected, most_likely_state(&posterior));
  }

  proptest! {
    #[test]
    fn test_prop_posterior_states_most_likely_state_ignores_differences_below_tie_tolerance(
      first in 0.34..0.5_f64,
      offset in 0.0..1e-10_f64,
    ) {
      let posterior = array![first, first + offset, 1.0 - 2.0 * first - offset];
      prop_assert_eq!(Some(0), most_likely_state(&posterior));
    }
  }
}
