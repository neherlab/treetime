#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::partition::marginal::posterior_states::{AmbiguityThreshold, most_likely_state};
  use eyre::Report;
  use itertools::Itertools;
  use ndarray::{Array1, array};
  use pretty_assertions::assert_eq;
  use proptest::prelude::*;
  use rstest::rstest;
  use std::str::FromStr;
  use treetime_primitives::AsciiChar;
  use treetime_utils::{assert_error, pretty_assert_ulps_eq};

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

  #[rustfmt::skip]
  #[rstest]
  #[case::zero(  0.0)]
  #[case::raxml( 0.125)]
  #[case::one(   1.0)]
  #[trace]
  fn test_posterior_states_threshold_accepts_probability_difference(#[case] value: f64) -> Result<(), Report> {
    pretty_assert_ulps_eq!(value, AmbiguityThreshold::new(value)?.get(), max_ulps = 0);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::negative( "-0.1", "the ambiguity threshold must be a probability difference from 0 to 1, got -0.1")]
  #[case::above_one("1.5",  "the ambiguity threshold must be a probability difference from 0 to 1, got 1.5")]
  #[case::nan(      "NaN",  "the ambiguity threshold must be a probability difference from 0 to 1, got NaN")]
  #[case::text(     "half", "'half' is not a number: invalid float literal")]
  #[trace]
  fn test_posterior_states_threshold_rejects_invalid_text(#[case] text: &str, #[case] expected: &str) {
    assert_error!(AmbiguityThreshold::from_str(text).map_err(Report::new), expected);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::tie_only(          (array![0.45, 0.1, 0.45],     0.0),   vec![0, 2])]
  #[case::rounding_tie(      (array![0.5 - 4e-11, 0.5 + 4e-11], 0.0), vec![0, 1])]
  #[case::within_threshold(  (array![0.5, 0.4, 0.1],       0.125), vec![0, 1])]
  #[case::outside_threshold( (array![0.6, 0.3, 0.1],       0.125), vec![0])]
  #[case::boundary_included( (array![0.5, 0.375, 0.125],   0.125), vec![0, 1])]
  #[trace]
  fn test_posterior_states_plausible_states_are_within_threshold_of_maximum(
    #[case] (posterior, threshold): (Array1<f64>, f64),
    #[case] expected: Vec<usize>,
  ) -> Result<(), Report> {
    let actual = AmbiguityThreshold::new(threshold)?.plausible_states(&posterior).collect_vec();
    assert_eq!(expected, actual);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::nuc_purine(       AlphabetName::Nuc, array![0.5, 0.0, 0.5, 0.0],         Some(b'R'))]
  #[case::nuc_pyrimidine(   AlphabetName::Nuc, array![0.0, 0.5, 0.0, 0.5],         Some(b'Y'))]
  #[case::nuc_strong(       AlphabetName::Nuc, array![0.0, 0.5, 0.5, 0.0],         Some(b'S'))]
  #[case::nuc_weak(         AlphabetName::Nuc, array![0.5, 0.0, 0.0, 0.5],         Some(b'W'))]
  #[case::nuc_keto(         AlphabetName::Nuc, array![0.0, 0.0, 0.5, 0.5],         Some(b'K'))]
  #[case::nuc_amino(        AlphabetName::Nuc, array![0.5, 0.5, 0.0, 0.0],         Some(b'M'))]
  #[case::nuc_not_a(        AlphabetName::Nuc, array![0.0, 1.0, 1.0, 1.0] / 3.0,   Some(b'B'))]
  #[case::nuc_not_c(        AlphabetName::Nuc, array![1.0, 0.0, 1.0, 1.0] / 3.0,   Some(b'D'))]
  #[case::nuc_not_g(        AlphabetName::Nuc, array![1.0, 1.0, 0.0, 1.0] / 3.0,   Some(b'H'))]
  #[case::nuc_not_t(        AlphabetName::Nuc, array![1.0, 1.0, 1.0, 0.0] / 3.0,   Some(b'V'))]
  #[case::nuc_any(          AlphabetName::Nuc, array![0.25, 0.25, 0.25, 0.25],     Some(b'N'))]
  #[case::nuc_resolved(     AlphabetName::Nuc, array![0.0, 0.0, 1.0, 0.0],         None)]
  #[case::aa_asx(           AlphabetName::Aa,  helpers::aa_tie(b"DN"),       Some(b'B'))]
  #[case::aa_glx(           AlphabetName::Aa,  helpers::aa_tie(b"EQ"),       Some(b'Z'))]
  #[case::aa_xle(           AlphabetName::Aa,  helpers::aa_tie(b"IL"),       Some(b'J'))]
  #[case::aa_without_code(  AlphabetName::Aa,  helpers::aa_tie(b"AC"),       Some(b'X'))]
  #[case::aa_three_states(  AlphabetName::Aa,  helpers::aa_tie(b"DNE"),  Some(b'X'))]
  #[trace]
  fn test_posterior_states_ambiguity_char_is_the_code_of_the_tied_states(
    #[case] alphabet: AlphabetName,
    #[case] posterior: Array1<f64>,
    #[case] expected: Option<u8>,
  ) -> Result<(), Report> {
    let alphabet = Alphabet::new(alphabet)?;
    let actual = AmbiguityThreshold::new(0.0)?.ambiguity_char(&posterior, &alphabet);
    assert_eq!(expected.map(AsciiChar::from_byte_unchecked), actual);
    Ok(())
  }

  proptest! {
    #[test]
    fn test_prop_posterior_states_most_likely_state_is_plausible(
      weights in prop::collection::vec(0.0..1.0_f64, 2..8),
      threshold in 0.0..1.0_f64,
    ) {
      let total: f64 = weights.iter().sum();
      prop_assume!(total > 0.0);
      let posterior = Array1::from_iter(weights) / total;
      let threshold = AmbiguityThreshold::new(threshold).unwrap();
      let most_likely = most_likely_state(&posterior).unwrap();
      prop_assert!(threshold.plausible_states(&posterior).contains(&most_likely));
    }

    #[test]
    fn test_prop_posterior_states_most_likely_state_ignores_differences_below_tie_tolerance(
      first in 0.34..0.5_f64,
      offset in 0.0..1e-10_f64,
    ) {
      let posterior = array![first, first + offset, 1.0 - 2.0 * first - offset];
      prop_assert_eq!(Some(0), most_likely_state(&posterior));
    }
  }

  mod helpers {
    use crate::alphabet::alphabet::{Alphabet, AlphabetName};
    use ndarray::Array1;
    use treetime_utils::cast::ToF64Exact;

    pub(super) fn aa_tie(states: &[u8]) -> Array1<f64> {
      let alphabet = Alphabet::new(AlphabetName::Aa).expect("the amino-acid alphabet is valid");
      let mut posterior = Array1::zeros(alphabet.n_canonical());
      for &state in states {
        let index = alphabet.index(state).expect("the state is canonical");
        posterior[index] = 1.0 / states.len().to_f64_exact();
      }
      posterior
    }
  }
}
