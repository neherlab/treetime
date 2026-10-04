#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::Alphabet;
  use crate::partition::storage::sparse::{FitchNodeData, SparseNodeObs};
  use crate::seq::composition::Composition;
  use crate::seq::find_char_ranges::find_letter_ranges;
  use proptest::prelude::*;
  use std::collections::BTreeMap;
  use treetime_primitives::{AlphabetLike, Seq, StateSet};
  use treetime_utils::interval::range_union::range_union;

  proptest! {
    #[test]
    fn test_prop_sparse_obs_matches_separate_scans(seq in generators::seq()) {
      let alphabet = Alphabet::default();
      let obs = SparseNodeObs::new(&seq, &alphabet);
      let fitch = FitchNodeData::new(seq.clone(), &alphabet).seq;

      let expected = helpers::separate_scans(&seq, &alphabet);
      let actual_obs = (obs.unknown, obs.gaps, obs.non_char, obs.composition, obs.fitch.variable);
      let actual_fitch = (fitch.unknown, fitch.gaps, fitch.non_char, fitch.composition, fitch.fitch.variable);

      prop_assert_eq!(&expected, &actual_obs);
      prop_assert_eq!(&expected, &actual_fitch);
      prop_assert_eq!(&seq, &fitch.sequence);
    }
  }

  mod generators {
    use proptest::collection::vec;
    use proptest::prelude::*;
    use treetime_primitives::{AsciiChar, Seq};

    pub(super) fn seq() -> impl Strategy<Value = Seq> {
      vec(prop::sample::select(b"ACGTNNNN----RYKM".to_vec()), 0..200)
        .prop_map(|bytes| bytes.into_iter().map(AsciiChar::from_byte_unchecked).collect())
    }
  }

  mod helpers {
    use super::*;

    type Observation = (
      Vec<(usize, usize)>,
      Vec<(usize, usize)>,
      Vec<(usize, usize)>,
      Composition,
      BTreeMap<usize, StateSet>,
    );

    pub(super) fn separate_scans(seq: &Seq, alphabet: &Alphabet) -> Observation {
      let unknown = find_letter_ranges(seq, alphabet.unknown());
      let gaps = find_letter_ranges(seq, alphabet.gap());
      let non_char = range_union(&[unknown.clone(), gaps.clone()]);
      let composition = Composition::with_seq(seq, alphabet.chars(), alphabet.gap());
      let variable = seq
        .iter()
        .enumerate()
        .filter(|&(_, c)| alphabet.is_ambiguous(*c))
        .map(|(pos, &c)| (pos, alphabet.char_to_set(c)))
        .collect();
      (unknown, gaps, non_char, composition, variable)
    }
  }
}
