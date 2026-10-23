#[cfg(test)]
mod tests {
  use crate::__tests__::test_support::tests::helpers::three_tip_tree;
  use crate::observe::{TipObservation, UNKNOWN, codes_containing, observe_tips};
  use helpers::leaves_with;
  use itertools::{Itertools, izip};
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use proptest::prelude::{ProptestConfig, any, prop_assert, proptest};
  use rstest::rstest;
  use std::ops::RangeInclusive;
  use treetime_primitives::{AsciiChar, Seq};
  use treetime_utils::assert_error;
  use treetime_utils::sync::random::get_random_number_generator;

  #[test]
  fn test_observe_tips_without_coarsening_returns_the_leaf_sequences() {
    let tree = three_tip_tree();
    let sequences = btreemap! {
      tree.root  => Seq::try_from_str("ACGT-ACGT").unwrap(),
      tree.inner => Seq::try_from_str("ACGT-ACGT").unwrap(),
      tree.tip_c => Seq::try_from_str("ACGT-ACGA").unwrap(),
      tree.tip_a => Seq::try_from_str("AC---ACGT").unwrap(),
      tree.tip_b => Seq::try_from_str("TCGT-ACGT").unwrap(),
    };
    let observation = TipObservation {
      ambiguity_rate: 0.0,
      unknown_run_rate: 0.0,
      unknown_run_lengths: 1..=1,
    };

    let actual = observe_tips(
      &tree.graph,
      &sequences,
      &observation,
      &mut get_random_number_generator(1),
    )
    .unwrap();

    let expected = btreemap! {
      tree.tip_a => sequences[&tree.tip_a].clone(),
      tree.tip_b => sequences[&tree.tip_b].clone(),
      tree.tip_c => sequences[&tree.tip_c].clone(),
    };
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_observe_tips_with_one_run_per_tip_puts_exactly_that_run_of_n_into_each_tip() {
    let tree = three_tip_tree();
    let observation = TipObservation {
      ambiguity_rate: 0.0,
      unknown_run_rate: 1.0,
      unknown_run_lengths: 4..=4,
    };

    let observed = observe_tips(
      &tree.graph,
      &leaves_with("ACGTACGTACGT"),
      &observation,
      &mut get_random_number_generator(2),
    )
    .unwrap();

    let run_lengths = observed
      .values()
      .map(|seq| {
        seq
          .iter()
          .map(|&state| state == UNKNOWN)
          .dedup_with_count()
          .filter(|&(_, unknown)| unknown)
          .map(|(count, _)| count)
          .collect_vec()
      })
      .collect_vec();
    assert_eq!(vec![vec![4]; 3], run_lengths);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::ambiguity_rate_negative(TipObservation { ambiguity_rate: -0.1, unknown_run_rate: 0.0, unknown_run_lengths: 1..=1 }, "Ambiguity rate must be in [0, 1], got -0.1")]
  #[case::run_rate_above_one(     TipObservation { ambiguity_rate:  0.0, unknown_run_rate: 1.5, unknown_run_lengths: 1..=1 }, "Rate of runs of N must be in [0, 1], got 1.5")]
  #[case::run_length_zero(        TipObservation { ambiguity_rate:  0.0, unknown_run_rate: 0.5, unknown_run_lengths: 0..=2 }, "Lengths of runs of N must be a non-empty range of positive lengths, got 0..=2")]
  #[case::run_lengths_empty(      TipObservation { ambiguity_rate:  0.0, unknown_run_rate: 0.5, unknown_run_lengths: RangeInclusive::new(3, 2) }, "Lengths of runs of N must be a non-empty range of positive lengths, got 3..=2")]
  #[case::run_longer_than_tip(    TipObservation { ambiguity_rate:  0.0, unknown_run_rate: 1.0, unknown_run_lengths: 5..=9 }, "Runs of N up to 9 sites do not fit into a sequence of 8 sites")]
  #[trace]
  fn test_observe_tips_rejects_invalid_settings(#[case] observation: TipObservation, #[case] expected: &str) {
    let tree = three_tip_tree();
    let sequences = leaves_with("ACGTACGT");

    assert_error!(observe_tips(&tree.graph, &sequences, &observation, &mut get_random_number_generator(3)), expected);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::a(b'A', "DHMRVW")]
  #[case::c(b'C', "BHMSVY")]
  #[case::g(b'G', "BDKRSV")]
  #[case::t(b'T', "BDHKWY")]
  #[trace]
  fn test_observe_codes_containing_lists_every_iupac_code_with_the_base(#[case] base: u8, #[case] expected: &str) {
    let actual = codes_containing(AsciiChar::from_byte_unchecked(base))
      .unwrap()
      .iter()
      .map(|code| char::from(code.inner()))
      .sorted()
      .collect::<String>();
    assert_eq!(expected, actual);
  }

  proptest! {
    #![proptest_config(ProptestConfig::with_cases(64))]

    #[test]
    fn test_prop_observe_tips_keeps_gaps_and_codes_only_states_that_contain_the_true_base(
      seed in any::<u64>(),
      ambiguity_rate in 0.0..=1.0_f64,
    ) {
      let tree = three_tip_tree();
      let sequences = leaves_with("ACGT--ACGTACGTTTGCA-");
      let observation = TipObservation { ambiguity_rate, unknown_run_rate: 0.0, unknown_run_lengths: 1..=1 };

      let observed = observe_tips(&tree.graph, &sequences, &observation, &mut get_random_number_generator(seed)).unwrap();

      let unexplained = observed
        .iter()
        .flat_map(|(node, seq)| izip!(sequences[node].iter().copied(), seq.iter().copied()))
        .filter(|&(truth, state)| state != truth && !codes_containing(truth).is_some_and(|codes| codes.contains(&state)))
        .collect_vec();
      prop_assert!(unexplained.is_empty(), "observed states that do not contain the true state: {unexplained:?}");
    }
  }

  mod helpers {
    use crate::__tests__::test_support::tests::helpers::three_tip_tree;
    use std::collections::BTreeMap;
    use treetime_graph::node::GraphNodeKey;
    use treetime_primitives::Seq;

    pub(super) fn leaves_with(seq: &str) -> BTreeMap<GraphNodeKey, Seq> {
      let tree = three_tip_tree();
      [tree.tip_a, tree.tip_b, tree.tip_c]
        .into_iter()
        .map(|leaf| (leaf, Seq::try_from_str(seq).unwrap()))
        .collect()
    }
  }
}
