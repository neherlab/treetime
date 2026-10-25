#[cfg(test)]
mod tests {
  use self::helpers::acgt_sequences;
  use crate::__tests__::test_support::tests::helpers::three_tip_tree;
  use crate::gaps::{GapEvent, apply_gap_events, disjoint_gap_events};
  use itertools::Itertools;
  use treetime_graph::keymap;

  use pretty_assertions::assert_eq;
  use proptest::prelude::{ProptestConfig, any, prop_assert, prop_assert_eq, proptest};
  use rstest::rstest;
  use treetime_primitives::Seq;
  use treetime_utils::assert_error;
  use treetime_utils::cast::{ToF64Exact, f64_to_usize_saturating};
  use treetime_utils::sync::random::get_random_number_generator;

  #[test]
  fn test_gaps_apply_deletes_the_sites_in_the_child_and_every_node_below() {
    let tree = three_tip_tree();
    let sequences = keymap! {
      tree.root  => Seq::try_from_str("ACGTACGTAC").unwrap(),
      tree.inner => Seq::try_from_str("ACGTACGTAC").unwrap(),
      tree.tip_c => Seq::try_from_str("ACGTACGTAC").unwrap(),
      tree.tip_a => Seq::try_from_str("ACGTACGTAC").unwrap(),
      tree.tip_b => Seq::try_from_str("ACGTACGTAC").unwrap(),
    };
    let events = [
      GapEvent {
        edge: tree.root_inner,
        start: 2,
        length: 3,
      },
      GapEvent {
        edge: tree.inner_a,
        start: 6,
        length: 2,
      },
    ];

    let actual = apply_gap_events(&tree.graph, &sequences, &events).unwrap();

    let expected = keymap! {
      tree.root  => Seq::try_from_str("ACGTACGTAC").unwrap(),
      tree.inner => Seq::try_from_str("AC---CGTAC").unwrap(),
      tree.tip_c => Seq::try_from_str("ACGTACGTAC").unwrap(),
      tree.tip_a => Seq::try_from_str("AC---C--AC").unwrap(),
      tree.tip_b => Seq::try_from_str("AC---CGTAC").unwrap(),
    };
    assert_eq!(expected, actual);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::past_the_end((3, 2), "Gap event on edge 2 covers sites 3..5, outside the 4 sites of the sequence")]
  #[case::empty(       (1, 0), "Gap event on edge 2 covers sites 1..1, outside the 4 sites of the sequence")]
  #[trace]
  fn test_gaps_apply_rejects_an_event_outside_the_sequence(#[case] (start, length): (usize, usize), #[case] expected: &str) {
    let tree = three_tip_tree();
    let event = GapEvent { edge: tree.inner_a, start, length };

    assert_error!(apply_gap_events(&tree.graph, &acgt_sequences(&tree), &[event]), expected);
  }

  #[test]
  fn test_gaps_apply_rejects_an_event_on_an_edge_of_another_tree() {
    let tree = three_tip_tree();
    let mut other = three_tip_tree();
    let extra = other.graph.add_node();
    let foreign = other.graph.add_edge(other.tip_c, extra).unwrap();
    let event = GapEvent {
      edge: foreign,
      start: 0,
      length: 1,
    };

    assert_error!(
      apply_gap_events(&tree.graph, &acgt_sequences(&tree), &[event]),
      "Gap event at site 0 is on edge 4, which the tree has not"
    );
  }

  #[test]
  fn test_gaps_disjoint_events_reject_more_gaps_than_the_sequence_holds() {
    let tree = three_tip_tree();
    let durations = keymap! { tree.root_inner => 1.0, tree.root_c => 1.0, tree.inner_a => 1.0, tree.inner_b => 1.0 };

    assert_error!(
      disjoint_gap_events(&durations, 3, 3, 10, &mut get_random_number_generator(1)),
      "A sequence of 10 sites holds at most 2 gaps of 3 sites with a site between them, got 3 gaps"
    );
  }

  #[test]
  fn test_gaps_disjoint_events_reject_edges_without_duration() {
    let tree = three_tip_tree();
    let durations = keymap! { tree.root_inner => 0.0, tree.root_c => 0.0, tree.inner_a => 0.0, tree.inner_b => 0.0 };

    assert_error!(
      disjoint_gap_events(&durations, 1, 3, 10, &mut get_random_number_generator(1)),
      "Gap events need an edge with a positive duration: Not enough weights > zero"
    );
  }

  proptest! {
    #![proptest_config(ProptestConfig::with_cases(64))]

    #[test]
    fn test_prop_gaps_disjoint_events_leave_a_site_between_gaps_on_edges_with_duration(
      seed in any::<u64>(),
      gap_length in 1..=5_usize,
      sequence_length in 10..=200_usize,
      count_fraction in 0.0..=1.0_f64,
    ) {
      let tree = three_tip_tree();
      let durations = keymap! { tree.root_inner => 0.5, tree.root_c => 0.0, tree.inner_a => 2.0, tree.inner_b => 1.0 };
      #[expect(clippy::integer_division, reason = "the slot count is the floor of the sites over the slot width")]
      let max_count = (sequence_length + 1) / (gap_length + 1);
      let count = f64_to_usize_saturating(count_fraction * max_count.to_f64_exact());

      let events = disjoint_gap_events(&durations, count, gap_length, sequence_length, &mut get_random_number_generator(seed)).unwrap();

      prop_assert_eq!(count, events.len());
      prop_assert!(events.iter().all(|event| event.length == gap_length && event.end() <= sequence_length), "{events:?}");
      prop_assert!(events.iter().tuple_windows().all(|(left, right)| left.end() < right.start), "{events:?}");
      prop_assert!(events.iter().all(|event| event.edge != tree.root_c), "{events:?}");
    }
  }

  mod helpers {
    use crate::__tests__::test_support::tests::helpers::ThreeTipTree;
    use treetime_graph::key_map::NodeMap;

    use treetime_primitives::Seq;

    pub(super) fn acgt_sequences(tree: &ThreeTipTree) -> NodeMap<Seq> {
      [tree.root, tree.inner, tree.tip_a, tree.tip_b, tree.tip_c]
        .into_iter()
        .map(|node| (node, Seq::try_from_str("ACGT").unwrap()))
        .collect()
    }
  }
}
