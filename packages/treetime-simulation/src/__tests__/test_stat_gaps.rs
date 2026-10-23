#[cfg(test)]
mod tests {
  use crate::__tests__::test_support::tests::helpers::{chisquare_p_value, three_tip_tree};
  use crate::FALSE_FAIL_RATE;
  use crate::gaps::disjoint_gap_events;
  use itertools::Itertools;
  use maplit::btreemap;
  use std::iter::repeat_with;
  use treetime_utils::sync::random::get_random_number_generator;

  const ROUNDS: usize = 4_000;
  const EVENTS_PER_ROUND: usize = 5;

  #[test]
  fn test_stat_gaps_disjoint_events_choose_edges_in_proportion_to_their_duration() {
    let tree = three_tip_tree();
    let durations = btreemap! { tree.root_inner => 0.5, tree.root_c => 1.0, tree.inner_a => 2.0, tree.inner_b => 0.5 };
    let mut rng = get_random_number_generator(51);

    let edges = repeat_with(|| disjoint_gap_events(&durations, EVENTS_PER_ROUND, 3, 100, &mut rng).unwrap())
      .take(ROUNDS)
      .flatten()
      .map(|event| event.edge)
      .counts();

    let observed = durations
      .keys()
      .map(|edge| edges.get(edge).copied().unwrap_or(0))
      .collect_vec();
    let probabilities = durations.values().map(|duration| duration / 4.0).collect_vec();
    let p_value = chisquare_p_value(&observed, &probabilities);
    assert!(
      p_value > FALSE_FAIL_RATE,
      "chi-square test of the edges of {ROUNDS} rounds of {EVENTS_PER_ROUND} gap events against the edge durations: counts {observed:?}, p = {p_value}"
    );
  }
}
