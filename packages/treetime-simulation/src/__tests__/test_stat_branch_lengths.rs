#[cfg(test)]
mod tests {
  use crate::FALSE_FAIL_RATE;
  use crate::branch_lengths::perturbed;
  use itertools::Itertools;
  use statrs::distribution::Normal;
  use statrs::stats_tests::NaNPolicy;
  use statrs::stats_tests::ks_test::{KSOneSampleAlternativeMethod, ks_onesample};
  use treetime_graph::key_map::EdgeMap;

  use std::iter::repeat_with;
  use treetime_graph::graph::Graph;
  use treetime_utils::sync::random::get_random_number_generator;

  const EDGE_COUNT: usize = 10_000;

  #[test]
  fn test_stat_branch_lengths_perturbed_log_ratios_are_normal() {
    let sigma = 0.3;
    let mut graph = Graph::new();
    let mut parent = graph.add_node();
    let lengths: EdgeMap<f64> = repeat_with(|| {
      let child = graph.add_node();
      let edge = graph.add_edge(parent, child).unwrap();
      parent = child;
      (edge, 0.02)
    })
    .take(EDGE_COUNT)
    .collect();
    let perturbed_lengths = perturbed(&lengths, sigma, &mut get_random_number_generator(40)).unwrap();

    let log_ratios = perturbed_lengths
      .iter()
      .map(|(edge, length)| (length / lengths[edge]).ln())
      .collect_vec();
    let (_, p_value) = ks_onesample(
      log_ratios,
      &Normal::new(0.0, sigma).unwrap(),
      KSOneSampleAlternativeMethod::TwoSidedAsymptotic,
      NaNPolicy::Error,
    )
    .unwrap();
    assert!(
      p_value > FALSE_FAIL_RATE,
      "KS test of {EDGE_COUNT} log ratios against Normal(0, {sigma}): p = {p_value}"
    );
  }
}
