#[cfg(test)]
mod tests {
  use crate::__tests__::test_support::tests::helpers::{chisquare_p_value, three_tip_tree};
  use crate::FALSE_FAIL_RATE;
  use crate::evolve::{evolve_sequences, evolve_site};
  use crate::substitution_model::{NUCLEOTIDES, SubstitutionModel};
  use helpers::{f81_transition, state_counts};
  use itertools::{Itertools, iproduct, izip};
  use treetime_graph::keymap;

  use ndarray::{Array1, array};
  use rstest::rstest;
  use statrs::distribution::{Discrete, Poisson};
  use statrs::stats_tests::ttest_onesample::ttest_onesample;
  use statrs::stats_tests::{Alternative, NaNPolicy};
  use std::iter::{once, repeat_with};
  use treetime_utils::cast::{ToF64Exact, usize_to_u64};
  use treetime_utils::o;
  use treetime_utils::sync::random::get_random_number_generator;

  const DRAW_COUNT: usize = 100_000;
  const SITE_COUNT: usize = 200_000;

  #[rustfmt::skip]
  #[rstest]
  #[case::jc69_from_a(     array![0.25, 0.25, 0.25, 0.25], 0, 0.3)]
  #[case::f81_from_c(      array![0.1, 0.2, 0.3, 0.4],     1, 0.5)]
  #[case::f81_from_t_long( array![0.1, 0.2, 0.3, 0.4],     3, 2.0)]
  #[trace]
  fn test_stat_evolve_site_end_states_match_f81_transition_probabilities(
    #[case] pi: Array1<f64>,
    #[case] start: usize,
    #[case] duration: f64,
  ) {
    let model = SubstitutionModel::f81(pi.clone(), 1.0).unwrap();
    let mut rng = get_random_number_generator(30);
    let end_states = repeat_with(|| evolve_site(start, duration, &model, &mut rng).unwrap().state)
      .take(DRAW_COUNT)
      .collect_vec();

    let probabilities = (0..4).map(|end| f81_transition(&pi, start, end, duration)).collect_vec();
    let p_value = chisquare_p_value(&state_counts(&end_states), &probabilities);
    assert!(
      p_value > FALSE_FAIL_RATE,
      "chi-square test of {DRAW_COUNT} end states against the F81 transition probabilities: p = {p_value}"
    );
  }

  #[test]
  fn test_stat_evolve_site_jc69_jump_counts_are_poisson() {
    let (mu, duration) = (1.0, 2.0);
    let model = SubstitutionModel::jc69(mu).unwrap();
    let mut rng = get_random_number_generator(31);
    let jump_counts = repeat_with(|| evolve_site(0, duration, &model, &mut rng).unwrap().jumps)
      .take(DRAW_COUNT)
      .counts();

    let bin_count = jump_counts.keys().max().unwrap() + 2;
    let poisson = Poisson::new(mu * duration).unwrap();
    let mut probabilities = (0..bin_count - 1)
      .map(|jumps| poisson.pmf(usize_to_u64(jumps)))
      .collect_vec();
    probabilities.push(1.0 - probabilities.iter().sum::<f64>());
    let observed = (0..bin_count)
      .map(|jumps| jump_counts.get(&jumps).copied().unwrap_or(0))
      .collect_vec();
    let p_value = chisquare_p_value(&observed, &probabilities);
    assert!(
      p_value > FALSE_FAIL_RATE,
      "chi-square test of {DRAW_COUNT} JC69 jump counts against Poisson({}): p = {p_value}",
      mu * duration
    );
  }

  #[test]
  fn test_stat_evolve_site_gtr_mean_jump_count_from_pi_is_rate_times_duration() {
    let mut rng = get_random_number_generator(32);
    let (mu, duration) = (1.5, 0.8);
    let model = SubstitutionModel::random(mu, &mut rng).unwrap();
    let jump_counts = repeat_with(|| {
      let start = model.sample_state(&mut rng);
      evolve_site(start, duration, &model, &mut rng)
        .unwrap()
        .jumps
        .to_f64_exact()
    })
    .take(DRAW_COUNT)
    .collect_vec();

    let (_, p_value) = ttest_onesample(jump_counts, mu * duration, Alternative::TwoSided, NaNPolicy::Error).unwrap();
    assert!(
      p_value > FALSE_FAIL_RATE,
      "t-test of the mean of {DRAW_COUNT} GTR jump counts against {}: p = {p_value}",
      mu * duration
    );
  }

  #[test]
  fn test_stat_evolve_site_gtr_keeps_pi_stationary() {
    let mut rng = get_random_number_generator(33);
    let model = SubstitutionModel::random(1.0, &mut rng).unwrap();
    let end_states = repeat_with(|| {
      let start = model.sample_state(&mut rng);
      evolve_site(start, 0.7, &model, &mut rng).unwrap().state
    })
    .take(DRAW_COUNT)
    .collect_vec();

    let p_value = chisquare_p_value(&state_counts(&end_states), &model.pi().to_vec());
    assert!(
      p_value > FALSE_FAIL_RATE,
      "chi-square test of {DRAW_COUNT} GTR end states against pi {}: p = {p_value}",
      model.pi()
    );
  }

  #[test]
  fn test_stat_evolve_sequences_tip_patterns_and_edge_transitions_match_jc69() {
    let tree = three_tip_tree();
    let durations = keymap! { tree.root_inner => 0.2, tree.root_c => 0.6, tree.inner_a => 0.1, tree.inner_b => 0.9 };
    let model = SubstitutionModel::jc69(1.0).unwrap();
    let simulated = evolve_sequences(
      &tree.graph,
      &durations,
      &model,
      SITE_COUNT,
      &mut get_random_number_generator(34),
    )
    .unwrap();

    let state_of = |node| {
      simulated.sequences[&node]
        .iter()
        .map(|character| {
          NUCLEOTIDES
            .iter()
            .position(|nucleotide| nucleotide == character)
            .unwrap()
        })
        .collect_vec()
    };
    let pi = model.pi().clone();
    let transition = |from, to, duration| f81_transition(&pi, from, to, duration);

    let pattern_counts = izip!(state_of(tree.tip_a), state_of(tree.tip_b), state_of(tree.tip_c)).counts();
    let observed = iproduct!(0..4, 0..4, 0..4)
      .map(|pattern| pattern_counts.get(&pattern).copied().unwrap_or(0))
      .collect_vec();
    let probabilities = iproduct!(0..4, 0..4, 0..4)
      .map(|(a, b, c)| {
        iproduct!(0..4, 0..4)
          .map(|(root, inner)| {
            pi[root]
              * transition(root, inner, durations[&tree.root_inner])
              * transition(inner, a, durations[&tree.inner_a])
              * transition(inner, b, durations[&tree.inner_b])
              * transition(root, c, durations[&tree.root_c])
          })
          .sum::<f64>()
      })
      .collect_vec();
    let tip_check = (o!("tip patterns"), chisquare_p_value(&observed, &probabilities));

    let edge_checks = tree.graph.get_edges().map(|edge| {
      let duration = durations[&edge.key()];
      let pair_counts = izip!(state_of(edge.source()), state_of(edge.target())).counts();
      let observed = iproduct!(0..4, 0..4)
        .map(|pair| pair_counts.get(&pair).copied().unwrap_or(0))
        .collect_vec();
      let probabilities = iproduct!(0..4, 0..4)
        .map(|(parent, child)| pi[parent] * transition(parent, child, duration))
        .collect_vec();
      (
        format!("state pairs of edge {} of duration {duration}", edge.key()),
        chisquare_p_value(&observed, &probabilities),
      )
    });

    let p_values = once(tip_check).chain(edge_checks).collect_vec();
    let threshold = FALSE_FAIL_RATE / p_values.len().to_f64_exact();
    let failed = p_values
      .iter()
      .filter(|(_, p_value)| *p_value <= threshold)
      .collect_vec();
    assert!(
      failed.is_empty(),
      "chi-square tests over {SITE_COUNT} sites against the JC69 probabilities with p <= {threshold}: {failed:?}"
    );
  }

  mod helpers {
    use itertools::Itertools;
    use ndarray::Array1;

    pub(super) fn f81_transition(pi: &Array1<f64>, from: usize, to: usize, duration: f64) -> f64 {
      let beta = 1.0 / (1.0 - pi.dot(pi));
      let decay = (-beta * duration).exp();
      let stay = if from == to { decay } else { 0.0 };
      stay + (1.0 - decay) * pi[to]
    }

    pub(super) fn state_counts(states: &[usize]) -> Vec<usize> {
      let counts = states.iter().counts();
      (0..4)
        .map(|state| counts.get(&state).copied().unwrap_or(0))
        .collect_vec()
    }
  }
}
