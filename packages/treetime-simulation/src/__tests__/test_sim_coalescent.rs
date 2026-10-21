#[cfg(test)]
mod tests {
  use crate::FALSE_FAIL_RATE;
  use crate::coalescent::{SampleBatch, simulate_coalescent};
  use crate::tc_schedule::TcSchedule;
  use eyre::Report;
  use helpers::{first_merger_pair, merger_hazards};
  use itertools::Itertools;
  use rstest::rstest;
  use statrs::distribution::Uniform;
  use statrs::stats_tests::NaNPolicy;
  use statrs::stats_tests::chisquare::chisquare;
  use statrs::stats_tests::ks_test::{KSOneSampleAlternativeMethod, ks_onesample};
  use std::iter::repeat_with;
  use treetime_utils::sync::random::get_random_number_generator;

  const TREE_COUNT: usize = 1000;
  const PAIR_TREE_COUNT: usize = 12_000;

  #[rustfmt::skip]
  #[rstest]
  #[case::constant(   TcSchedule::constant(10.0))]
  #[case::exponential(TcSchedule::exponential(10.0, 5.0))]
  #[case::piecewise(  TcSchedule::piecewise(vec![5.0, 6.0], vec![10.0, 1.0, 3.0]))]
  #[trace]
  fn test_sim_coalescent_merger_hazards_are_uniform_after_time_change(#[case] schedule: Result<TcSchedule, Report>) {
    let schedule = schedule.unwrap();
    let sampling = [
      SampleBatch { tau: 2.0, count: 4 },
      SampleBatch { tau: 0.0, count: 5 },
      SampleBatch { tau: 6.0, count: 3 },
    ];
    let mut rng = get_random_number_generator(20);
    let uniforms = (0..TREE_COUNT)
      .flat_map(|_| {
        let tree = simulate_coalescent(&schedule, &sampling, &mut rng).unwrap();
        merger_hazards(&tree, &schedule, &sampling)
      })
      .map(|hazard| -(-hazard).exp_m1())
      .collect_vec();

    let sample_size = uniforms.len();
    let (_, p_value) = ks_onesample(
      uniforms,
      &Uniform::new(0.0, 1.0).unwrap(),
      KSOneSampleAlternativeMethod::TwoSidedAsymptotic,
      NaNPolicy::Error,
    )
    .unwrap();
    assert!(
      p_value > FALSE_FAIL_RATE,
      "KS test of {sample_size} time-changed merger waiting times: p = {p_value}"
    );
  }

  #[test]
  fn test_sim_coalescent_first_merger_joins_each_pair_equally_often() {
    let schedule = TcSchedule::constant(1.0).unwrap();
    let sampling = [SampleBatch { tau: 0.0, count: 4 }];
    let mut rng = get_random_number_generator(21);
    let pair_counts = repeat_with(|| first_merger_pair(&simulate_coalescent(&schedule, &sampling, &mut rng).unwrap()))
      .take(PAIR_TREE_COUNT)
      .counts();
    let observed = (0..4)
      .tuple_combinations()
      .map(|pair: (usize, usize)| pair_counts.get(&pair).copied().unwrap_or(0))
      .collect_vec();

    let (_, p_value) = chisquare(&observed, None, None).unwrap();
    assert!(
      p_value > FALSE_FAIL_RATE,
      "chi-square test of first merger pairs {observed:?} over {PAIR_TREE_COUNT} trees: p = {p_value}"
    );
  }

  mod helpers {
    use crate::coalescent::{SampleBatch, SimulatedTree};
    use crate::tc_schedule::TcSchedule;
    use itertools::Itertools;
    use std::iter::once;
    use treetime_utils::cast::ToF64Exact;

    pub(super) fn merger_hazards(tree: &SimulatedTree, schedule: &TcSchedule, sampling: &[SampleBatch]) -> Vec<f64> {
      let merger_times = tree
        .graph
        .get_internal_nodes()
        .map(|node| tree.tau[&node.key()])
        .sorted_by(f64::total_cmp)
        .collect_vec();
      let batch_times = sampling
        .iter()
        .map(|batch| batch.tau)
        .sorted_by(f64::total_cmp)
        .collect_vec();
      let arrived_by = |time: f64| -> usize {
        sampling
          .iter()
          .filter(|batch| batch.tau <= time)
          .map(|batch| batch.count)
          .sum()
      };

      let mut start = batch_times[0];
      merger_times
        .iter()
        .enumerate()
        .map(|(merged_before, &end)| {
          let bounds = once(start)
            .chain(batch_times.iter().copied().filter(|&time| start < time && time < end))
            .chain([end])
            .collect_vec();
          let hazard = bounds
            .iter()
            .tuple_windows()
            .map(|(&from, &to)| {
              let lineages = (arrived_by(from) - merged_before).to_f64_exact();
              lineages * (lineages - 1.0) * 0.5 * (schedule.cumulative_hazard(to) - schedule.cumulative_hazard(from))
            })
            .sum();
          start = end;
          hazard
        })
        .collect()
    }

    pub(super) fn first_merger_pair(tree: &SimulatedTree) -> (usize, usize) {
      let first_merger = tree
        .graph
        .get_internal_nodes()
        .min_by(|left, right| tree.tau[&left.key()].total_cmp(&tree.tau[&right.key()]))
        .unwrap();
      let tips = tree.graph.get_leaves().map(|leaf| leaf.key()).collect_vec();
      tree
        .graph
        .children_of(first_merger)
        .map(|(child, _)| tips.iter().position(|&tip| tip == child.key()).unwrap())
        .sorted()
        .collect_tuple()
        .unwrap()
    }
  }
}
