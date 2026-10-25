#[cfg(test)]
mod tests {
  use helpers::{SAMPLING_LAYOUTS, TIME_SCALES, normalized_factor_sum};
  use proptest::arbitrary::any;
  use proptest::prelude::{ProptestConfig, TestCaseError, proptest};
  use proptest::sample::select;
  use treetime_utils::prop_assert_abs_diff_eq;

  proptest! {
    #![proptest_config(ProptestConfig::with_cases(48))]

    #[test]
    fn test_prop_mean_field_normalized_point_prior_matches_the_kingman_density(
      seed in any::<u64>(),
      layout in select(SAMPLING_LAYOUTS.to_vec()),
      tc in select(TIME_SCALES.to_vec()),
    ) {
      let sums = normalized_factor_sum(seed, layout, tc).map_err(|error| TestCaseError::fail(format!("{error:?}")))?;

      prop_assert_abs_diff_eq!(sums.kingman, sums.normalized_factors, epsilon = 1e-9);
    }
  }

  mod helpers {
    use crate::test_utils::{kingman_neg_log_density, point_coalescent_prior_by_key, simulated_node_ages};
    use eyre::Report;
    use treetime_graph::key_map::NodeMap;

    use treetime_simulation::coalescent::{SampleBatch, simulate_coalescent};
    use treetime_simulation::tc_schedule::TcSchedule;
    use treetime_utils::sync::random::get_random_number_generator;

    pub(super) const SAMPLING_LAYOUTS: [&[(f64, usize)]; 4] = [
      &[(0.0, 2)],
      &[(0.0, 7)],
      &[(0.0, 3), (0.4, 2), (1.5, 4)],
      &[(0.0, 2), (1.0, 2), (2.0, 2), (3.0, 2), (4.0, 2), (5.0, 2)],
    ];
    pub(super) const TIME_SCALES: [f64; 3] = [0.5, 3.0, 20.0];

    pub(super) struct FactorSums {
      pub kingman: f64,
      pub normalized_factors: f64,
    }

    pub(super) fn normalized_factor_sum(seed: u64, layout: &[(f64, usize)], tc: f64) -> Result<FactorSums, Report> {
      let sampling = layout
        .iter()
        .map(|&(tau, count)| SampleBatch { tau, count })
        .collect::<Vec<_>>();
      let tree = simulate_coalescent(
        &TcSchedule::constant(tc)?,
        &sampling,
        &mut get_random_number_generator(seed),
      )?;
      let dates: NodeMap<_> = tree
        .graph
        .get_nodes()
        .map(|node| (node.key(), tree.date(node.key())))
        .collect();
      let prior = point_coalescent_prior_by_key(&tree.graph, &dates, tc)?;
      let factors = tree
        .graph
        .get_nodes()
        .map(|node| prior.node_contribution(node.key(), node.outbound().len(), dates[&node.key()]))
        .sum::<Result<f64, Report>>()?;
      Ok(FactorSums {
        kingman: kingman_neg_log_density(&simulated_node_ages(&tree), tc)?,
        normalized_factors: factors - prior.log_normalization(),
      })
    }
  }
}
