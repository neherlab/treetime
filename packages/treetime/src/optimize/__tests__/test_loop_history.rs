#[cfg(test)]
mod tests {
  use crate::optimize::loop_history::{BestIteration, ConvergenceReason, TopologyHistory};
  use helpers::record_all;
  use pretty_assertions::assert_eq;
  use rstest::rstest;

  #[rustfmt::skip]
  #[rstest]
  #[case::converged_across_a_change(  &[-100.0, -90.0], &[-90.0005],         None                              )]
  #[case::oscillating_across_a_change(&[-100.0, -90.0], &[-100.0005],        None                              )]
  #[case::worsened_across_a_change(   &[-100.0, -90.0], &[-95.0],            None                              )]
  #[case::converged_within_topology(  &[],              &[-90.0, -90.0005],  Some(ConvergenceReason::Converged))]
  #[trace]
  fn test_loop_history_judges_the_stop_only_on_likelihoods_of_the_current_topology(
    #[case] before_change: &[f64],
    #[case] after_change: &[f64],
    #[case] expected: Option<ConvergenceReason>,
  ) {
    let mut history = TopologyHistory::new(1e-3);
    record_all(&mut history, 0, before_change);
    history.topology_changed();

    let actual = record_all(&mut history, before_change.len(), after_change);

    assert_eq!(expected, actual);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::a_later_iteration_of_the_new_topology(&[-50.0, -60.0], Some(3))]
  #[case::only_the_iteration_after_the_change(  &[-50.0],        None   )]
  #[trace]
  fn test_loop_history_never_keeps_the_iteration_after_a_topology_change_as_best(
    #[case] after_change: &[f64],
    #[case] expected: Option<usize>,
  ) {
    let mut history = TopologyHistory::new(1e-3);
    record_all(&mut history, 0, &[-100.0, -90.0]);
    history.topology_changed();

    record_all(&mut history, 2, after_change);

    assert_eq!(expected, history.into_best().as_ref().map(BestIteration::iteration));
  }

  mod helpers {
    use crate::optimize::indel::EdgeIndels;
    use crate::optimize::loop_history::{ConvergenceReason, TopologyHistory};
    use crate::partition::set::PartitionSet;
    use treetime_graph::key_map::EdgeMap;
    use treetime_primitives::LogLh;

    pub(super) fn record_all(
      history: &mut TopologyHistory,
      first_iteration: usize,
      log_lhs: &[f64],
    ) -> Option<ConvergenceReason> {
      let branch_lengths = EdgeMap::new();
      let partitions = PartitionSet::Sparse(vec![]);
      log_lhs.iter().enumerate().fold(None, |_, (i, &log_lh)| {
        let stop = history.record(
          first_iteration + i,
          LogLh::new(log_lh),
          &branch_lengths,
          &partitions,
          &EdgeIndels::Excluded,
        );
        history.lengths_optimized();
        stop
      })
    }
  }
}
