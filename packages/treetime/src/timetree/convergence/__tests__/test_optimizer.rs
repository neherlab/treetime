#[cfg(test)]
mod tests {
  use crate::progress::NoopProgress;
  use crate::test_utils::empty_time_inference;
  use crate::timetree::branch_model::BranchModel;
  use crate::timetree::convergence::metrics::{IterationClock, NODE_TIME_TOLERANCE_YEARS};
  use crate::timetree::convergence::node_times::NodeTimeChange;
  use crate::timetree::convergence::optimizer::TimetreeOptimizer;
  use eyre::Report;
  use parking_lot::Mutex;
  use pretty_assertions::assert_eq;
  use std::collections::BTreeMap;
  use std::sync::Arc;

  #[test]
  fn test_optimizer_converges_when_n_diff_zero() -> Result<(), Report> {
    let graph = helpers::empty_graph();
    let inference = empty_time_inference(&graph);
    let mut optimizer = TimetreeOptimizer::new(5);

    assert!(optimizer.next_iter(&NoopProgress).is_some());
    optimizer.record(
      0,
      0,
      NodeTimeChange::default(),
      &graph,
      &BranchModel::Input,
      &inference,
      None,
      helpers::fixed_clock(),
      &BTreeMap::new(),
      &NoopProgress,
    )?;

    assert!(optimizer.next_iter(&NoopProgress).is_none());
    assert_eq!(1, optimizer.i);
    assert_eq!(1, optimizer.trace.len());

    Ok(())
  }

  #[test]
  fn test_optimizer_continues_while_node_times_move() -> Result<(), Report> {
    let graph = helpers::empty_graph();
    let inference = empty_time_inference(&graph);
    let mut optimizer = TimetreeOptimizer::new(5);

    assert!(optimizer.next_iter(&NoopProgress).is_some());
    optimizer.record(
      0,
      0,
      helpers::moved_by(10.0 * NODE_TIME_TOLERANCE_YEARS),
      &graph,
      &BranchModel::Input,
      &inference,
      None,
      helpers::fixed_clock(),
      &BTreeMap::new(),
      &NoopProgress,
    )?;

    assert!(optimizer.next_iter(&NoopProgress).is_some());
    optimizer.record(
      0,
      0,
      helpers::moved_by(0.1 * NODE_TIME_TOLERANCE_YEARS),
      &graph,
      &BranchModel::Input,
      &inference,
      None,
      helpers::fixed_clock(),
      &BTreeMap::new(),
      &NoopProgress,
    )?;

    assert!(optimizer.next_iter(&NoopProgress).is_none());
    assert_eq!(2, optimizer.i);

    Ok(())
  }

  #[test]
  fn test_optimizer_settled_times_do_not_converge_while_polytomies_resolve() -> Result<(), Report> {
    let graph = helpers::empty_graph();
    let inference = empty_time_inference(&graph);
    let mut optimizer = TimetreeOptimizer::new(5);
    let settled = helpers::moved_by(0.1 * NODE_TIME_TOLERANCE_YEARS);

    assert!(optimizer.next_iter(&NoopProgress).is_some());
    optimizer.record(
      0,
      2,
      settled,
      &graph,
      &BranchModel::Input,
      &inference,
      None,
      helpers::fixed_clock(),
      &BTreeMap::new(),
      &NoopProgress,
    )?;

    assert!(optimizer.next_iter(&NoopProgress).is_some());
    optimizer.record(
      0,
      0,
      settled,
      &graph,
      &BranchModel::Input,
      &inference,
      None,
      helpers::fixed_clock(),
      &BTreeMap::new(),
      &NoopProgress,
    )?;

    assert!(optimizer.next_iter(&NoopProgress).is_none());
    assert_eq!(2, optimizer.i);

    Ok(())
  }

  #[test]
  fn test_optimizer_continues_when_n_diff_positive() -> Result<(), Report> {
    let graph = helpers::empty_graph();
    let inference = empty_time_inference(&graph);
    let mut optimizer = TimetreeOptimizer::new(5);

    assert!(optimizer.next_iter(&NoopProgress).is_some());
    optimizer.record(
      10,
      0,
      NodeTimeChange::default(),
      &graph,
      &BranchModel::Input,
      &inference,
      None,
      helpers::fixed_clock(),
      &BTreeMap::new(),
      &NoopProgress,
    )?;

    assert!(optimizer.next_iter(&NoopProgress).is_some());
    optimizer.record(
      3,
      0,
      NodeTimeChange::default(),
      &graph,
      &BranchModel::Input,
      &inference,
      None,
      helpers::fixed_clock(),
      &BTreeMap::new(),
      &NoopProgress,
    )?;

    assert!(optimizer.next_iter(&NoopProgress).is_some());
    optimizer.record(
      0,
      0,
      NodeTimeChange::default(),
      &graph,
      &BranchModel::Input,
      &inference,
      None,
      helpers::fixed_clock(),
      &BTreeMap::new(),
      &NoopProgress,
    )?;

    assert!(optimizer.next_iter(&NoopProgress).is_none());
    assert_eq!(3, optimizer.i);

    let trace = &optimizer.trace;
    assert_eq!(3, trace.len());
    assert_eq!(10, trace[0].n_diff);
    assert_eq!(3, trace[1].n_diff);
    assert_eq!(0, trace[2].n_diff);

    Ok(())
  }

  #[test]
  fn test_optimizer_stops_at_max_iterations() -> Result<(), Report> {
    let graph = helpers::empty_graph();
    let inference = empty_time_inference(&graph);
    let mut optimizer = TimetreeOptimizer::new(3);

    for _ in 0..3 {
      assert!(optimizer.next_iter(&NoopProgress).is_some());
      optimizer.record(
        10,
        0,
        NodeTimeChange::default(),
        &graph,
        &BranchModel::Input,
        &inference,
        None,
        helpers::fixed_clock(),
        &BTreeMap::new(),
        &NoopProgress,
      )?;
    }

    assert!(optimizer.next_iter(&NoopProgress).is_none());
    assert_eq!(3, optimizer.i);
    assert_eq!(3, optimizer.trace.len());

    Ok(())
  }

  #[test]
  fn test_optimizer_n_resolved_prevents_convergence() -> Result<(), Report> {
    let graph = helpers::empty_graph();
    let inference = empty_time_inference(&graph);
    let mut optimizer = TimetreeOptimizer::new(5);

    assert!(optimizer.next_iter(&NoopProgress).is_some());
    optimizer.record(
      0,
      3,
      NodeTimeChange::default(),
      &graph,
      &BranchModel::Input,
      &inference,
      None,
      helpers::fixed_clock(),
      &BTreeMap::new(),
      &NoopProgress,
    )?;

    assert!(optimizer.next_iter(&NoopProgress).is_some());
    optimizer.record(
      0,
      0,
      NodeTimeChange::default(),
      &graph,
      &BranchModel::Input,
      &inference,
      None,
      helpers::fixed_clock(),
      &BTreeMap::new(),
      &NoopProgress,
    )?;

    assert!(optimizer.next_iter(&NoopProgress).is_none());
    assert_eq!(2, optimizer.i);

    Ok(())
  }

  #[test]
  fn test_optimizer_trace_sink_receives_each_iteration_with_its_clock() -> Result<(), Report> {
    let graph = helpers::empty_graph();
    let inference = empty_time_inference(&graph);
    let records = Arc::new(Mutex::new(vec![]));
    let mut sink = helpers::RecordingSink(Arc::clone(&records));
    let mut optimizer = TimetreeOptimizer::new(3).with_trace_sink(&mut sink);

    assert!(optimizer.next_iter(&NoopProgress).is_some());
    optimizer.record(
      5,
      1,
      NodeTimeChange::default(),
      &graph,
      &BranchModel::Input,
      &inference,
      None,
      IterationClock {
        clock_rate: 2e-3,
        r_squared: Some(0.25),
      },
      &BTreeMap::new(),
      &NoopProgress,
    )?;

    assert!(optimizer.next_iter(&NoopProgress).is_some());
    optimizer.record(
      0,
      0,
      NodeTimeChange::default(),
      &graph,
      &BranchModel::Input,
      &inference,
      None,
      helpers::fixed_clock(),
      &BTreeMap::new(),
      &NoopProgress,
    )?;

    let summary: Vec<(usize, usize, usize, IterationClock)> = records
      .lock()
      .iter()
      .map(|record| {
        (
          record.iteration,
          record.metrics.n_diff,
          record.metrics.n_resolved,
          record.clock,
        )
      })
      .collect();
    assert_eq!(
      vec![
        (
          1,
          5,
          1,
          IterationClock {
            clock_rate: 2e-3,
            r_squared: Some(0.25)
          }
        ),
        (2, 0, 0, helpers::fixed_clock()),
      ],
      summary
    );
    Ok(())
  }

  mod helpers {
    use crate::timetree::convergence::metrics::{IterationClock, IterationRecord};
    use crate::timetree::convergence::node_times::NodeTimeChange;
    use crate::timetree::convergence::optimizer::TraceSink;
    use eyre::Report;
    use parking_lot::Mutex;
    use std::sync::Arc;
    use treetime_graph::graph::Graph;

    pub(super) struct RecordingSink(pub(crate) Arc<Mutex<Vec<IterationRecord>>>);

    impl TraceSink for RecordingSink {
      fn emit(&mut self, record: &IterationRecord) -> Result<(), Report> {
        self.0.lock().push(record.clone());
        Ok(())
      }
    }

    pub(super) const fn fixed_clock() -> IterationClock {
      IterationClock {
        clock_rate: 1e-3,
        r_squared: None,
      }
    }

    pub(super) fn moved_by(years: f64) -> NodeTimeChange {
      NodeTimeChange {
        max: Some(years),
        rms: Some(years),
      }
    }

    pub(super) fn empty_graph() -> Graph {
      let mut graph = Graph::new();
      graph.add_node();
      graph.build().expect("build graph");
      graph
    }
  }
}
