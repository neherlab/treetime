#[cfg(test)]
mod tests {
  use crate::timetree::convergence::metrics::NODE_TIME_TOLERANCE_YEARS;
  use crate::timetree::convergence::node_times::NodeTimeChange;
  use crate::timetree::convergence::optimizer::TimetreeOptimizer;
  use crate::timetree::timetree_state::TimetreeState;
  use eyre::Report;
  use pretty_assertions::assert_eq;
  use std::sync::Arc;
  use std::sync::atomic::{AtomicUsize, Ordering};

  #[test]
  fn test_optimizer_converges_when_n_diff_zero() -> Result<(), Report> {
    let graph = helpers::empty_graph();
    let state = TimetreeState::new(&graph);
    let mut optimizer = TimetreeOptimizer::new(5, false);

    assert!(optimizer.next_iter().is_some());
    optimizer.record(0, 0, NodeTimeChange::default(), &graph, &[], &state, None)?;

    assert!(optimizer.next_iter().is_none());
    assert_eq!(1, optimizer.i);
    assert_eq!(1, optimizer.trace.len());

    Ok(())
  }

  #[test]
  fn test_optimizer_continues_while_node_times_move() -> Result<(), Report> {
    let graph = helpers::empty_graph();
    let state = TimetreeState::new(&graph);
    let mut optimizer = TimetreeOptimizer::new(5, false);

    assert!(optimizer.next_iter().is_some());
    optimizer.record(
      0,
      0,
      helpers::moved_by(10.0 * NODE_TIME_TOLERANCE_YEARS),
      &graph,
      &[],
      &state,
      None,
    )?;

    assert!(optimizer.next_iter().is_some());
    optimizer.record(
      0,
      0,
      helpers::moved_by(0.1 * NODE_TIME_TOLERANCE_YEARS),
      &graph,
      &[],
      &state,
      None,
    )?;

    assert!(optimizer.next_iter().is_none());
    assert_eq!(2, optimizer.i);

    Ok(())
  }

  #[test]
  fn test_optimizer_settled_times_do_not_converge_while_polytomies_resolve() -> Result<(), Report> {
    let graph = helpers::empty_graph();
    let state = TimetreeState::new(&graph);
    let mut optimizer = TimetreeOptimizer::new(5, false);
    let settled = helpers::moved_by(0.1 * NODE_TIME_TOLERANCE_YEARS);

    assert!(optimizer.next_iter().is_some());
    optimizer.record(0, 2, settled, &graph, &[], &state, None)?;

    assert!(optimizer.next_iter().is_some());
    optimizer.record(0, 0, settled, &graph, &[], &state, None)?;

    assert!(optimizer.next_iter().is_none());
    assert_eq!(2, optimizer.i);

    Ok(())
  }

  #[test]
  fn test_optimizer_continues_when_n_diff_positive() -> Result<(), Report> {
    let graph = helpers::empty_graph();
    let state = TimetreeState::new(&graph);
    let mut optimizer = TimetreeOptimizer::new(5, false);

    assert!(optimizer.next_iter().is_some());
    optimizer.record(10, 0, NodeTimeChange::default(), &graph, &[], &state, None)?;

    assert!(optimizer.next_iter().is_some());
    optimizer.record(3, 0, NodeTimeChange::default(), &graph, &[], &state, None)?;

    assert!(optimizer.next_iter().is_some());
    optimizer.record(0, 0, NodeTimeChange::default(), &graph, &[], &state, None)?;

    assert!(optimizer.next_iter().is_none());
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
    let state = TimetreeState::new(&graph);
    let mut optimizer = TimetreeOptimizer::new(3, false);

    for _ in 0..3 {
      assert!(optimizer.next_iter().is_some());
      optimizer.record(10, 0, NodeTimeChange::default(), &graph, &[], &state, None)?;
    }

    assert!(optimizer.next_iter().is_none());
    assert_eq!(3, optimizer.i);
    assert_eq!(3, optimizer.trace.len());

    Ok(())
  }

  #[test]
  fn test_optimizer_n_resolved_prevents_convergence() -> Result<(), Report> {
    let graph = helpers::empty_graph();
    let state = TimetreeState::new(&graph);
    let mut optimizer = TimetreeOptimizer::new(5, false);

    assert!(optimizer.next_iter().is_some());
    optimizer.record(0, 3, NodeTimeChange::default(), &graph, &[], &state, None)?;

    assert!(optimizer.next_iter().is_some());
    optimizer.record(0, 0, NodeTimeChange::default(), &graph, &[], &state, None)?;

    assert!(optimizer.next_iter().is_none());
    assert_eq!(2, optimizer.i);

    Ok(())
  }

  #[test]
  fn test_optimizer_trace_sink_receives_each_iteration() -> Result<(), Report> {
    let graph = helpers::empty_graph();
    let state = TimetreeState::new(&graph);
    let count = Arc::new(AtomicUsize::new(0));
    let mut optimizer =
      TimetreeOptimizer::new(3, false).with_trace_sink(Box::new(helpers::CountingSink(Arc::clone(&count))));

    assert!(optimizer.next_iter().is_some());
    optimizer.record(5, 1, NodeTimeChange::default(), &graph, &[], &state, None)?;

    assert!(optimizer.next_iter().is_some());
    optimizer.record(0, 0, NodeTimeChange::default(), &graph, &[], &state, None)?;

    assert_eq!(2, optimizer.trace.len());
    assert_eq!(2, count.load(Ordering::Relaxed));
    Ok(())
  }

  mod helpers {
    use crate::timetree::convergence::metrics::ConvergenceMetrics;
    use crate::timetree::convergence::node_times::NodeTimeChange;
    use crate::timetree::convergence::optimizer::TraceSink;
    use eyre::Report;
    use std::sync::Arc;
    use std::sync::atomic::{AtomicUsize, Ordering};
    use treetime_graph::graph::Graph;

    pub(super) struct CountingSink(pub(crate) Arc<AtomicUsize>);

    impl TraceSink for CountingSink {
      fn emit(&mut self, _metric: &ConvergenceMetrics) -> Result<(), Report> {
        self.0.fetch_add(1, Ordering::Relaxed);
        Ok(())
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
