#[cfg(test)]
mod tests {
  use crate::timetree::convergence::metrics::NODE_TIME_TOLERANCE_YEARS;
  use crate::timetree::convergence::node_times::NodeTimeChange;
  use crate::timetree::convergence::optimizer::TimetreeOptimizer;
  use eyre::Report;
  use pretty_assertions::assert_eq;

  /// No node is dated on both sides of the iteration, so no movement is measurable and the
  /// criterion falls back to the ancestral-sequence count.
  #[test]
  fn test_optimizer_converges_when_n_diff_zero() -> Result<(), Report> {
    let graph = helpers::empty_graph();
    let mut optimizer = TimetreeOptimizer::new(5, false);

    assert!(optimizer.next_iter().is_some());
    optimizer.record(0, 0, NodeTimeChange::default(), &graph, &[], None)?;

    assert!(optimizer.next_iter().is_none());
    assert_eq!(1, optimizer.iteration_count());
    assert_eq!(1, optimizer.trace().len());

    Ok(())
  }

  /// Node times are the primary signal: a round that moved them is not converged even though the
  /// reconstructed sequences are identical, which is the case `n_diff` alone cannot see.
  #[test]
  fn test_optimizer_continues_while_node_times_move() -> Result<(), Report> {
    let graph = helpers::empty_graph();
    let mut optimizer = TimetreeOptimizer::new(5, false);

    assert!(optimizer.next_iter().is_some());
    optimizer.record(
      0,
      0,
      helpers::moved_by(10.0 * NODE_TIME_TOLERANCE_YEARS),
      &graph,
      &[],
      None,
    )?;

    assert!(optimizer.next_iter().is_some());
    optimizer.record(
      0,
      0,
      helpers::moved_by(0.1 * NODE_TIME_TOLERANCE_YEARS),
      &graph,
      &[],
      None,
    )?;

    assert!(optimizer.next_iter().is_none());
    assert_eq!(2, optimizer.iteration_count());

    Ok(())
  }

  /// Settled times are not enough on their own: a round that resolved a polytomy changed the
  /// tree, so the next round has to run whatever the times did.
  #[test]
  fn test_optimizer_settled_times_do_not_converge_while_polytomies_resolve() -> Result<(), Report> {
    let graph = helpers::empty_graph();
    let mut optimizer = TimetreeOptimizer::new(5, false);
    let settled = helpers::moved_by(0.1 * NODE_TIME_TOLERANCE_YEARS);

    assert!(optimizer.next_iter().is_some());
    optimizer.record(0, 2, settled, &graph, &[], None)?;

    assert!(optimizer.next_iter().is_some());
    optimizer.record(0, 0, settled, &graph, &[], None)?;

    assert!(optimizer.next_iter().is_none());
    assert_eq!(2, optimizer.iteration_count());

    Ok(())
  }

  #[allow(clippy::missing_asserts_for_indexing)]
  #[test]
  fn test_optimizer_continues_when_n_diff_positive() -> Result<(), Report> {
    let graph = helpers::empty_graph();
    let mut optimizer = TimetreeOptimizer::new(5, false);

    assert!(optimizer.next_iter().is_some());
    optimizer.record(10, 0, NodeTimeChange::default(), &graph, &[], None)?;

    assert!(optimizer.next_iter().is_some());
    optimizer.record(3, 0, NodeTimeChange::default(), &graph, &[], None)?;

    assert!(optimizer.next_iter().is_some());
    optimizer.record(0, 0, NodeTimeChange::default(), &graph, &[], None)?;

    assert!(optimizer.next_iter().is_none());
    assert_eq!(3, optimizer.iteration_count());

    let trace = optimizer.trace();
    assert_eq!(3, trace.len());
    assert_eq!(10, trace[0].n_diff);
    assert_eq!(3, trace[1].n_diff);
    assert_eq!(0, trace[2].n_diff);

    Ok(())
  }

  #[test]
  fn test_optimizer_stops_at_max_iterations() -> Result<(), Report> {
    let graph = helpers::empty_graph();
    let mut optimizer = TimetreeOptimizer::new(3, false);

    for _ in 0..3 {
      assert!(optimizer.next_iter().is_some());
      optimizer.record(10, 0, NodeTimeChange::default(), &graph, &[], None)?;
    }

    assert!(optimizer.next_iter().is_none());
    assert_eq!(3, optimizer.iteration_count());
    assert_eq!(3, optimizer.trace().len());

    Ok(())
  }

  #[test]
  fn test_optimizer_n_resolved_prevents_convergence() -> Result<(), Report> {
    let graph = helpers::empty_graph();
    let mut optimizer = TimetreeOptimizer::new(5, false);

    assert!(optimizer.next_iter().is_some());
    optimizer.record(0, 3, NodeTimeChange::default(), &graph, &[], None)?;

    assert!(optimizer.next_iter().is_some());
    optimizer.record(0, 0, NodeTimeChange::default(), &graph, &[], None)?;

    assert!(optimizer.next_iter().is_none());
    assert_eq!(2, optimizer.iteration_count());

    Ok(())
  }

  #[test]
  fn test_optimizer_tracelog_writes_csv() -> Result<(), Report> {
    let graph = helpers::empty_graph();
    let buf = Vec::<u8>::new();
    let mut optimizer = TimetreeOptimizer::new(3, false).with_tracelog(buf)?;

    assert!(optimizer.next_iter().is_some());
    optimizer.record(5, 1, NodeTimeChange::default(), &graph, &[], None)?;

    assert!(optimizer.next_iter().is_some());
    optimizer.record(0, 0, NodeTimeChange::default(), &graph, &[], None)?;

    assert_eq!(2, optimizer.trace().len());
    Ok(())
  }

  mod helpers {
    use crate::partition::timetree::partition::GraphTimetree;
    use crate::payload::timetree::NodeTimetree;
    use crate::timetree::convergence::node_times::NodeTimeChange;

    /// A measured movement of `years`, as a single node moving that far would produce.
    pub fn moved_by(years: f64) -> NodeTimeChange {
      NodeTimeChange {
        max: Some(years),
        rms: Some(years),
      }
    }

    pub fn empty_graph() -> GraphTimetree {
      let mut graph = GraphTimetree::new();
      graph.add_node(NodeTimetree::default());
      graph.build().expect("build graph");
      graph
    }
  }
}
