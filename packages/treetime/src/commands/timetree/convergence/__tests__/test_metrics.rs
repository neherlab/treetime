#[cfg(test)]
mod tests {
  use crate::commands::timetree::convergence::metrics::TimetreeOptimizer;
  use crate::representation::partition::timetree::GraphTimetree;
  use crate::representation::payload::timetree::NodeTimetree;
  use eyre::Report;
  use pretty_assertions::assert_eq;
  use std::path::PathBuf;

  #[test]
  fn test_metrics_optimizer_converges_when_n_diff_zero() -> Result<(), Report> {
    let graph = empty_graph();
    let mut optimizer = TimetreeOptimizer::new(5, None::<PathBuf>)?;

    assert!(optimizer.next_iter().is_some());
    optimizer.record(0, 0, &graph, &[])?;

    assert!(optimizer.next_iter().is_none());
    assert_eq!(1, optimizer.iteration_count());
    assert_eq!(1, optimizer.trace().len());

    Ok(())
  }

  #[test]
  fn test_metrics_optimizer_continues_when_n_diff_positive() -> Result<(), Report> {
    let graph = empty_graph();
    let mut optimizer = TimetreeOptimizer::new(5, None::<PathBuf>)?;

    // Iteration 1: sequences still changing
    assert!(optimizer.next_iter().is_some());
    optimizer.record(10, 0, &graph, &[])?;

    // Iteration 2: fewer changes
    assert!(optimizer.next_iter().is_some());
    optimizer.record(3, 0, &graph, &[])?;

    // Iteration 3: converged
    assert!(optimizer.next_iter().is_some());
    optimizer.record(0, 0, &graph, &[])?;

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
  fn test_metrics_optimizer_stops_at_max_iterations() -> Result<(), Report> {
    let graph = empty_graph();
    let mut optimizer = TimetreeOptimizer::new(3, None::<PathBuf>)?;

    for _ in 0..3 {
      assert!(optimizer.next_iter().is_some());
      optimizer.record(10, 0, &graph, &[])?;
    }

    assert!(optimizer.next_iter().is_none());
    assert_eq!(3, optimizer.iteration_count());
    assert_eq!(3, optimizer.trace().len());

    Ok(())
  }

  #[test]
  fn test_metrics_optimizer_n_resolved_prevents_convergence() -> Result<(), Report> {
    let graph = empty_graph();
    let mut optimizer = TimetreeOptimizer::new(5, None::<PathBuf>)?;

    // n_diff=0 but polytomies resolved: not converged
    assert!(optimizer.next_iter().is_some());
    optimizer.record(0, 3, &graph, &[])?;

    // Both zero: converged
    assert!(optimizer.next_iter().is_some());
    optimizer.record(0, 0, &graph, &[])?;

    assert!(optimizer.next_iter().is_none());
    assert_eq!(2, optimizer.iteration_count());

    Ok(())
  }

  fn empty_graph() -> GraphTimetree {
    let mut graph = GraphTimetree::new();
    graph.add_node(NodeTimetree::default());
    graph.build().expect("build graph");
    graph
  }
}
