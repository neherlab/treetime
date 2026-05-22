#[cfg(test)]
mod tests {
  use crate::timetree::convergence::optimizer::TimetreeOptimizer;
  use eyre::Report;
  use pretty_assertions::assert_eq;

  #[test]
  fn test_optimizer_converges_when_n_diff_zero() -> Result<(), Report> {
    let graph = helpers::empty_graph();
    let mut optimizer = TimetreeOptimizer::new(5, false);

    assert!(optimizer.next_iter().is_some());
    optimizer.record(0, 0, &graph, &[], None)?;

    assert!(optimizer.next_iter().is_none());
    assert_eq!(1, optimizer.iteration_count());
    assert_eq!(1, optimizer.trace().len());

    Ok(())
  }

  #[allow(clippy::missing_asserts_for_indexing)]
  #[test]
  fn test_optimizer_continues_when_n_diff_positive() -> Result<(), Report> {
    let graph = helpers::empty_graph();
    let mut optimizer = TimetreeOptimizer::new(5, false);

    assert!(optimizer.next_iter().is_some());
    optimizer.record(10, 0, &graph, &[], None)?;

    assert!(optimizer.next_iter().is_some());
    optimizer.record(3, 0, &graph, &[], None)?;

    assert!(optimizer.next_iter().is_some());
    optimizer.record(0, 0, &graph, &[], None)?;

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
      optimizer.record(10, 0, &graph, &[], None)?;
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
    optimizer.record(0, 3, &graph, &[], None)?;

    assert!(optimizer.next_iter().is_some());
    optimizer.record(0, 0, &graph, &[], None)?;

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
    optimizer.record(5, 1, &graph, &[], None)?;

    assert!(optimizer.next_iter().is_some());
    optimizer.record(0, 0, &graph, &[], None)?;

    assert_eq!(2, optimizer.trace().len());
    Ok(())
  }

  mod helpers {
    use crate::partition::timetree::GraphTimetree;
    use crate::payload::timetree::NodeTimetree;

    pub fn empty_graph() -> GraphTimetree {
      let mut graph = GraphTimetree::new();
      graph.add_node(NodeTimetree::default());
      graph.build().expect("build graph");
      graph
    }
  }
}
