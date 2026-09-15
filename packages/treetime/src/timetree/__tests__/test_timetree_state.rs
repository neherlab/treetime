#[cfg(test)]
mod tests {
  use crate::timetree::timetree_state::TimetreeState;
  use eyre::Report;
  use pretty_assertions::assert_eq;
  use std::sync::Arc;
  use treetime_distribution::Distribution;
  use treetime_io::nwk::nwk_read_str;

  /// A topology change clears the value-resident branch-length distribution and backward message on
  /// every edge, so the next branch-distribution build starts each surviving edge from scratch.
  #[test]
  fn test_timetree_state_reset_date_edges_clears_distribution_and_message() -> Result<(), Report> {
    let graph = nwk_read_str("((A:1.0,B:1.0)I:1.0)root;")?.graph;
    let mut state = TimetreeState::new(&graph);
    for edge_ref in graph.get_edges() {
      let key = edge_ref.key();
      let entry = state.edge_mut(key);
      entry.branch_length_distribution = Some(Arc::new(Distribution::point(1.0, 0.0)));
      entry.msg_to_parent = Some(Arc::new(Distribution::point(2.0, 0.0)));
    }

    state.reset_date_edges_for_topology_change(&graph);

    for edge_ref in graph.get_edges() {
      let key = edge_ref.key();
      let entry = state.edge(key);
      assert_eq!(None, entry.branch_length_distribution);
      assert_eq!(None, entry.msg_to_parent);
    }

    Ok(())
  }

  /// Rebuilding the maps keeps the value-resident branch-length distribution and backward message for
  /// edges already in the state, so they carry across passes.
  #[test]
  fn test_timetree_state_reseed_preserves_distribution_and_message() -> Result<(), Report> {
    let graph = nwk_read_str("((A:1.0,B:1.0)I:1.0)root;")?.graph;
    let mut state = TimetreeState::new(&graph);
    let key = graph
      .get_edges()
      .collect::<Vec<_>>()
      .first()
      .expect("tree has at least one edge")
      .key();
    let dist = Arc::new(Distribution::point(3.0, 0.0));
    let msg = Arc::new(Distribution::point(4.0, 0.0));
    {
      let entry = state.edge_mut(key);
      entry.branch_length_distribution = Some(Arc::clone(&dist));
      entry.msg_to_parent = Some(Arc::clone(&msg));
    }

    state.reseed_from_values(&graph);

    let entry = state.edge(key);
    assert_eq!(Some(dist), entry.branch_length_distribution);
    assert_eq!(Some(msg), entry.msg_to_parent);

    Ok(())
  }
}
