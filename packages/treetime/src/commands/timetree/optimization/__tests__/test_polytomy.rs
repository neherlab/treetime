#[cfg(test)]
mod tests {
  use crate::commands::timetree::optimization::polytomy::{
    resolve_polytomies_with_options, DEFAULT_RESOLUTION_THRESHOLD,
  };
  use crate::representation::partition::timetree::GraphTimetree;
  use crate::test_utils::find_node_key_by_name;
  use eyre::Report;
  use std::sync::Arc;
  use treetime_distribution::Distribution;
  use treetime_graph::node::Named;
  use treetime_io::nwk::nwk_read_str;

  /// Create a tree with a polytomy (node with 3+ children).
  /// Tree: ((A,B,C)ABC)root where ABC has 3 children
  fn create_polytomy_tree() -> Result<GraphTimetree, Report> {
    // Newick with 3 children at ABC node
    let graph: GraphTimetree = nwk_read_str("((A:0.1,B:0.2,C:0.15)ABC:0.05)root;")?;

    // Set times for all nodes (tips have explicit times, internals derived)
    let tip_times = [("A", 2020.0), ("B", 2015.0), ("C", 2018.0)];
    for (name, time) in tip_times {
      let key = find_node_key_by_name(&graph, name).ok_or_else(|| eyre::eyre!("{name} not found"))?;
      let node = graph.get_node(key).expect("Node must exist");
      node.write_arc().payload().write_arc().time = Some(time);
    }

    // Set internal node times
    if let Some(abc_key) = find_node_key_by_name(&graph, "ABC") {
      let node = graph.get_node(abc_key).expect("Node must exist");
      node.write_arc().payload().write_arc().time = Some(2010.0);
    }

    if let Some(root_key) = find_node_key_by_name(&graph, "root") {
      let node = graph.get_node(root_key).expect("Node must exist");
      node.write_arc().payload().write_arc().time = Some(2000.0);
    }

    // Set branch length distributions on edges
    for edge in graph.get_edges() {
      let edge = edge.write_arc();
      let mut payload = edge.payload().write_arc();
      // Simple uniform distribution for branch lengths
      payload.branch_length_distribution = Some(Arc::new(Distribution::point(0.1, 1.0)));
    }

    Ok(graph)
  }

  /// Create a binary tree (no polytomies).
  fn create_binary_tree() -> Result<GraphTimetree, Report> {
    let graph: GraphTimetree = nwk_read_str("((A:0.1,B:0.2)AB:0.05,(C:0.15,D:0.1)CD:0.08)root;")?;

    let tip_times = [("A", 2020.0), ("B", 2015.0), ("C", 2018.0), ("D", 2012.0)];
    for (name, time) in tip_times {
      if let Some(key) = find_node_key_by_name(&graph, name) {
        let node = graph.get_node(key).expect("Node must exist");
        node.write_arc().payload().write_arc().time = Some(time);
      }
    }

    Ok(graph)
  }

  #[test]
  fn test_find_polytomy_nodes_detects_multifurcation() -> Result<(), Report> {
    let graph = create_polytomy_tree()?;

    // Count nodes with >2 children
    let polytomy_count = graph
      .get_nodes()
      .into_iter()
      .filter(|n| n.read_arc().degree_out() > 2)
      .count();

    assert_eq!(polytomy_count, 1, "Should detect exactly one polytomy (ABC node)");
    Ok(())
  }

  #[test]
  fn test_find_polytomy_nodes_returns_empty_for_binary_tree() -> Result<(), Report> {
    let graph = create_binary_tree()?;

    let polytomy_count = graph
      .get_nodes()
      .into_iter()
      .filter(|n| n.read_arc().degree_out() > 2)
      .count();

    assert_eq!(polytomy_count, 0, "Binary tree should have no polytomies");
    Ok(())
  }

  #[test]
  fn test_resolve_polytomies_reduces_children_count() -> Result<(), Report> {
    let mut graph = create_polytomy_tree()?;

    // Find ABC node and verify it starts with 3 children
    let abc_key = find_node_key_by_name(&graph, "ABC").ok_or_else(|| eyre::eyre!("ABC not found"))?;
    let initial_children = graph.get_node(abc_key).expect("Node must exist").read_arc().degree_out();
    assert_eq!(initial_children, 3, "ABC should start with 3 children");

    // Resolve with very low threshold (should merge all possible pairs)
    let partitions = vec![];
    let n_resolved = resolve_polytomies_with_options(&mut graph, &partitions, -1000.0, false)?;

    // After resolution, ABC should have 2 children (one merge happened)
    let final_children = graph.get_node(abc_key).expect("Node must exist").read_arc().degree_out();
    assert_eq!(final_children, 2, "ABC should have 2 children after resolution");
    assert_eq!(n_resolved, 1, "Should have created 1 new node");

    Ok(())
  }

  #[test]
  fn test_resolve_polytomies_no_change_for_binary_tree() -> Result<(), Report> {
    let mut graph = create_binary_tree()?;

    let initial_node_count = graph.get_nodes().len();
    let partitions = vec![];
    let n_resolved = resolve_polytomies_with_options(&mut graph, &partitions, DEFAULT_RESOLUTION_THRESHOLD, false)?;

    assert_eq!(n_resolved, 0, "Binary tree should have no resolutions");
    assert_eq!(
      graph.get_nodes().len(),
      initial_node_count,
      "Node count should remain unchanged"
    );

    Ok(())
  }

  #[test]
  fn test_resolve_polytomies_respects_threshold() -> Result<(), Report> {
    let mut graph = create_polytomy_tree()?;

    // With very high threshold, no merges should occur
    let partitions = vec![];
    let n_resolved = resolve_polytomies_with_options(&mut graph, &partitions, 1000.0, false)?;

    assert_eq!(n_resolved, 0, "Very high threshold should prevent any merges");

    // ABC should still have 3 children
    let abc_key = find_node_key_by_name(&graph, "ABC").ok_or_else(|| eyre::eyre!("ABC not found"))?;
    let children = graph.get_node(abc_key).expect("Node must exist").read_arc().degree_out();
    assert_eq!(children, 3, "ABC should still have 3 children");

    Ok(())
  }

  #[test]
  fn test_resolve_polytomies_new_node_has_correct_time() -> Result<(), Report> {
    let mut graph = create_polytomy_tree()?;

    let abc_key = find_node_key_by_name(&graph, "ABC").ok_or_else(|| eyre::eyre!("ABC not found"))?;
    let abc_time = graph
      .get_node(abc_key)
      .expect("Node must exist")
      .read_arc()
      .payload()
      .read_arc()
      .time
      .unwrap_or(0.0);

    // Get children times before resolution
    let children_times: Vec<f64> = graph
      .get_node(abc_key)
      .expect("Node must exist")
      .read_arc()
      .outbound()
      .iter()
      .filter_map(|&edge_key| {
        let edge = graph.get_edge(edge_key)?;
        let child_key = edge.read_arc().target();
        let child = graph.get_node(child_key)?;
        child.read_arc().payload().read_arc().time
      })
      .collect();

    // Resolve polytomy
    let partitions = vec![];
    resolve_polytomies_with_options(&mut graph, &partitions, -1000.0, false)?;

    // Find the new internal node (not ABC, not a leaf, not root)
    // New nodes created by polytomy resolution don't have a name set
    let new_nodes: Vec<_> = graph
      .get_nodes()
      .into_iter()
      .filter(|n| {
        let n = n.read_arc();
        let payload = n.payload().read_arc();
        let has_name = payload.name().is_some();
        !n.is_leaf() && !n.is_root() && !has_name
      })
      .collect();

    assert_eq!(new_nodes.len(), 1, "Should have exactly one new unnamed internal node");

    let new_node_time = new_nodes[0].read_arc().payload().read_arc().time.unwrap_or(0.0);

    // In calendar time: parent (ABC=2010) < new_node < min(children)
    // The new node must be between the parent and the closest child
    let min_child_time = children_times.iter().copied().fold(f64::INFINITY, f64::min);

    assert!(
      new_node_time > abc_time,
      "New node time ({new_node_time}) should be > parent time ({abc_time})"
    );
    assert!(
      new_node_time < min_child_time,
      "New node time ({new_node_time}) should be < min child time ({min_child_time})"
    );

    Ok(())
  }

  #[test]
  fn test_resolve_polytomies_large_polytomy() -> Result<(), Report> {
    // Create a tree with a 5-way polytomy
    let mut graph: GraphTimetree = nwk_read_str("((A:0.1,B:0.2,C:0.15,D:0.12,E:0.18)ABCDE:0.05)root;")?;

    // Set times
    let tip_times = [
      ("A", 2020.0),
      ("B", 2015.0),
      ("C", 2018.0),
      ("D", 2012.0),
      ("E", 2019.0),
    ];
    for (name, time) in tip_times {
      if let Some(key) = find_node_key_by_name(&graph, name) {
        let node = graph.get_node(key).expect("Node must exist");
        node.write_arc().payload().write_arc().time = Some(time);
      }
    }

    if let Some(key) = find_node_key_by_name(&graph, "ABCDE") {
      let node = graph.get_node(key).expect("Node must exist");
      node.write_arc().payload().write_arc().time = Some(2005.0);
    }

    if let Some(key) = find_node_key_by_name(&graph, "root") {
      let node = graph.get_node(key).expect("Node must exist");
      node.write_arc().payload().write_arc().time = Some(2000.0);
    }

    // Set branch length distributions
    for edge in graph.get_edges() {
      let edge = edge.write_arc();
      let mut payload = edge.payload().write_arc();
      payload.branch_length_distribution = Some(Arc::new(Distribution::point(0.1, 1.0)));
    }

    let abcde_key = find_node_key_by_name(&graph, "ABCDE").ok_or_else(|| eyre::eyre!("ABCDE not found"))?;
    let initial_children = graph.get_node(abcde_key).expect("Node must exist").read_arc().degree_out();
    assert_eq!(initial_children, 5, "ABCDE should start with 5 children");

    // Resolve with very low threshold
    let partitions = vec![];
    let n_resolved = resolve_polytomies_with_options(&mut graph, &partitions, -1000.0, false)?;

    // 5-way polytomy needs 3 merges to become binary (5->4->3->2)
    assert_eq!(n_resolved, 3, "Should create 3 new nodes to resolve 5-way polytomy");

    // ABCDE should now have 2 children
    let final_children = graph.get_node(abcde_key).expect("Node must exist").read_arc().degree_out();
    assert_eq!(final_children, 2, "ABCDE should have 2 children after full resolution");

    Ok(())
  }
}
