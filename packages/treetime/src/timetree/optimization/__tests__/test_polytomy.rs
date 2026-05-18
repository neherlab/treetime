#[cfg(test)]
mod tests {
  use crate::representation::partition::timetree::GraphTimetree;
  use crate::test_utils::find_node_key_by_name;
  use crate::timetree::optimization::polytomy::{
    DEFAULT_RESOLUTION_THRESHOLD, prepare_tree_after_topology_change, resolve_polytomies_with_options,
  };
  use eyre::Report;
  use ndarray::Array1;
  use std::sync::Arc;
  use treetime_distribution::Distribution;
  use treetime_graph::edge::HasBranchLength;
  use treetime_graph::node::Named;
  use treetime_io::nwk::nwk_read_str;
  use treetime_utils::make_report;

  const TEST_CLOCK_RATE: f64 = 0.001;

  #[test]
  fn test_find_polytomy_nodes_detects_multifurcation() -> Result<(), Report> {
    let graph = helpers::create_polytomy_tree()?;

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
    let graph = helpers::create_binary_tree()?;

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
    let mut graph = helpers::create_polytomy_tree()?;

    // Find ABC node and verify it starts with 3 children
    let abc_key = find_node_key_by_name(&graph, "ABC").ok_or_else(|| make_report!("ABC not found"))?;
    let initial_children = graph
      .get_node(abc_key)
      .expect("Node must exist")
      .read_arc()
      .degree_out();
    assert_eq!(initial_children, 3, "ABC should start with 3 children");

    // Resolve with very low threshold (should merge all possible pairs)
    let partitions = vec![];
    let n_resolved = resolve_polytomies_with_options(&mut graph, &partitions, -1000.0, 10.0, TEST_CLOCK_RATE, false)?;

    // After resolution, ABC should have 2 children (one merge happened)
    let final_children = graph
      .get_node(abc_key)
      .expect("Node must exist")
      .read_arc()
      .degree_out();
    assert_eq!(final_children, 2, "ABC should have 2 children after resolution");
    assert_eq!(n_resolved, 1, "Should have created 1 new node");

    Ok(())
  }

  #[test]
  fn test_resolve_polytomies_no_change_for_binary_tree() -> Result<(), Report> {
    let mut graph = helpers::create_binary_tree()?;

    let initial_node_count = graph.get_nodes().len();
    let partitions = vec![];
    let n_resolved = resolve_polytomies_with_options(
      &mut graph,
      &partitions,
      DEFAULT_RESOLUTION_THRESHOLD,
      10.0,
      TEST_CLOCK_RATE,
      false,
    )?;

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
    let mut graph = helpers::create_polytomy_tree()?;

    // With very high threshold, no merges should occur
    let partitions = vec![];
    let n_resolved = resolve_polytomies_with_options(&mut graph, &partitions, 1000.0, 10.0, TEST_CLOCK_RATE, false)?;

    assert_eq!(n_resolved, 0, "Very high threshold should prevent any merges");

    // ABC should still have 3 children
    let abc_key = find_node_key_by_name(&graph, "ABC").ok_or_else(|| make_report!("ABC not found"))?;
    let children = graph
      .get_node(abc_key)
      .expect("Node must exist")
      .read_arc()
      .degree_out();
    assert_eq!(children, 3, "ABC should still have 3 children");

    Ok(())
  }

  #[test]
  fn test_resolve_polytomies_new_node_has_correct_time() -> Result<(), Report> {
    let mut graph = helpers::create_polytomy_tree()?;

    let abc_key = find_node_key_by_name(&graph, "ABC").ok_or_else(|| make_report!("ABC not found"))?;
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
    resolve_polytomies_with_options(&mut graph, &partitions, -1000.0, 10.0, TEST_CLOCK_RATE, false)?;

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

    for edge in graph.get_edges() {
      let edge = edge.write_arc();
      let mut payload = edge.payload().write_arc();
      payload.branch_length_distribution = Some(Arc::new(Distribution::point(0.1, 1.0)));
      payload.time_length = payload.branch_length();
      payload.set_branch_length(Some(0.0));
    }

    let abcde_key = find_node_key_by_name(&graph, "ABCDE").ok_or_else(|| make_report!("ABCDE not found"))?;
    let initial_children = graph
      .get_node(abcde_key)
      .expect("Node must exist")
      .read_arc()
      .degree_out();
    assert_eq!(initial_children, 5, "ABCDE should start with 5 children");

    // Resolve with very low threshold
    let partitions = vec![];
    let n_resolved = resolve_polytomies_with_options(&mut graph, &partitions, -1000.0, 10.0, TEST_CLOCK_RATE, false)?;

    // 5-way polytomy needs 3 merges to become binary (5->4->3->2)
    assert_eq!(n_resolved, 3, "Should create 3 new nodes to resolve 5-way polytomy");

    // ABCDE should now have 2 children
    let final_children = graph
      .get_node(abcde_key)
      .expect("Node must exist")
      .read_arc()
      .degree_out();
    assert_eq!(final_children, 2, "ABCDE should have 2 children after full resolution");

    Ok(())
  }

  /// Higher zero_branch_slope increases the penalty for zero-mutation branches,
  /// reducing cost gain and preventing merges that succeed with lower slopes.
  #[test]
  fn test_resolve_polytomies_zero_branch_slope_affects_merge_gain() -> Result<(), Report> {
    let partitions = vec![];

    // Large slope: penalty dominates branch length improvement → no merge
    let mut graph_high = helpers::create_polytomy_tree_with_realistic_distributions()?;
    let n_high = resolve_polytomies_with_options(
      &mut graph_high,
      &partitions,
      DEFAULT_RESOLUTION_THRESHOLD,
      1000.0,
      TEST_CLOCK_RATE,
      false,
    )?;

    // Small slope: penalty negligible, branch redistribution drives merge
    let mut graph_low = helpers::create_polytomy_tree_with_realistic_distributions()?;
    let n_low = resolve_polytomies_with_options(
      &mut graph_low,
      &partitions,
      DEFAULT_RESOLUTION_THRESHOLD,
      0.01,
      TEST_CLOCK_RATE,
      false,
    )?;

    assert!(
      n_low > n_high,
      "Lower zero_branch_slope ({n_low} merges) should allow more merges than higher slope ({n_high} merges)"
    );

    Ok(())
  }

  /// With zero_branch_slope = 0, penalty vanishes entirely (no sequence data).
  #[test]
  fn test_resolve_polytomies_zero_slope_no_penalty() -> Result<(), Report> {
    let mut graph = helpers::create_polytomy_tree_with_realistic_distributions()?;
    let partitions = vec![];
    let n_resolved = resolve_polytomies_with_options(
      &mut graph,
      &partitions,
      DEFAULT_RESOLUTION_THRESHOLD,
      0.0,
      TEST_CLOCK_RATE,
      false,
    )?;

    assert_eq!(n_resolved, 1, "Zero slope should allow merge");
    Ok(())
  }

  #[test]
  fn test_resolve_polytomies_compressed_children_skipped_by_default() -> Result<(), Report> {
    let mut graph = helpers::create_compressed_polytomy_tree()?;
    let partitions = vec![];
    let n_resolved = resolve_polytomies_with_options(
      &mut graph,
      &partitions,
      DEFAULT_RESOLUTION_THRESHOLD,
      0.01,
      TEST_CLOCK_RATE,
      false,
    )?;

    assert_eq!(
      n_resolved, 0,
      "Compressed children should not be merged when merge_compressed=false"
    );
    Ok(())
  }

  #[test]
  fn test_resolve_polytomies_compressed_children_merged_when_enabled() -> Result<(), Report> {
    let mut graph = helpers::create_compressed_polytomy_tree()?;
    let partitions = vec![];
    let n_resolved = resolve_polytomies_with_options(&mut graph, &partitions, -1000.0, 0.01, TEST_CLOCK_RATE, true)?;

    assert!(
      n_resolved > 0,
      "Compressed children should be merged when merge_compressed=true"
    );
    Ok(())
  }

  #[test]
  fn test_resolve_polytomies_stretched_children_merged_by_default() -> Result<(), Report> {
    let mut graph = helpers::create_polytomy_tree()?;
    let partitions = vec![];
    let n_resolved = resolve_polytomies_with_options(&mut graph, &partitions, -1000.0, 10.0, TEST_CLOCK_RATE, false)?;

    assert!(n_resolved > 0, "Stretched children should be merged by default");
    Ok(())
  }

  #[test]
  fn test_prepare_tree_after_topology_change_preserves_leaf_state() -> Result<(), Report> {
    let graph = helpers::create_polytomy_tree()?;

    // Set time_distribution on leaves (simulating load_date_constraints)
    let leaf_a_key = find_node_key_by_name(&graph, "A").ok_or_else(|| make_report!("A not found"))?;
    let leaf_b_key = find_node_key_by_name(&graph, "B").ok_or_else(|| make_report!("B not found"))?;
    let leaf_c_key = find_node_key_by_name(&graph, "C").ok_or_else(|| make_report!("C not found"))?;

    let date_dist = Arc::new(Distribution::point(2020.0, 1.0));
    for &key in &[leaf_a_key, leaf_b_key, leaf_c_key] {
      let node = graph.get_node(key).expect("Node must exist");
      node.read_arc().payload().write_arc().time_distribution = Some(Arc::clone(&date_dist));
    }

    // Mark leaf B as bad_branch (simulating outlier detection)
    {
      let node = graph.get_node(leaf_b_key).expect("Node must exist");
      node.read_arc().payload().write_arc().bad_branch = true;
    }

    // Set time_distribution on internal node ABC
    let abc_key = find_node_key_by_name(&graph, "ABC").ok_or_else(|| make_report!("ABC not found"))?;
    {
      let node = graph.get_node(abc_key).expect("Node must exist");
      node.read_arc().payload().write_arc().time_distribution = Some(Arc::new(Distribution::point(2010.0, 1.0)));
      node.read_arc().payload().write_arc().bad_branch = true;
    }

    prepare_tree_after_topology_change(&graph)?;

    // Leaf date constraints must be preserved
    for &key in &[leaf_a_key, leaf_b_key, leaf_c_key] {
      let node = graph.get_node(key).expect("Node must exist");
      let payload = node.read_arc().payload().read_arc();
      assert!(
        payload.time_distribution.is_some(),
        "Leaf time_distribution must survive topology change"
      );
    }

    // Leaf B bad_branch flag must be preserved
    {
      let node = graph.get_node(leaf_b_key).expect("Node must exist");
      let payload = node.read_arc().payload().read_arc();
      assert!(payload.bad_branch, "Leaf bad_branch flag must survive topology change");
    }

    // Internal node ABC must be cleared
    {
      let node = graph.get_node(abc_key).expect("Node must exist");
      let payload = node.read_arc().payload().read_arc();
      assert!(
        payload.time_distribution.is_none(),
        "Internal node time_distribution must be cleared"
      );
      assert!(!payload.bad_branch, "Internal node bad_branch must be cleared");
    }

    // Edge distributions must be cleared
    for edge in graph.get_edges() {
      let payload = edge.read_arc().payload().read_arc();
      assert!(
        payload.branch_length_distribution.is_none(),
        "Edge branch_length_distribution must be cleared"
      );
      assert!(payload.msg_to_parent.is_none(), "Edge msg_to_parent must be cleared");
    }

    Ok(())
  }

  mod helpers {
    use super::*;

    /// Create a tree with a polytomy (node with 3+ children).
    /// Tree: ((A,B,C)ABC)root where ABC has 3 children
    ///
    /// All children are "stretched" by default: branch_length (mutation_length) is
    /// set lower than time_length * clock_rate so the split classifies them as
    /// stretched. This matches the most common test scenario.
    pub fn create_polytomy_tree() -> Result<GraphTimetree, Report> {
      let graph: GraphTimetree = nwk_read_str("((A:0.1,B:0.2,C:0.15)ABC:0.05)root;")?;

      let tip_times = [("A", 2020.0), ("B", 2015.0), ("C", 2018.0)];
      for (name, time) in tip_times {
        let key = find_node_key_by_name(&graph, name).ok_or_else(|| make_report!("{name} not found"))?;
        let node = graph.get_node(key).expect("Node must exist");
        node.write_arc().payload().write_arc().time = Some(time);
      }

      if let Some(abc_key) = find_node_key_by_name(&graph, "ABC") {
        let node = graph.get_node(abc_key).expect("Node must exist");
        node.write_arc().payload().write_arc().time = Some(2010.0);
      }

      if let Some(root_key) = find_node_key_by_name(&graph, "root") {
        let node = graph.get_node(root_key).expect("Node must exist");
        node.write_arc().payload().write_arc().time = Some(2000.0);
      }

      for edge in graph.get_edges() {
        let edge = edge.write_arc();
        let mut payload = edge.payload().write_arc();
        payload.branch_length_distribution = Some(Arc::new(Distribution::point(0.1, 1.0)));
        // Ensure children are "stretched": mutation_length < time_length * clock_rate.
        // Newick branch lengths become time_length. With TEST_CLOCK_RATE=0.001,
        // clock_length = time_length * 0.001 which is tiny. Set mutation_length even
        // smaller to guarantee stretched.
        payload.time_length = payload.branch_length();
        payload.set_branch_length(Some(0.0));
      }

      Ok(graph)
    }

    /// Create a binary tree (no polytomies).
    pub fn create_binary_tree() -> Result<GraphTimetree, Report> {
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

    /// Create a polytomy tree with branch length distributions that produce genuine
    /// cost improvement from splitting (exponential decay: short branches more probable).
    /// All children are stretched (mutation_length=0 < clock_length).
    pub fn create_polytomy_tree_with_realistic_distributions() -> Result<GraphTimetree, Report> {
      let graph: GraphTimetree = nwk_read_str("((A:0.1,B:0.2,C:0.15)ABC:0.05)root;")?;

      let tip_times = [("A", 2020.0), ("B", 2015.0), ("C", 2018.0)];
      for (name, time) in tip_times {
        let key = find_node_key_by_name(&graph, name).ok_or_else(|| make_report!("{name} not found"))?;
        let node = graph.get_node(key).expect("Node must exist");
        node.write_arc().payload().write_arc().time = Some(time);
      }

      if let Some(key) = find_node_key_by_name(&graph, "ABC") {
        let node = graph.get_node(key).expect("Node must exist");
        node.write_arc().payload().write_arc().time = Some(2010.0);
      }

      if let Some(key) = find_node_key_by_name(&graph, "root") {
        let node = graph.get_node(key).expect("Node must exist");
        node.write_arc().payload().write_arc().time = Some(2000.0);
      }

      let x = Array1::linspace(0.0, 25.0, 200);
      let y = x.mapv(|t: f64| (-0.5 * t).exp());
      let dist = Arc::new(Distribution::function(x, y)?);

      for edge in graph.get_edges() {
        let edge = edge.write_arc();
        let mut payload = edge.payload().write_arc();
        payload.branch_length_distribution = Some(Arc::clone(&dist));
        payload.time_length = payload.branch_length();
        payload.set_branch_length(Some(0.0));
      }

      Ok(graph)
    }

    /// Create a polytomy tree where all children are "compressed"
    /// (mutation_length >= clock_length).
    pub fn create_compressed_polytomy_tree() -> Result<GraphTimetree, Report> {
      let graph: GraphTimetree = nwk_read_str("((A:0.1,B:0.2,C:0.15)ABC:0.05)root;")?;

      let tip_times = [("A", 2020.0), ("B", 2015.0), ("C", 2018.0)];
      for (name, time) in tip_times {
        let key = find_node_key_by_name(&graph, name).ok_or_else(|| make_report!("{name} not found"))?;
        let node = graph.get_node(key).expect("Node must exist");
        node.write_arc().payload().write_arc().time = Some(time);
      }

      if let Some(key) = find_node_key_by_name(&graph, "ABC") {
        let node = graph.get_node(key).expect("Node must exist");
        node.write_arc().payload().write_arc().time = Some(2010.0);
      }

      if let Some(key) = find_node_key_by_name(&graph, "root") {
        let node = graph.get_node(key).expect("Node must exist");
        node.write_arc().payload().write_arc().time = Some(2000.0);
      }

      let x = Array1::linspace(0.0, 25.0, 200);
      let y = x.mapv(|t: f64| (-0.5 * t).exp());
      let dist = Arc::new(Distribution::function(x, y)?);

      for edge in graph.get_edges() {
        let edge = edge.write_arc();
        let mut payload = edge.payload().write_arc();
        payload.branch_length_distribution = Some(Arc::clone(&dist));
        payload.time_length = payload.branch_length();
        // High mutation_length makes children "compressed"
        payload.set_branch_length(Some(1.0));
      }

      Ok(graph)
    }
  }
}
