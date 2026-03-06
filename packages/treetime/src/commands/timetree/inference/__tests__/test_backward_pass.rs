#[cfg(test)]
mod tests {
  use crate::commands::timetree::inference::backward_pass::propagate_distributions_backward;
  use crate::representation::payload::timetree::{EdgeTimetree, NodeTimetree};
  use crate::test_utils::find_node_key_by_name;
  use approx::assert_ulps_eq;
  use eyre::Report;
  use std::sync::Arc;
  use treetime_distribution::Distribution;
  use treetime_graph::edge::BranchDistribution;
  use treetime_graph::node::{GraphNodeKey, TimeConstraint};
  use treetime_io::nwk::nwk_read_str;

  /// Test backward pass on a simple 2-leaf tree:
  /// ((A:2.5)I:1.0)root;
  /// A has time 2013.0, I should get 2013.0 - 2.5 = 2010.5
  #[test]
  fn test_backward_pass_computes_internal_node_time() -> Result<(), Report> {
    let graph = nwk_read_str::<NodeTimetree, EdgeTimetree, ()>("((A:2.5)I:1.0)root;")?;

    // Set leaf A's time distribution to point at 2013.0
    let leaf_key = find_node_key_by_name(&graph, "A").expect("leaf A not found");
    {
      let leaf_node = graph.get_node(leaf_key).expect("leaf A exists");
      let mut payload = leaf_node.read_arc().payload().write_arc();
      payload.time_distribution = Some(Arc::new(Distribution::point(2013.0, 1.0)));
    }

    // Set branch length distribution on edge from I to A (branch length 2.5 years)
    for edge in graph.get_edges() {
      let edge_read = edge.read_arc();
      if edge_read.target() == leaf_key {
        let mut payload = edge_read.payload().write_arc();
        payload.set_branch_length_distribution(Some(Arc::new(Distribution::point(2.5, 1.0))));
      }
    }

    // Run backward pass
    propagate_distributions_backward(&graph, None)?;

    // Check internal node I has time distribution centered at 2013.0 - 2.5 = 2010.5
    let internal_key = find_node_key_by_name(&graph, "I").expect("internal node I not found");
    let internal_node = graph.get_node(internal_key).expect("internal node I exists");
    let payload = internal_node.read_arc().payload().read_arc();

    let time_dist = payload
      .time_distribution()
      .as_ref()
      .expect("internal node should have time distribution after backward pass");
    let likely_time = time_dist
      .likely_time()
      .expect("time distribution should have likely_time");

    assert_ulps_eq!(likely_time, 2010.5, max_ulps = 4);
    assert!(likely_time < 2013.0, "Parent should be older than child");

    Ok(())
  }

  /// Test backward pass with two children: messages are multiplied (intersection).
  /// Tree: ((A:3.0,B:2.0)I:1.0)root;
  /// A at 2015.0, B at 2014.0
  /// I gets messages: from A -> 2015-3=2012, from B -> 2014-2=2012
  /// Both agree, so I should be at 2012.0
  #[test]
  fn test_backward_pass_multiplies_child_messages() -> Result<(), Report> {
    let graph = nwk_read_str::<NodeTimetree, EdgeTimetree, ()>("((A:3.0,B:2.0)I:1.0)root;")?;

    let leaf_a_key = find_node_key_by_name(&graph, "A").expect("leaf A not found");
    let leaf_b_key = find_node_key_by_name(&graph, "B").expect("leaf B not found");

    // Set time distributions on leaves
    {
      let node_a = graph.get_node(leaf_a_key).expect("leaf A exists");
      let mut payload = node_a.read_arc().payload().write_arc();
      payload.time_distribution = Some(Arc::new(Distribution::point(2015.0, 1.0)));
    }
    {
      let node_b = graph.get_node(leaf_b_key).expect("leaf B exists");
      let mut payload = node_b.read_arc().payload().write_arc();
      payload.time_distribution = Some(Arc::new(Distribution::point(2014.0, 1.0)));
    }

    // Set branch length distributions
    for edge in graph.get_edges() {
      let edge_read = edge.read_arc();
      let target = edge_read.target();
      let branch_length = if target == leaf_a_key {
        3.0
      } else if target == leaf_b_key {
        2.0
      } else {
        continue;
      };
      let mut payload = edge_read.payload().write_arc();
      payload.set_branch_length_distribution(Some(Arc::new(Distribution::point(branch_length, 1.0))));
    }

    propagate_distributions_backward(&graph, None)?;

    // Internal node I should have time at 2012.0 (both children agree)
    let internal_key = find_node_key_by_name(&graph, "I").expect("internal node I not found");
    let internal_node = graph.get_node(internal_key).expect("internal node I exists");
    let payload = internal_node.read_arc().payload().read_arc();

    let time_dist = payload
      .time_distribution()
      .as_ref()
      .expect("internal node should have time distribution");
    let likely_time = time_dist
      .likely_time()
      .expect("time distribution should have likely_time");

    assert_ulps_eq!(likely_time, 2012.0, max_ulps = 4);

    Ok(())
  }

  /// Test that backward pass preserves leaf time_distribution when coalescent is provided.
  /// This is a regression test for a bug where coalescent contributions overwrote leaf dates,
  /// causing subsequent clock regression to fail with "No variation in sampling dates".
  #[test]
  fn test_backward_pass_preserves_leaf_time_distribution_with_coalescent() -> Result<(), Report> {
    use indexmap::IndexMap;
    use treetime_distribution::DistributionNegLog;

    let graph = nwk_read_str::<NodeTimetree, EdgeTimetree, ()>("((A:3.0,B:2.0)I:1.0)root;")?;

    let leaf_a_key = find_node_key_by_name(&graph, "A").expect("leaf A not found");
    let leaf_b_key = find_node_key_by_name(&graph, "B").expect("leaf B not found");

    let date_a = 2015.0;
    let date_b = 2014.0;

    // Set time distributions on leaves (date constraints)
    {
      let node_a = graph.get_node(leaf_a_key).expect("leaf A exists");
      let mut payload = node_a.read_arc().payload().write_arc();
      payload.time_distribution = Some(Arc::new(Distribution::point(date_a, 1.0)));
    }
    {
      let node_b = graph.get_node(leaf_b_key).expect("leaf B exists");
      let mut payload = node_b.read_arc().payload().write_arc();
      payload.time_distribution = Some(Arc::new(Distribution::point(date_b, 1.0)));
    }

    // Set branch length distributions
    for edge in graph.get_edges() {
      let edge_read = edge.read_arc();
      let target = edge_read.target();
      let branch_length = if target == leaf_a_key {
        3.0
      } else if target == leaf_b_key {
        2.0
      } else {
        continue;
      };
      let mut payload = edge_read.payload().write_arc();
      payload.set_branch_length_distribution(Some(Arc::new(Distribution::point(branch_length, 1.0))));
    }

    // Create coalescent contributions for all nodes (simulating --coalescent option)
    let mut coalescent_contributions: IndexMap<_, Arc<DistributionNegLog>> = IndexMap::new();
    for node in graph.get_nodes() {
      let key = node.read_arc().key();
      // Use a constant coalescent contribution
      coalescent_contributions.insert(key, Arc::new(DistributionNegLog::constant(0.01)));
    }

    // Run backward pass WITH coalescent contributions
    propagate_distributions_backward(&graph, Some(&coalescent_contributions))?;

    // Verify leaf A still has its original date
    {
      let node_a = graph.get_node(leaf_a_key).expect("leaf A exists");
      let payload = node_a.read_arc().payload().read_arc();
      let time_dist = payload
        .time_distribution()
        .as_ref()
        .expect("leaf A should have time distribution");
      let likely_time = time_dist
        .likely_time()
        .expect("leaf A time distribution should have likely_time");
      assert_ulps_eq!(likely_time, date_a, max_ulps = 4);
    }

    // Verify leaf B still has its original date
    {
      let node_b = graph.get_node(leaf_b_key).expect("leaf B exists");
      let payload = node_b.read_arc().payload().read_arc();
      let time_dist = payload
        .time_distribution()
        .as_ref()
        .expect("leaf B should have time distribution");
      let likely_time = time_dist
        .likely_time()
        .expect("leaf B time distribution should have likely_time");
      assert_ulps_eq!(likely_time, date_b, max_ulps = 4);
    }

    Ok(())
  }

  /// Test that backward pass stores msg_to_parent on edges.
  #[test]
  fn test_backward_pass_sets_edge_messages() -> Result<(), Report> {
    let graph = nwk_read_str::<NodeTimetree, EdgeTimetree, ()>("((A:2.5)I:1.0)root;")?;

    let leaf_key = find_node_key_by_name(&graph, "A").expect("leaf A not found");

    // Set leaf time distribution
    {
      let leaf_node = graph.get_node(leaf_key).expect("leaf A exists");
      let mut payload = leaf_node.read_arc().payload().write_arc();
      payload.time_distribution = Some(Arc::new(Distribution::point(2013.0, 1.0)));
    }

    // Set branch length distribution
    for edge in graph.get_edges() {
      let edge_read = edge.read_arc();
      if edge_read.target() == leaf_key {
        let mut payload = edge_read.payload().write_arc();
        payload.set_branch_length_distribution(Some(Arc::new(Distribution::point(2.5, 1.0))));
      }
    }

    propagate_distributions_backward(&graph, None)?;

    // Check edge from I to A has msg_to_parent set
    for edge in graph.get_edges() {
      let edge_read = edge.read_arc();
      if edge_read.target() == leaf_key {
        let payload = edge_read.payload().read_arc();
        let msg = payload
          .msg_to_parent()
          .as_ref()
          .expect("edge should have msg_to_parent after backward pass");
        let msg_time = msg.likely_time().expect("message should have likely_time");
        // Message should be parent time: 2013.0 - 2.5 = 2010.5
        assert_ulps_eq!(msg_time, 2010.5, max_ulps = 4);
      }
    }

    Ok(())
  }

  /// Test backward pass skips children with bad_branch=true.
  /// Tree: ((A:3.0,B:2.0)I:1.0)root;
  /// A at 2015.0, B at 2014.0. Mark B as bad_branch.
  /// I should get time from A only: 2015.0 - 3.0 = 2012.0
  #[test]
  fn test_backward_pass_skips_bad_branch_children() -> Result<(), Report> {
    let graph = nwk_read_str::<NodeTimetree, EdgeTimetree, ()>("((A:3.0,B:2.0)I:1.0)root;")?;

    let leaf_a_key = find_node_key_by_name(&graph, "A").expect("leaf A not found");
    let leaf_b_key = find_node_key_by_name(&graph, "B").expect("leaf B not found");

    // Set time distributions on leaves
    set_leaf_time(&graph, leaf_a_key, 2015.0);
    set_leaf_time(&graph, leaf_b_key, 2014.0);

    // Mark B as bad_branch
    {
      let node_b = graph.get_node(leaf_b_key).expect("leaf B exists");
      node_b.read_arc().payload().write_arc().bad_branch = true;
    }

    // Set branch length distributions on both edges
    set_edge_branch_dist(&graph, leaf_a_key, 3.0);
    set_edge_branch_dist(&graph, leaf_b_key, 2.0);

    propagate_distributions_backward(&graph, None)?;

    // I should get time only from A: 2015.0 - 3.0 = 2012.0
    let internal_key = find_node_key_by_name(&graph, "I").expect("internal node I not found");
    let internal_node = graph.get_node(internal_key).expect("internal node I exists");
    let payload = internal_node.read_arc().payload().read_arc();

    let time_dist = payload
      .time_distribution()
      .as_ref()
      .expect("internal node should have time distribution");
    let likely_time = time_dist
      .likely_time()
      .expect("time distribution should have likely_time");

    assert_ulps_eq!(likely_time, 2012.0, max_ulps = 4);

    Ok(())
  }

  /// Test that marking a leaf as bad_branch produces identical parent time
  /// to a tree where that leaf is absent.
  #[test]
  fn test_backward_pass_bad_branch_equivalent_to_removal() -> Result<(), Report> {
    // Reference tree: only A
    let ref_graph = nwk_read_str::<NodeTimetree, EdgeTimetree, ()>("((A:3.0)I:1.0)root;")?;
    let ref_a_key = find_node_key_by_name(&ref_graph, "A").expect("leaf A not found");
    set_leaf_time(&ref_graph, ref_a_key, 2015.0);
    set_edge_branch_dist(&ref_graph, ref_a_key, 3.0);
    propagate_distributions_backward(&ref_graph, None)?;

    let ref_internal_key = find_node_key_by_name(&ref_graph, "I").expect("internal I not found");
    let ref_time = ref_graph
      .get_node(ref_internal_key)
      .expect("I exists")
      .read_arc()
      .payload()
      .read_arc()
      .time_distribution()
      .as_ref()
      .expect("should have time dist")
      .likely_time()
      .expect("should have likely_time");

    // Test tree: A + B(bad)
    let test_graph = nwk_read_str::<NodeTimetree, EdgeTimetree, ()>("((A:3.0,B:2.0)I:1.0)root;")?;
    let test_a_key = find_node_key_by_name(&test_graph, "A").expect("leaf A not found");
    let test_b_key = find_node_key_by_name(&test_graph, "B").expect("leaf B not found");
    set_leaf_time(&test_graph, test_a_key, 2015.0);
    set_leaf_time(&test_graph, test_b_key, 2014.0);
    set_edge_branch_dist(&test_graph, test_a_key, 3.0);
    set_edge_branch_dist(&test_graph, test_b_key, 2.0);

    // Mark B as bad
    test_graph
      .get_node(test_b_key)
      .expect("B exists")
      .read_arc()
      .payload()
      .write_arc()
      .bad_branch = true;

    propagate_distributions_backward(&test_graph, None)?;

    let test_internal_key = find_node_key_by_name(&test_graph, "I").expect("internal I not found");
    let test_time = test_graph
      .get_node(test_internal_key)
      .expect("I exists")
      .read_arc()
      .payload()
      .read_arc()
      .time_distribution()
      .as_ref()
      .expect("should have time dist")
      .likely_time()
      .expect("should have likely_time");

    assert_ulps_eq!(ref_time, test_time, max_ulps = 4);

    Ok(())
  }

  mod helpers {
    use super::*;
    use treetime_graph::graph::Graph;

    pub(super) fn set_leaf_time(graph: &Graph<NodeTimetree, EdgeTimetree, ()>, key: GraphNodeKey, time: f64) {
      let node = graph.get_node(key).expect("node exists");
      node.read_arc().payload().write_arc().time_distribution = Some(Arc::new(Distribution::point(time, 1.0)));
    }

    pub(super) fn set_edge_branch_dist(graph: &Graph<NodeTimetree, EdgeTimetree, ()>, target_key: GraphNodeKey, bl: f64) {
      for edge in graph.get_edges() {
        let edge_read = edge.read_arc();
        if edge_read.target() == target_key {
          edge_read
            .payload()
            .write_arc()
            .set_branch_length_distribution(Some(Arc::new(Distribution::point(bl, 1.0))));
        }
      }
    }
  }

  use helpers::*;
}
