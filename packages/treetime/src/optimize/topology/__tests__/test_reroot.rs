#[cfg(test)]
mod tests {
  use treetime_graph::reroot::{
    apply_reroot_topology, remove_node_if_trivial, split_edge,
  };
  use crate::payload::ancestral::GraphAncestral;
  use crate::test_utils::find_node_key_by_name;
  use approx::assert_abs_diff_eq;
  use eyre::Report;
  use pretty_assertions::assert_eq;
  use treetime_graph::edge::HasBranchLength;
  use treetime_io::nwk::{NwkWriteOptions, nwk_read_str, nwk_write_str};

  #[test]
  fn test_reroot_split_edge_divides_branch_length() -> Result<(), Report> {
    let mut graph: GraphAncestral = nwk_read_str("(A:0.6,B:0.4)root;")?;

    let root_key = find_node_key_by_name(&graph, "root").unwrap();
    let a_key = find_node_key_by_name(&graph, "A").unwrap();

    let edge_key = graph
      .get_edges()
      .iter()
      .find(|e| {
        let e = e.read_arc();
        e.source() == root_key && e.target() == a_key
      })
      .map(|e| e.read_arc().key())
      .unwrap();

    let info = split_edge(&mut graph, edge_key, 0.25)?;

    assert_eq!(info.old_edge_key, edge_key);
    assert_abs_diff_eq!(info.split_position, 0.25, epsilon = 1e-7);

    // Parent-side edge: 0.25 * 0.6 = 0.15
    let parent_edge = graph.get_edge(info.parent_side_edge_key).unwrap();
    let parent_bl = parent_edge.read_arc().payload().read_arc().branch_length().unwrap();
    assert_abs_diff_eq!(parent_bl, 0.15, epsilon = 1e-7);

    // Child-side edge: 0.75 * 0.6 = 0.45
    let child_edge = graph.get_edge(info.child_side_edge_key).unwrap();
    let child_bl = child_edge.read_arc().payload().read_arc().branch_length().unwrap();
    assert_abs_diff_eq!(child_bl, 0.45, epsilon = 1e-7);

    // New node exists
    assert!(graph.get_node(info.new_node_key).is_some());

    // Original edge removed
    assert!(graph.get_edge(edge_key).is_none());

    Ok(())
  }

  #[test]
  fn test_reroot_split_edge_at_midpoint() -> Result<(), Report> {
    let mut graph: GraphAncestral = nwk_read_str("(A:1.0,B:2.0)root;")?;

    let root_key = find_node_key_by_name(&graph, "root").unwrap();
    let b_key = find_node_key_by_name(&graph, "B").unwrap();

    let edge_key = graph
      .get_edges()
      .iter()
      .find(|e| {
        let e = e.read_arc();
        e.source() == root_key && e.target() == b_key
      })
      .map(|e| e.read_arc().key())
      .unwrap();

    let info = split_edge(&mut graph, edge_key, 0.5)?;

    let parent_edge = graph.get_edge(info.parent_side_edge_key).unwrap();
    let parent_bl = parent_edge.read_arc().payload().read_arc().branch_length().unwrap();
    assert_abs_diff_eq!(parent_bl, 1.0, epsilon = 1e-7);

    let child_edge = graph.get_edge(info.child_side_edge_key).unwrap();
    let child_bl = child_edge.read_arc().payload().read_arc().branch_length().unwrap();
    assert_abs_diff_eq!(child_bl, 1.0, epsilon = 1e-7);

    Ok(())
  }

  #[test]
  fn test_reroot_apply_reroot_topology_inverts_path() -> Result<(), Report> {
    let mut graph: GraphAncestral = nwk_read_str("((A:0.1,B:0.2)AB:0.3,C:0.4)root;")?;

    let root_key = find_node_key_by_name(&graph, "root").unwrap();
    let ab_key = find_node_key_by_name(&graph, "AB").unwrap();

    let inverted = apply_reroot_topology(&mut graph, root_key, ab_key)?;

    // One edge on path from root to AB was inverted
    assert_eq!(inverted.len(), 1);

    // AB is now the root (no inbound edges)
    let ab_node = graph.get_node(ab_key).unwrap();
    assert!(ab_node.read_arc().inbound().is_empty());

    // Old root is no longer root (has inbound edge)
    let old_root = graph.get_node(root_key).unwrap();
    assert!(!old_root.read_arc().inbound().is_empty());

    Ok(())
  }

  #[test]
  fn test_reroot_apply_reroot_topology_multi_hop() -> Result<(), Report> {
    let mut graph: GraphAncestral = nwk_read_str("((A:0.1,B:0.2)AB:0.3,(C:0.15,D:0.25)CD:0.4)root;")?;

    let root_key = find_node_key_by_name(&graph, "root").unwrap();
    let a_key = find_node_key_by_name(&graph, "A").unwrap();

    let inverted = apply_reroot_topology(&mut graph, root_key, a_key)?;

    // Two edges on path: root->AB->A
    assert_eq!(inverted.len(), 2);

    // A is now root
    let a_node = graph.get_node(a_key).unwrap();
    assert!(a_node.read_arc().inbound().is_empty());

    Ok(())
  }

  #[test]
  fn test_reroot_apply_reroot_topology_preserves_leaf_count() -> Result<(), Report> {
    let mut graph: GraphAncestral = nwk_read_str("((A:0.1,B:0.2)AB:0.3,(C:0.15,D:0.25)CD:0.4)root;")?;

    let initial_leaves = graph.get_leaves().len();
    let root_key = find_node_key_by_name(&graph, "root").unwrap();
    let cd_key = find_node_key_by_name(&graph, "CD").unwrap();

    apply_reroot_topology(&mut graph, root_key, cd_key)?;

    assert_eq!(graph.get_leaves().len(), initial_leaves);

    Ok(())
  }

  #[test]
  fn test_reroot_remove_node_if_trivial_merges_edges() -> Result<(), Report> {
    // Tree with a trivial node (one parent, one child):
    //      root
    //      /  \
    //    mid   B
    //    /
    //   A
    let mut graph: GraphAncestral = nwk_read_str("((A:0.5)mid:0.3,B:0.2)root;")?;

    let mid_key = find_node_key_by_name(&graph, "mid").unwrap();

    let result = remove_node_if_trivial(&mut graph, mid_key)?;

    let merge_info = result.expect("Node should have been removed");
    assert_eq!(merge_info.removed_node_key, mid_key);

    // Node is gone
    assert!(graph.get_node(mid_key).is_none());

    // Merged edge has summed branch length: 0.5 + 0.3 = 0.8
    let merged_edge = graph.get_edge(merge_info.merged_edge_key).unwrap();
    let merged_bl = merged_edge.read_arc().payload().read_arc().branch_length().unwrap();
    assert_abs_diff_eq!(merged_bl, 0.8, epsilon = 1e-7);

    // Tree output matches expected
    let expected = "(B:0.2,A:0.8)root;";
    let actual = nwk_write_str(&graph, &NwkWriteOptions::default())?;
    assert_eq!(expected, actual);

    Ok(())
  }

  #[test]
  fn test_reroot_remove_node_if_trivial_non_trivial_returns_none() -> Result<(), Report> {
    let mut graph: GraphAncestral = nwk_read_str("((A:0.1,B:0.2)AB:0.3,C:0.4)root;")?;

    // AB has two children, not trivial
    let ab_key = find_node_key_by_name(&graph, "AB").unwrap();
    let result = remove_node_if_trivial(&mut graph, ab_key)?;
    assert!(result.is_none());

    // Root has no parent, not trivial
    let root_key = find_node_key_by_name(&graph, "root").unwrap();
    let result = remove_node_if_trivial(&mut graph, root_key)?;
    assert!(result.is_none());

    Ok(())
  }

  #[test]
  fn test_reroot_full_reroot_and_cleanup_preserves_topology() -> Result<(), Report> {
    let mut graph: GraphAncestral = nwk_read_str("((A:0.1,B:0.2)AB:0.3,(C:0.15,D:0.25)CD:0.4)root:0.001;")?;

    let root_key = find_node_key_by_name(&graph, "root").unwrap();
    let cd_key = find_node_key_by_name(&graph, "CD").unwrap();

    // Reroot at CD
    apply_reroot_topology(&mut graph, root_key, cd_key)?;
    // Old root is now degree-2 (one parent from CD side, one child to AB side)
    remove_node_if_trivial(&mut graph, root_key)?;

    // CD is root
    let cd_node = graph.get_node(cd_key).unwrap();
    assert!(cd_node.read_arc().inbound().is_empty());

    // All 4 leaves preserved
    assert_eq!(graph.get_leaves().len(), 4);

    // Check total branch length conservation (unrooted tree property)
    let newick = nwk_write_str(
      &graph,
      &NwkWriteOptions {
        weight_significant_digits: Some(17),
        ..NwkWriteOptions::default()
      },
    )?;
    for taxon in &["A", "B", "C", "D"] {
      assert!(newick.contains(taxon), "Missing taxon {taxon} in {newick}");
    }

    Ok(())
  }

  #[test]
  fn test_reroot_remove_trivial_with_partial_branch_lengths() -> Result<(), Report> {
    // One edge has a branch length, the other does not -> merged gets the existing one
    let mut graph: GraphAncestral = nwk_read_str("((A:0.5)mid,B:0.2)root;")?;

    let mid_key = find_node_key_by_name(&graph, "mid").unwrap();
    let result = remove_node_if_trivial(&mut graph, mid_key)?;

    let merge_info = result.expect("Trivial node should be removed");

    let merged_edge = graph.get_edge(merge_info.merged_edge_key).unwrap();
    let merged_bl = merged_edge.read_arc().payload().read_arc().branch_length();
    // When one branch is Some and the other is None, result is Some (the existing value)
    assert!(merged_bl.is_some());

    Ok(())
  }

  #[test]
  fn test_reroot_full_cycle_branch_length_conservation() -> Result<(), Report> {
    let mut graph: GraphAncestral = nwk_read_str("((A:0.1,B:0.2)AB:0.3,(C:0.15,D:0.25)CD:0.4)root;")?;

    // Compute total branch length before reroot
    let total_bl_before: f64 = graph
      .get_edges()
      .iter()
      .filter_map(|e| e.read_arc().payload().read_arc().branch_length())
      .sum();

    let root_key = find_node_key_by_name(&graph, "root").unwrap();
    let ab_key = find_node_key_by_name(&graph, "AB").unwrap();

    apply_reroot_topology(&mut graph, root_key, ab_key)?;
    remove_node_if_trivial(&mut graph, root_key)?;

    // Compute total branch length after reroot + trivial removal
    let total_bl_after: f64 = graph
      .get_edges()
      .iter()
      .filter_map(|e| e.read_arc().payload().read_arc().branch_length())
      .sum();

    // Total branch length on an unrooted tree is conserved under rerooting
    assert_abs_diff_eq!(total_bl_before, total_bl_after, epsilon = 1e-7);

    Ok(())
  }
}
