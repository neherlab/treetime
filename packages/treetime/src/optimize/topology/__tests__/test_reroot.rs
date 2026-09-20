#[cfg(test)]
mod tests {
  use crate::test_utils::find_node_key_by_name;
  use approx::assert_abs_diff_eq;
  use eyre::Report;
  use pretty_assertions::assert_eq;
  use treetime_graph::graph::Graph;
  use treetime_graph::reroot::{
    apply_reroot_topology, record_merge, remove_node_if_trivial, split_edge, trivial_node_branch_lengths,
  };
  use treetime_io::nwk::{NwkWriteOptions, nwk_read_str, nwk_write_str};

  #[test]
  fn test_reroot_split_edge_divides_branch_length() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("(A:0.6,B:0.4)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let mut graph: Graph = graph;

    let root_key = find_node_key_by_name(&graph, &names, "root").unwrap();
    let a_key = find_node_key_by_name(&graph, &names, "A").unwrap();

    let edge_key = graph
      .get_edges()
      .find(|e| e.source() == root_key && e.target() == a_key)
      .map(|e| e.key())
      .unwrap();

    let info = split_edge(&mut graph, edge_key, 0.25, branch_lengths[&edge_key])?;

    assert_eq!(info.old_edge_key, edge_key);
    assert_abs_diff_eq!(info.split_position, 0.25, epsilon = 1e-7);

    let parent_bl = info.parent_side_length.unwrap();
    assert_abs_diff_eq!(parent_bl, 0.15, epsilon = 1e-7);

    let child_bl = info.child_side_length.unwrap();
    assert_abs_diff_eq!(child_bl, 0.45, epsilon = 1e-7);

    assert!(graph.get_node(info.new_node_key).is_some());

    assert!(graph.get_edge(edge_key).is_none());

    Ok(())
  }

  #[test]
  fn test_reroot_split_edge_at_midpoint() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("(A:1.0,B:2.0)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let mut graph: Graph = graph;

    let root_key = find_node_key_by_name(&graph, &names, "root").unwrap();
    let b_key = find_node_key_by_name(&graph, &names, "B").unwrap();

    let edge_key = graph
      .get_edges()
      .find(|e| e.source() == root_key && e.target() == b_key)
      .map(|e| e.key())
      .unwrap();

    let info = split_edge(&mut graph, edge_key, 0.5, branch_lengths[&edge_key])?;

    let parent_bl = info.parent_side_length.unwrap();
    assert_abs_diff_eq!(parent_bl, 1.0, epsilon = 1e-7);

    let child_bl = info.child_side_length.unwrap();
    assert_abs_diff_eq!(child_bl, 1.0, epsilon = 1e-7);

    Ok(())
  }

  #[test]
  fn test_reroot_apply_reroot_topology_inverts_path() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("((A:0.1,B:0.2)AB:0.3,C:0.4)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let mut graph: Graph = graph;

    let root_key = find_node_key_by_name(&graph, &names, "root").unwrap();
    let ab_key = find_node_key_by_name(&graph, &names, "AB").unwrap();

    let inverted = apply_reroot_topology(&mut graph, root_key, ab_key)?;

    assert_eq!(inverted.len(), 1);

    let ab_node = graph.get_node(ab_key).unwrap();
    assert!(ab_node.inbound().is_empty());

    let old_root = graph.get_node(root_key).unwrap();
    assert!(!old_root.inbound().is_empty());

    Ok(())
  }

  #[test]
  fn test_reroot_apply_reroot_topology_multi_hop() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("((A:0.1,B:0.2)AB:0.3,(C:0.15,D:0.25)CD:0.4)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let mut graph: Graph = graph;

    let root_key = find_node_key_by_name(&graph, &names, "root").unwrap();
    let a_key = find_node_key_by_name(&graph, &names, "A").unwrap();

    let inverted = apply_reroot_topology(&mut graph, root_key, a_key)?;

    assert_eq!(inverted.len(), 2);

    let a_node = graph.get_node(a_key).unwrap();
    assert!(a_node.inbound().is_empty());

    Ok(())
  }

  #[test]
  fn test_reroot_apply_reroot_topology_preserves_leaf_count() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("((A:0.1,B:0.2)AB:0.3,(C:0.15,D:0.25)CD:0.4)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let mut graph: Graph = graph;

    let initial_leaves = graph.get_leaves().count();
    let root_key = find_node_key_by_name(&graph, &names, "root").unwrap();
    let cd_key = find_node_key_by_name(&graph, &names, "CD").unwrap();

    apply_reroot_topology(&mut graph, root_key, cd_key)?;

    assert_eq!(graph.get_leaves().count(), initial_leaves);

    Ok(())
  }

  #[test]
  fn test_reroot_remove_node_if_trivial_merges_edges() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("((A:0.5)mid:0.3,B:0.2)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
    let mut graph: Graph = graph;

    let mid_key = find_node_key_by_name(&graph, &names, "mid").unwrap();

    let (mid_parent, mid_child) = trivial_node_branch_lengths(&graph, mid_key, &branch_lengths);
    let result = remove_node_if_trivial(&mut graph, mid_key, mid_parent, mid_child)?;

    let merge_info = result.expect("Node should have been removed");
    record_merge(&mut branch_lengths, &merge_info);
    assert_eq!(merge_info.removed_node_key, mid_key);

    assert!(graph.get_node(mid_key).is_none());

    let merged_bl = merge_info.merged_branch_length.unwrap();
    assert_abs_diff_eq!(merged_bl, 0.8, epsilon = 1e-7);

    let expected = "(B:0.2,A:0.8)root;";
    let actual = nwk_write_str(&graph, &names, &branch_lengths, &NwkWriteOptions::default())?;
    assert_eq!(expected, actual);

    Ok(())
  }

  #[test]
  fn test_reroot_remove_node_if_trivial_non_trivial_returns_none() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("((A:0.1,B:0.2)AB:0.3,C:0.4)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let mut graph: Graph = graph;

    let ab_key = find_node_key_by_name(&graph, &names, "AB").unwrap();
    let (ab_parent, ab_child) = trivial_node_branch_lengths(&graph, ab_key, &branch_lengths);
    let result = remove_node_if_trivial(&mut graph, ab_key, ab_parent, ab_child)?;
    assert!(result.is_none());

    let root_key = find_node_key_by_name(&graph, &names, "root").unwrap();
    let (root_parent, root_child) = trivial_node_branch_lengths(&graph, root_key, &branch_lengths);
    let result = remove_node_if_trivial(&mut graph, root_key, root_parent, root_child)?;
    assert!(result.is_none());

    Ok(())
  }

  #[test]
  fn test_reroot_full_reroot_and_cleanup_preserves_topology() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("((A:0.1,B:0.2)AB:0.3,(C:0.15,D:0.25)CD:0.4)root:0.001;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
    let mut graph: Graph = graph;

    let root_key = find_node_key_by_name(&graph, &names, "root").unwrap();
    let cd_key = find_node_key_by_name(&graph, &names, "CD").unwrap();

    apply_reroot_topology(&mut graph, root_key, cd_key)?;
    let (root_parent, root_child) = trivial_node_branch_lengths(&graph, root_key, &branch_lengths);
    if let Some(info) = remove_node_if_trivial(&mut graph, root_key, root_parent, root_child)? {
      record_merge(&mut branch_lengths, &info);
    }

    let cd_node = graph.get_node(cd_key).unwrap();
    assert!(cd_node.inbound().is_empty());

    assert_eq!(graph.get_leaves().count(), 4);

    let newick = nwk_write_str(
      &graph,
      &names,
      &branch_lengths,
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
    let nwk_parsed = nwk_read_str("((A:0.5)mid,B:0.2)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let mut graph: Graph = graph;

    let mid_key = find_node_key_by_name(&graph, &names, "mid").unwrap();
    let (mid_parent, mid_child) = trivial_node_branch_lengths(&graph, mid_key, &branch_lengths);
    let result = remove_node_if_trivial(&mut graph, mid_key, mid_parent, mid_child)?;

    let merge_info = result.expect("Trivial node should be removed");

    let merged_bl = merge_info.merged_branch_length;
    assert!(merged_bl.is_some());

    Ok(())
  }

  #[test]
  fn test_reroot_full_cycle_branch_length_conservation() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("((A:0.1,B:0.2)AB:0.3,(C:0.15,D:0.25)CD:0.4)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
    let mut graph: Graph = graph;

    let total_bl_before: f64 = graph
      .get_edges()
      .filter_map(|e| branch_lengths.get(&e.key()).copied().flatten())
      .sum();

    let root_key = find_node_key_by_name(&graph, &names, "root").unwrap();
    let ab_key = find_node_key_by_name(&graph, &names, "AB").unwrap();

    apply_reroot_topology(&mut graph, root_key, ab_key)?;
    let (root_parent, root_child) = trivial_node_branch_lengths(&graph, root_key, &branch_lengths);
    if let Some(info) = remove_node_if_trivial(&mut graph, root_key, root_parent, root_child)? {
      record_merge(&mut branch_lengths, &info);
    }

    let total_bl_after: f64 = graph
      .get_edges()
      .filter_map(|e| branch_lengths.get(&e.key()).copied().flatten())
      .sum();

    assert_abs_diff_eq!(total_bl_before, total_bl_after, epsilon = 1e-7);

    Ok(())
  }
}
