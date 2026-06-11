#[cfg(test)]
mod tests {
  use crate::test_utils::{TestEdge, TestNode};
  use approx::assert_abs_diff_eq;
  use eyre::Report;
  use pretty_assertions::assert_eq;
  use std::collections::BTreeMap;
  use std::path::PathBuf;
  use treetime_graph::edge::HasBranchLength;
  use treetime_graph::graph::Graph;
  use treetime_graph::node::Named;
  use treetime_io::nwk::{NwkWriteOptions, nwk_read_file, nwk_read_str, nwk_write_str};

  #[test]
  fn test_nwk_roundtrip_binary_tree() -> Result<(), Report> {
    let input = "((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root;";
    let graph = nwk_read_str::<TestNode, TestEdge, ()>(input)?;

    // Verify graph structure
    assert_eq!(graph.get_nodes().len(), 7, "Should have 7 nodes");
    assert_eq!(graph.get_edges().len(), 6, "Should have 6 edges");
    assert_eq!(graph.get_leaves().len(), 4, "Should have 4 leaves");

    // Verify root
    let root = graph.get_exactly_one_root()?;
    let root_name = root
      .read_arc()
      .payload()
      .read_arc()
      .name()
      .map(|n| n.as_ref().to_owned());
    assert_eq!(root_name.as_deref(), Some("root"));

    // Verify roundtrip
    let output = nwk_write_str(&graph, &NwkWriteOptions::default())?;
    assert_eq!(input, output);
    Ok(())
  }

  #[test]
  fn test_nwk_parse_no_branch_lengths() -> Result<(), Report> {
    let input = "((A,B)AB,(C,D)CD)root;";
    let graph = nwk_read_str::<TestNode, TestEdge, ()>(input)?;

    // Verify structure is correct
    assert_eq!(graph.get_nodes().len(), 7);
    assert_eq!(graph.get_edges().len(), 6);
    assert_eq!(graph.get_leaves().len(), 4);

    for edge in graph.get_edges() {
      let branch_len = edge.read_arc().payload().read_arc().branch_length();
      assert!(branch_len.is_none(), "Missing branch length should be parsed as None");
    }

    Ok(())
  }

  #[test]
  fn test_nwk_roundtrip_single_node() -> Result<(), Report> {
    let input = "A;";
    let graph = nwk_read_str::<TestNode, TestEdge, ()>(input)?;

    assert_eq!(graph.get_nodes().len(), 1, "Should have 1 node");
    assert_eq!(graph.get_edges().len(), 0, "Should have 0 edges");
    assert_eq!(graph.get_leaves().len(), 1, "Should have 1 leaf");

    let root = graph.get_exactly_one_root()?;
    let root_name = root
      .read_arc()
      .payload()
      .read_arc()
      .name()
      .map(|n| n.as_ref().to_owned());
    assert_eq!(root_name.as_deref(), Some("A"));

    let output = nwk_write_str(&graph, &NwkWriteOptions::default())?;
    assert_eq!(input, output);
    Ok(())
  }

  #[test]
  fn test_nwk_roundtrip_polytomy() -> Result<(), Report> {
    // Tree with polytomy: root has 3 children
    let input = "(A:0.1,B:0.2,C:0.3)root;";
    let graph = nwk_read_str::<TestNode, TestEdge, ()>(input)?;

    assert_eq!(graph.get_nodes().len(), 4, "Should have 4 nodes");
    assert_eq!(graph.get_edges().len(), 3, "Should have 3 edges");
    assert_eq!(graph.get_leaves().len(), 3, "Should have 3 leaves");

    // Verify root has 3 children (polytomy)
    let root = graph.get_exactly_one_root()?;
    let root_outbound = root.read_arc().outbound().len();
    assert_eq!(root_outbound, 3, "Root should have 3 children (polytomy)");

    let output = nwk_write_str(&graph, &NwkWriteOptions::default())?;
    assert_eq!(input, output);
    Ok(())
  }

  #[test]
  fn test_nwk_roundtrip_nested_polytomy() -> Result<(), Report> {
    // Nested polytomies: internal node also has >2 children
    let input = "((A:0.1,B:0.2,C:0.3)ABC:0.4,D:0.5,E:0.6)root;";
    let graph = nwk_read_str::<TestNode, TestEdge, ()>(input)?;

    assert_eq!(graph.get_nodes().len(), 7, "Should have 7 nodes");
    assert_eq!(graph.get_edges().len(), 6, "Should have 6 edges");
    assert_eq!(graph.get_leaves().len(), 5, "Should have 5 leaves");

    let output = nwk_write_str(&graph, &NwkWriteOptions::default())?;
    assert_eq!(input, output);
    Ok(())
  }

  #[test]
  fn test_nwk_roundtrip_zero_length_branches() -> Result<(), Report> {
    let input = "((A:0,B:0.2)AB:0,(C:0.2,D:0)CD:0.05)root;";
    let graph = nwk_read_str::<TestNode, TestEdge, ()>(input)?;

    assert_eq!(graph.get_nodes().len(), 7);
    assert_eq!(graph.get_edges().len(), 6);

    // Count zero-length branches
    let zero_branches: usize = graph
      .get_edges()
      .iter()
      .filter(|e| e.read_arc().payload().read_arc().branch_length() == Some(0.0))
      .count();
    assert_eq!(zero_branches, 3, "Should have 3 zero-length branches");

    let output = nwk_write_str(&graph, &NwkWriteOptions::default())?;
    assert_eq!(input, output);
    Ok(())
  }

  #[test]
  fn test_nwk_parse_verifies_branch_lengths() -> Result<(), Report> {
    let input = "((A:0.123,B:0.456)AB:0.789,C:1.5)root;";
    let graph = nwk_read_str::<TestNode, TestEdge, ()>(input)?;

    // Collect branch lengths by finding edges to specific nodes
    let mut branch_lengths = BTreeMap::new();
    for edge in graph.get_edges() {
      let edge_ref = edge.read_arc();
      let target_key = edge_ref.target();
      let target_node = graph.get_node(target_key).expect("Target node should exist");
      let target_name = target_node
        .read_arc()
        .payload()
        .read_arc()
        .name()
        .map(|n| n.as_ref().to_owned());
      if let Some(name) = target_name {
        let length = edge_ref.payload().read_arc().branch_length();
        branch_lengths.insert(name, length);
      }
    }

    assert_abs_diff_eq!(
      branch_lengths["A"].expect("A should have branch length"),
      0.123,
      epsilon = 1e-10
    );
    assert_abs_diff_eq!(
      branch_lengths["B"].expect("B should have branch length"),
      0.456,
      epsilon = 1e-10
    );
    assert_abs_diff_eq!(
      branch_lengths["AB"].expect("AB should have branch length"),
      0.789,
      epsilon = 1e-10
    );
    assert_abs_diff_eq!(
      branch_lengths["C"].expect("C should have branch length"),
      1.5,
      epsilon = 1e-10
    );

    Ok(())
  }

  #[test]
  fn test_nwk_parse_verifies_leaf_names() -> Result<(), Report> {
    let input = "((leaf_A:0.1,leaf_B:0.2)internal:0.1,leaf_C:0.3)root;";
    let graph = nwk_read_str::<TestNode, TestEdge, ()>(input)?;

    let leaf_names: Vec<String> = graph
      .get_leaves()
      .iter()
      .map(|n| {
        n.read_arc()
          .payload()
          .read_arc()
          .name()
          .map(|s| s.as_ref().to_owned())
          .unwrap_or_default()
      })
      .collect();

    assert!(leaf_names.contains(&"leaf_A".to_owned()));
    assert!(leaf_names.contains(&"leaf_B".to_owned()));
    assert!(leaf_names.contains(&"leaf_C".to_owned()));
    assert_eq!(leaf_names.len(), 3);

    Ok(())
  }

  #[test]
  fn test_nwk_rejects_enewick_hybrid() {
    let result = nwk_read_str::<TestNode, TestEdge, ()>("((A,(B)x#H1)c,(x#H1,C)d);");
    let err = result.unwrap_err();
    let msg = format!("{err:?}");
    assert!(
      msg.contains("eNewick") || msg.contains("reticulate") || msg.contains("hybrid"),
      "Error should mention eNewick/hybrid/reticulate, got: {msg}"
    );
  }

  #[test]
  fn test_nwk_parse_error_includes_position() {
    let result = nwk_read_str::<TestNode, TestEdge, ()>("(A:0.1,B:0.2");
    let err = result.unwrap_err();
    let msg = format!("{err:?}");
    assert!(
      msg.contains("line") || msg.contains("-->"),
      "Parse error should include position info, got: {msg}"
    );
  }

  #[test]
  fn test_nwk_parse_with_beast_comment_succeeds() -> Result<(), Report> {
    let input = "(A[&prob=0.95]:0.1,B:0.2);";
    let graph = nwk_read_str::<TestNode, TestEdge, ()>(input)?;

    assert_eq!(graph.get_nodes().len(), 3);
    assert_eq!(graph.get_edges().len(), 2);

    let mut branch_lengths = BTreeMap::new();
    for edge in graph.get_edges() {
      let edge_ref = edge.read_arc();
      let target_node = graph.get_node(edge_ref.target()).expect("target node");
      let name = target_node
        .read_arc()
        .payload()
        .read_arc()
        .name()
        .map(|n| n.as_ref().to_owned());
      if let Some(name) = name {
        branch_lengths.insert(name, edge_ref.payload().read_arc().branch_length());
      }
    }

    assert_abs_diff_eq!(
      branch_lengths["A"].expect("A should have branch length"),
      0.1,
      epsilon = 1e-10
    );

    Ok(())
  }

  #[test]
  fn test_nwk_parse_negative_branch_length() -> Result<(), Report> {
    let input = "(A:-0.001,B:0.2)root;";
    let graph = nwk_read_str::<TestNode, TestEdge, ()>(input)?;

    assert_eq!(graph.get_nodes().len(), 3);
    assert_eq!(graph.get_edges().len(), 2);

    let output = nwk_write_str(&graph, &NwkWriteOptions::default())?;
    assert_eq!(input, output);

    Ok(())
  }

  #[test]
  fn test_nwk_parse_all_datasets() -> Result<(), Report> {
    let data_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../data");
    let trees = helpers::find_files_recursive(&data_dir, "tree.nwk")?;
    assert!(!trees.is_empty(), "data/ directory should contain tree.nwk files");

    for path in &trees {
      let graph: Graph<TestNode, TestEdge, ()> = nwk_read_file(path)?;
      assert!(!graph.get_nodes().is_empty(), "Tree {path:?} should have nodes");
      assert!(!graph.get_edges().is_empty(), "Tree {path:?} should have edges");
    }

    Ok(())
  }

  mod helpers {
    use eyre::Report;
    use std::path::{Path, PathBuf};

    pub fn find_files_recursive(dir: &Path, filename: &str) -> Result<Vec<PathBuf>, Report> {
      let mut results = Vec::new();
      for entry in std::fs::read_dir(dir)? {
        let path = entry?.path();
        if path.is_dir() {
          results.extend(find_files_recursive(&path, filename)?);
        } else if path.file_name().is_some_and(|n| n == filename) {
          results.push(path);
        }
      }
      Ok(results)
    }
  }
}
