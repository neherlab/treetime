#[cfg(test)]
mod tests {
  use crate::commands::shared::output::{CommandKind, OutputCoreArgs, OutputSelection};
  use pretty_assertions::assert_eq;
  use std::path::PathBuf;
  use treetime_io::graph::TreeWriteKind;
  use treetime_io::nwk::NwkStyle;
  use treetime_utils::assert_error;

  fn make_test_graph() -> treetime_graph::graph::Graph<TestNode, TestEdge, ()> {
    use treetime_graph::graph::Graph;
    let mut graph: Graph<TestNode, TestEdge, ()> = Graph::new();
    let root = graph.add_node(TestNode {
      name: Some("root".to_owned()),
    });
    let leaf_a = graph.add_node(TestNode {
      name: Some("A".to_owned()),
    });
    let leaf_b = graph.add_node(TestNode {
      name: Some("B".to_owned()),
    });
    let _ = graph.add_edge(
      root,
      leaf_a,
      TestEdge {
        branch_length: Some(0.1),
      },
    );
    let _ = graph.add_edge(
      root,
      leaf_b,
      TestEdge {
        branch_length: Some(0.2),
      },
    );
    graph
  }

  #[test]
  fn test_resolve_output_all_default_selection() {
    let dir = tempfile::TempDir::new().unwrap();
    let graph = make_test_graph();
    let args = OutputCoreArgs {
      output_all: Some(dir.path().to_path_buf()),
      ..Default::default()
    };
    let resolved = args.resolve(CommandKind::Ancestral, &graph, &[], None).unwrap();

    assert!(
      resolved.tree_outputs.contains_key(&TreeWriteKind::nwk(NwkStyle::Plain)),
      "default selection should include Nwk"
    );
    assert!(
      resolved
        .tree_outputs
        .contains_key(&TreeWriteKind::nexus(NwkStyle::Plain)),
      "default selection should include Nexus"
    );
    assert!(
      !resolved.tree_outputs.contains_key(&TreeWriteKind::GraphJson),
      "GraphJson should not be in defaults"
    );
    assert!(
      !resolved.tree_outputs.contains_key(&TreeWriteKind::Dot),
      "Dot should not be in defaults"
    );
  }

  #[test]
  fn test_resolve_output_all_with_selection_filter() {
    let dir = tempfile::TempDir::new().unwrap();
    let graph = make_test_graph();
    let args = OutputCoreArgs {
      output_all: Some(dir.path().to_path_buf()),
      output_selection: vec![OutputSelection::Nwk],
      ..Default::default()
    };
    let resolved = args.resolve(CommandKind::Ancestral, &graph, &[], None).unwrap();

    assert_eq!(resolved.tree_outputs.len(), 1);
    assert!(resolved.tree_outputs.contains_key(&TreeWriteKind::nwk(NwkStyle::Plain)));
  }

  #[test]
  fn test_resolve_per_file_override_wins() {
    let dir = tempfile::TempDir::new().unwrap();
    let custom_path = dir.path().join("custom.nwk");
    let graph = make_test_graph();
    let args = OutputCoreArgs {
      output_all: Some(dir.path().to_path_buf()),
      output_tree_nwk: Some(custom_path.clone()),
      ..Default::default()
    };
    let resolved = args.resolve(CommandKind::Ancestral, &graph, &[], None).unwrap();

    let nwk_path = &resolved.tree_outputs[&TreeWriteKind::nwk(NwkStyle::Plain)];
    assert_eq!(nwk_path, &custom_path, "per-file override should take precedence");
  }

  #[test]
  fn test_resolve_per_file_outside_selection_honored() {
    let dir = tempfile::TempDir::new().unwrap();
    let nexus_path = dir.path().join("extra.nexus");
    let graph = make_test_graph();
    let args = OutputCoreArgs {
      output_all: Some(dir.path().to_path_buf()),
      output_selection: vec![OutputSelection::Nwk],
      output_tree_nexus: Some(nexus_path.clone()),
      ..Default::default()
    };
    let resolved = args.resolve(CommandKind::Ancestral, &graph, &[], None).unwrap();

    assert!(
      resolved
        .tree_outputs
        .contains_key(&TreeWriteKind::nexus(NwkStyle::Plain)),
      "explicit per-file nexus should be honored even outside selection"
    );
    assert_eq!(
      &resolved.tree_outputs[&TreeWriteKind::nexus(NwkStyle::Plain)],
      &nexus_path
    );
  }

  #[test]
  fn test_resolve_per_file_only_no_output_all() {
    let graph = make_test_graph();
    let args = OutputCoreArgs {
      output_tree_nwk: Some(PathBuf::from("/tmp/test.nwk")),
      ..Default::default()
    };
    let resolved = args.resolve(CommandKind::Ancestral, &graph, &[], None).unwrap();

    assert_eq!(resolved.tree_outputs.len(), 1);
    assert!(resolved.tree_outputs.contains_key(&TreeWriteKind::nwk(NwkStyle::Plain)));
  }

  #[test]
  fn test_resolve_no_outputs_errors() {
    let graph = make_test_graph();
    let args = OutputCoreArgs::default();
    let result = args.resolve(CommandKind::Ancestral, &graph, &[], None);
    assert_error!(
      result,
      "No output flags provided. At least one is required: --output-all or one of the --output-tree-* / --output-* flags"
    );
  }

  #[test]
  fn test_resolve_non_tree_outputs_gated() {
    let dir = tempfile::TempDir::new().unwrap();
    let graph = make_test_graph();
    let gtr_path = dir.path().join("my.gtr.json");
    let args = OutputCoreArgs {
      output_all: Some(dir.path().to_path_buf()),
      ..Default::default()
    };
    let resolved = args
      .resolve(
        CommandKind::Ancestral,
        &graph,
        &[(OutputSelection::Gtr, Some(gtr_path.as_path()))],
        None,
      )
      .unwrap();

    assert_eq!(
      &resolved.non_tree_outputs[&OutputSelection::Gtr],
      &gtr_path,
      "explicit non-tree path should appear in resolved outputs"
    );
  }

  #[test]
  fn test_resolve_non_tree_default_path_from_output_all() {
    let dir = tempfile::TempDir::new().unwrap();
    let graph = make_test_graph();
    let args = OutputCoreArgs {
      output_all: Some(dir.path().to_path_buf()),
      ..Default::default()
    };
    let resolved = args
      .resolve(CommandKind::Ancestral, &graph, &[(OutputSelection::Gtr, None)], None)
      .unwrap();

    let expected = dir.path().join("annotated_tree.gtr.json");
    assert_eq!(
      &resolved.non_tree_outputs[&OutputSelection::Gtr],
      &expected,
      "non-tree default path should use stem + extension"
    );
  }

  #[test]
  fn test_resolve_unavailable_non_tree_errors() {
    let graph = make_test_graph();
    let args = OutputCoreArgs {
      output_tree_nwk: Some(PathBuf::from("/tmp/out.nwk")),
      ..Default::default()
    };
    let result = args.resolve(
      CommandKind::Clock,
      &graph,
      &[(
        OutputSelection::AugurNodeData,
        Some(std::path::Path::new("/tmp/node.json")),
      )],
      None,
    );
    assert_error!(
      result,
      "Output format '--output-augur-node-data' is not available for this command. Available non-tree outputs: --output-clock-model, --output-clock-csv"
    );
  }

  #[test]
  fn test_resolve_output_selection_all_expands_to_defaults() {
    let dir = tempfile::TempDir::new().unwrap();
    let graph = make_test_graph();
    let args = OutputCoreArgs {
      output_all: Some(dir.path().to_path_buf()),
      output_selection: vec![OutputSelection::All],
      ..Default::default()
    };
    let resolved = args.resolve(CommandKind::Ancestral, &graph, &[], None).unwrap();

    assert!(
      resolved.tree_outputs.contains_key(&TreeWriteKind::nwk(NwkStyle::Plain)),
      "All should expand to include Nwk"
    );
    assert!(
      resolved
        .tree_outputs
        .contains_key(&TreeWriteKind::nexus(NwkStyle::Plain)),
      "All should expand to include Nexus"
    );
  }

  #[test]
  fn test_resolve_timetree_defaults_include_auspice() {
    let dir = tempfile::TempDir::new().unwrap();
    let graph = make_test_graph();
    let args = OutputCoreArgs {
      output_all: Some(dir.path().to_path_buf()),
      ..Default::default()
    };
    let resolved = args.resolve(CommandKind::Timetree, &graph, &[], None).unwrap();

    assert!(
      resolved.tree_outputs.contains_key(&TreeWriteKind::Auspice),
      "timetree defaults should include Auspice"
    );
  }

  #[test]
  fn test_resolve_default_paths_use_command_stem() {
    let dir = tempfile::TempDir::new().unwrap();
    let graph = make_test_graph();
    let args = OutputCoreArgs {
      output_all: Some(dir.path().to_path_buf()),
      output_selection: vec![OutputSelection::Nwk],
      ..Default::default()
    };

    let resolved = args.resolve(CommandKind::Timetree, &graph, &[], None).unwrap();
    let nwk_path = &resolved.tree_outputs[&TreeWriteKind::nwk(NwkStyle::Plain)];
    assert_eq!(nwk_path, &dir.path().join("timetree.nwk"));

    let args2 = OutputCoreArgs {
      output_all: Some(dir.path().to_path_buf()),
      output_selection: vec![OutputSelection::Nwk],
      ..Default::default()
    };
    let resolved2 = args2.resolve(CommandKind::Clock, &graph, &[], None).unwrap();
    let nwk_path2 = &resolved2.tree_outputs[&TreeWriteKind::nwk(NwkStyle::Plain)];
    assert_eq!(nwk_path2, &dir.path().join("rerooted.nwk"));
  }

  use treetime_graph::edge::{GraphEdge, HasBranchLength};
  use treetime_graph::node::{GraphNode, Named};
  use treetime_io::nwk::{EdgeFromNwk, NodeFromNwk};

  #[derive(Clone, Debug, Default, serde::Serialize)]
  struct TestNode {
    name: Option<String>,
  }

  impl GraphNode for TestNode {}

  impl Named for TestNode {
    fn name(&self) -> Option<impl AsRef<str>> {
      self.name.as_deref()
    }

    fn set_name(&mut self, name: Option<impl AsRef<str>>) {
      self.name = name.map(|n| n.as_ref().to_owned());
    }
  }

  impl NodeFromNwk for TestNode {
    fn from_nwk(
      name: Option<impl AsRef<str>>,
      _: &std::collections::BTreeMap<String, String>,
    ) -> Result<Self, eyre::Report> {
      Ok(Self {
        name: name.map(|n| n.as_ref().to_owned()),
      })
    }
  }

  #[derive(Clone, Debug, Default, serde::Serialize)]
  struct TestEdge {
    branch_length: Option<f64>,
  }

  impl GraphEdge for TestEdge {}

  impl HasBranchLength for TestEdge {
    fn branch_length(&self) -> Option<f64> {
      self.branch_length
    }

    fn set_branch_length(&mut self, bl: Option<f64>) {
      self.branch_length = bl;
    }
  }

  impl EdgeFromNwk for TestEdge {
    fn from_nwk(weight: Option<f64>) -> Result<Self, eyre::Report> {
      Ok(Self { branch_length: weight })
    }
  }
}
