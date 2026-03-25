#[cfg(test)]
mod tests {
  use crate::commands::timetree::output::auspice::write_auspice_json;
  use crate::commands::timetree::output::confidence::NodeConfidenceInterval;
  use crate::representation::partition::timetree::GraphTimetree;
  use crate::representation::payload::timetree::NodeTimetree;
  use approx::assert_relative_eq;
  use pretty_assertions::assert_eq;

  use treetime_graph::node::{GraphNodeKey, Named};
  use treetime_io::auspice_types::{AuspiceColoring, AuspiceTree};
  use treetime_io::json::json_read_file;

  #[test]
  fn test_auspice_metadata_structure() {
    let graph = build_simple_tree();
    let dir = tempfile::tempdir().unwrap();
    write_auspice_json(&graph, None, dir.path()).unwrap();

    let tree: AuspiceTree = json_read_file(dir.path().join("auspice_tree.json")).unwrap();
    assert_eq!(tree.data.version.as_deref(), Some("v2"));
    assert!(tree.data.meta.title.as_deref().is_some());
    assert!(tree.data.meta.updated.is_some());
    assert_eq!(tree.data.meta.panels, vec!["tree"]);
    assert_eq!(tree.data.meta.filters, vec!["bad_branch"]);

    // v0 parity: color_by=bad_branch, no distance_measure
    assert_eq!(tree.data.meta.display_defaults.color_by.as_deref(), Some("bad_branch"));
    assert!(tree.data.meta.display_defaults.distance_measure.is_none());

    // Two colorings: num_date (continuous) and bad_branch (categorical).
    // After JSON round-trip, #[serde(flatten)] other becomes empty object.
    let empty_obj = serde_json::json!({});
    let expected_colorings = vec![
      AuspiceColoring {
        title: "Date".to_owned(),
        type_: "continuous".to_owned(),
        key: "num_date".to_owned(),
        scale: vec![],
        other: empty_obj.clone(),
      },
      AuspiceColoring {
        title: "Excluded".to_owned(),
        type_: "categorical".to_owned(),
        key: "bad_branch".to_owned(),
        scale: vec![],
        other: empty_obj,
      },
    ];
    assert_eq!(tree.data.meta.colorings, expected_colorings);
  }

  #[test]
  fn test_auspice_root_node_attributes() {
    let graph = build_simple_tree();
    let dir = tempfile::tempdir().unwrap();
    write_auspice_json(&graph, None, dir.path()).unwrap();

    let tree: AuspiceTree = json_read_file(dir.path().join("auspice_tree.json")).unwrap();

    // Root has zero divergence (from NodeTimetree.div)
    assert_relative_eq!(tree.tree.node_attrs.div.unwrap(), 0.0);

    // Root has num_date from time field
    let num_date = tree.tree.node_attrs.num_date.as_ref().unwrap();
    assert_relative_eq!(num_date.value, 2000.0);
    assert!(num_date.confidence.is_none());

    // Root bad_branch = "No"
    assert_eq!(tree.tree.node_attrs.bad_branch.as_ref().unwrap().value, "No");
  }

  #[test]
  fn test_auspice_divergence_from_node_payload() {
    let graph = build_simple_tree();
    let dir = tempfile::tempdir().unwrap();
    write_auspice_json(&graph, None, dir.path()).unwrap();

    let tree: AuspiceTree = json_read_file(dir.path().join("auspice_tree.json")).unwrap();

    // Root div = 0.0 (set in make_named_node)
    assert_relative_eq!(tree.tree.node_attrs.div.unwrap(), 0.0);

    // Children read div directly from NodeTimetree.div (set in build_simple_tree)
    assert_eq!(tree.tree.children.len(), 2);
    assert_relative_eq!(tree.tree.children[0].node_attrs.div.unwrap(), 0.005);
    assert_relative_eq!(tree.tree.children[1].node_attrs.div.unwrap(), 0.010);
  }

  #[test]
  fn test_auspice_bad_branch_attribute() {
    let mut graph = GraphTimetree::new();
    let mut root = make_named_node("root", 2000.0, 0.0);
    root.bad_branch = false;
    let root_key = graph.add_node(root);

    let mut bad = make_named_node("bad_node", 2005.0, 0.003);
    bad.bad_branch = true;
    let bad_key = graph.add_node(bad);

    graph.add_edge(root_key, bad_key, make_edge(0.003)).unwrap();
    graph.build().unwrap();

    let dir = tempfile::tempdir().unwrap();
    write_auspice_json(&graph, None, dir.path()).unwrap();

    let tree: AuspiceTree = json_read_file(dir.path().join("auspice_tree.json")).unwrap();
    assert_eq!(tree.tree.node_attrs.bad_branch.as_ref().unwrap().value, "No");
    assert_eq!(
      tree.tree.children[0].node_attrs.bad_branch.as_ref().unwrap().value,
      "Yes"
    );
  }

  #[test]
  fn test_auspice_confidence_intervals_by_key() {
    let graph = build_simple_tree();

    // Use GraphNodeKey for CI lookup (keys are 0-indexed by insertion order)
    let intervals = vec![
      NodeConfidenceInterval {
        key: GraphNodeKey(0),
        name: "root".to_owned(),
        date: 2000.0,
        lower: 1998.0,
        upper: 2002.0,
      },
      NodeConfidenceInterval {
        key: GraphNodeKey(1),
        name: "child_a".to_owned(),
        date: 2005.0,
        lower: 2004.0,
        upper: 2006.0,
      },
    ];

    let dir = tempfile::tempdir().unwrap();
    write_auspice_json(&graph, Some(&intervals), dir.path()).unwrap();

    let tree: AuspiceTree = json_read_file(dir.path().join("auspice_tree.json")).unwrap();

    // Root gets CI (looked up by key, not name)
    let root_ci = tree.tree.node_attrs.num_date.as_ref().unwrap().confidence.unwrap();
    assert_relative_eq!(root_ci[0], 1998.0);
    assert_relative_eq!(root_ci[1], 2002.0);

    // child_a gets CI
    let child_a_ci = tree.tree.children[0]
      .node_attrs
      .num_date
      .as_ref()
      .unwrap()
      .confidence
      .unwrap();
    assert_relative_eq!(child_a_ci[0], 2004.0);
    assert_relative_eq!(child_a_ci[1], 2006.0);

    // child_b has no CI entry (key not in map)
    let child_b_ci = tree.tree.children[1].node_attrs.num_date.as_ref().unwrap().confidence;
    assert!(child_b_ci.is_none());

    // Structural invariant: for every node with confidence, lower <= value <= upper
    assert_ci_ordering_invariant(&tree.tree);
  }

  #[test]
  fn test_auspice_unnamed_node_gets_ci_by_key() {
    let mut graph = GraphTimetree::new();
    let root = NodeTimetree {
      time: Some(2000.0),
      ..NodeTimetree::default()
    };
    let root_key = graph.add_node(root);
    graph.build().unwrap();

    // CI for unnamed node, looked up by key
    let intervals = vec![NodeConfidenceInterval {
      key: root_key,
      name: String::new(),
      date: 2000.0,
      lower: 1999.0,
      upper: 2001.0,
    }];

    let dir = tempfile::tempdir().unwrap();
    write_auspice_json(&graph, Some(&intervals), dir.path()).unwrap();

    let tree: AuspiceTree = json_read_file(dir.path().join("auspice_tree.json")).unwrap();

    // Unnamed node gets fallback name
    assert!(tree.tree.name.starts_with("node_"));

    // CI is present despite unnamed node (key-based lookup)
    let ci = tree.tree.node_attrs.num_date.as_ref().unwrap().confidence.unwrap();
    assert_relative_eq!(ci[0], 1999.0);
    assert_relative_eq!(ci[1], 2001.0);
  }

  #[test]
  fn test_auspice_node_without_time_has_no_num_date() {
    let mut graph = GraphTimetree::new();
    let mut root = NodeTimetree::default();
    root.base.set_name(Some("root"));
    graph.add_node(root);
    graph.build().unwrap();

    let dir = tempfile::tempdir().unwrap();
    write_auspice_json(&graph, None, dir.path()).unwrap();

    let tree: AuspiceTree = json_read_file(dir.path().join("auspice_tree.json")).unwrap();
    assert!(tree.tree.node_attrs.num_date.is_none());
  }

  #[test]
  fn test_auspice_rejects_nan_div() {
    let mut graph = GraphTimetree::new();
    let node = make_named_node("bad", 2020.0, f64::NAN);
    graph.add_node(node);
    graph.build().unwrap();

    let dir = tempfile::tempdir().unwrap();
    let err = write_auspice_json(&graph, None, dir.path()).unwrap_err();
    let msg = format!("{err}");
    assert!(msg.contains("non-finite div"), "expected div error, got: {msg}");
    assert!(msg.contains("bad"), "expected node name in error, got: {msg}");
  }

  #[test]
  fn test_auspice_rejects_infinite_time() {
    let mut graph = GraphTimetree::new();
    let mut node = NodeTimetree::default();
    node.base.set_name(Some("inf_node"));
    node.time = Some(f64::INFINITY);
    graph.add_node(node);
    graph.build().unwrap();

    let dir = tempfile::tempdir().unwrap();
    let err = write_auspice_json(&graph, None, dir.path()).unwrap_err();
    let msg = format!("{err}");
    assert!(msg.contains("non-finite time"), "expected time error, got: {msg}");
    assert!(msg.contains("inf_node"), "expected node name in error, got: {msg}");
  }

  #[test]
  fn test_auspice_output_file_is_valid_json() {
    let graph = build_simple_tree();
    let dir = tempfile::tempdir().unwrap();
    write_auspice_json(&graph, None, dir.path()).unwrap();

    let path = dir.path().join("auspice_tree.json");
    assert!(path.exists());

    // Verify it parses as valid AuspiceTree
    let tree: AuspiceTree = json_read_file(&path).unwrap();
    assert!(!tree.tree.name.is_empty());
  }

  /// Build a 3-node tree: root -> child_a, root -> child_b
  /// Sets div on each node to simulate post-clock-regression state.
  fn build_simple_tree() -> GraphTimetree {
    let mut graph = GraphTimetree::new();
    let root_key = graph.add_node(make_named_node("root", 2000.0, 0.0));
    let child_a_key = graph.add_node(make_named_node("child_a", 2005.0, 0.005));
    let child_b_key = graph.add_node(make_named_node("child_b", 2010.0, 0.010));

    graph.add_edge(root_key, child_a_key, make_edge(0.005)).unwrap();
    graph.add_edge(root_key, child_b_key, make_edge(0.010)).unwrap();
    graph.build().unwrap();
    graph
  }

  /// Assert that every node with num_date.confidence satisfies lower <= value <= upper.
  fn assert_ci_ordering_invariant(node: &treetime_io::auspice_types::AuspiceTreeNode) {
    if let Some(nd) = &node.node_attrs.num_date {
      if let Some([lower, upper]) = nd.confidence {
        assert!(
          lower <= nd.value && nd.value <= upper,
          "CI ordering violated for '{}': {lower} <= {val} <= {upper}",
          node.name,
          val = nd.value,
        );
      }
    }
    for child in &node.children {
      assert_ci_ordering_invariant(child);
    }
  }

  mod helpers {
    use crate::representation::payload::timetree::{EdgeTimetree, NodeTimetree};
    use treetime_graph::node::Named;

    pub fn make_named_node(name: &str, time: f64, div: f64) -> NodeTimetree {
      let mut node = NodeTimetree::default();
      node.base.set_name(Some(name));
      node.time = Some(time);
      node.div = div;
      node
    }

    pub fn make_edge(branch_length: f64) -> EdgeTimetree {
      let mut edge = EdgeTimetree::default();
      edge.base.branch_length = Some(branch_length);
      edge
    }
  }

  use helpers::{make_edge, make_named_node};
}
