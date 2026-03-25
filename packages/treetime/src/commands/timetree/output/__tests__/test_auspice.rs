#[cfg(test)]
mod tests {
  use crate::commands::timetree::output::auspice::write_auspice_json;
  use crate::commands::timetree::output::confidence::NodeConfidenceInterval;
  use crate::representation::partition::timetree::GraphTimetree;
  use crate::representation::payload::timetree::NodeTimetree;
  use approx::assert_relative_eq;
  use pretty_assertions::assert_eq;
  use treetime_graph::node::Named;
  use treetime_io::auspice_types::AuspiceTree;
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
    assert_eq!(tree.data.meta.display_defaults.color_by.as_deref(), Some("num_date"));
    assert_eq!(
      tree.data.meta.display_defaults.distance_measure.as_deref(),
      Some("num_date")
    );

    // Two colorings: num_date (continuous) and bad_branch (categorical)
    assert_eq!(tree.data.meta.colorings.len(), 2);
    assert_eq!(tree.data.meta.colorings[0].key, "num_date");
    assert_eq!(tree.data.meta.colorings[0].type_, "continuous");
    assert_eq!(tree.data.meta.colorings[1].key, "bad_branch");
    assert_eq!(tree.data.meta.colorings[1].type_, "categorical");
  }

  #[test]
  fn test_auspice_root_node_attributes() {
    let graph = build_simple_tree();
    let dir = tempfile::tempdir().unwrap();
    write_auspice_json(&graph, None, dir.path()).unwrap();

    let tree: AuspiceTree = json_read_file(dir.path().join("auspice_tree.json")).unwrap();

    // Root has zero divergence
    assert_relative_eq!(tree.tree.node_attrs.div.unwrap(), 0.0);

    // Root has num_date from time field
    let num_date = tree.tree.node_attrs.num_date.as_ref().unwrap();
    assert_relative_eq!(num_date.value, 2000.0);
    assert!(num_date.confidence.is_none());

    // Root bad_branch = "No"
    assert_eq!(tree.tree.node_attrs.bad_branch.as_ref().unwrap().value, "No");
  }

  #[test]
  fn test_auspice_cumulative_divergence() {
    let graph = build_simple_tree();
    let dir = tempfile::tempdir().unwrap();
    write_auspice_json(&graph, None, dir.path()).unwrap();

    let tree: AuspiceTree = json_read_file(dir.path().join("auspice_tree.json")).unwrap();

    // Root div = 0.0
    assert_relative_eq!(tree.tree.node_attrs.div.unwrap(), 0.0);

    // Children: root_div + branch_length
    // Child A: 0.0 + 0.005 = 0.005
    // Child B: 0.0 + 0.010 = 0.010
    assert_eq!(tree.tree.children.len(), 2);
    let child_a = &tree.tree.children[0];
    let child_b = &tree.tree.children[1];

    // Children are in insertion order (pre-order traversal)
    assert_relative_eq!(child_a.node_attrs.div.unwrap(), 0.005);
    assert_relative_eq!(child_b.node_attrs.div.unwrap(), 0.010);
  }

  #[test]
  fn test_auspice_bad_branch_attribute() {
    let mut graph = GraphTimetree::new();
    let mut root = make_named_node("root", 2000.0);
    root.bad_branch = false;
    let root_key = graph.add_node(root);

    let mut bad = make_named_node("bad_node", 2005.0);
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
  fn test_auspice_confidence_intervals() {
    let graph = build_simple_tree();
    let intervals = vec![
      NodeConfidenceInterval {
        name: "root".to_owned(),
        date: 2000.0,
        lower: 1998.0,
        upper: 2002.0,
      },
      NodeConfidenceInterval {
        name: "child_a".to_owned(),
        date: 2005.0,
        lower: 2004.0,
        upper: 2006.0,
      },
    ];

    let dir = tempfile::tempdir().unwrap();
    write_auspice_json(&graph, Some(&intervals), dir.path()).unwrap();

    let tree: AuspiceTree = json_read_file(dir.path().join("auspice_tree.json")).unwrap();

    // Root gets CI
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

    // child_b has no CI entry
    let child_b_ci = tree.tree.children[1].node_attrs.num_date.as_ref().unwrap().confidence;
    assert!(child_b_ci.is_none());
  }

  #[test]
  fn test_auspice_unnamed_node_gets_fallback_name() {
    let mut graph = GraphTimetree::new();
    graph.add_node(NodeTimetree {
      time: Some(2000.0),
      ..NodeTimetree::default()
    });
    graph.build().unwrap();

    let dir = tempfile::tempdir().unwrap();
    write_auspice_json(&graph, None, dir.path()).unwrap();

    let tree: AuspiceTree = json_read_file(dir.path().join("auspice_tree.json")).unwrap();
    // Fallback name is "node_<key>"
    assert!(tree.tree.name.starts_with("node_"));
  }

  #[test]
  fn test_auspice_node_without_time_has_no_num_date() {
    let mut graph = GraphTimetree::new();
    let mut root = NodeTimetree::default();
    root.base.set_name(Some("root"));
    root.time = None;
    graph.add_node(root);
    graph.build().unwrap();

    let dir = tempfile::tempdir().unwrap();
    write_auspice_json(&graph, None, dir.path()).unwrap();

    let tree: AuspiceTree = json_read_file(dir.path().join("auspice_tree.json")).unwrap();
    assert!(tree.tree.node_attrs.num_date.is_none());
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
  fn build_simple_tree() -> GraphTimetree {
    let mut graph = GraphTimetree::new();
    let root_key = graph.add_node(make_named_node("root", 2000.0));
    let child_a_key = graph.add_node(make_named_node("child_a", 2005.0));
    let child_b_key = graph.add_node(make_named_node("child_b", 2010.0));

    graph.add_edge(root_key, child_a_key, make_edge(0.005)).unwrap();
    graph.add_edge(root_key, child_b_key, make_edge(0.010)).unwrap();
    graph.build().unwrap();
    graph
  }

  mod helpers {
    use crate::representation::payload::timetree::{EdgeTimetree, NodeTimetree};
    use treetime_graph::node::Named;

    pub fn make_named_node(name: &str, time: f64) -> NodeTimetree {
      let mut node = NodeTimetree::default();
      node.base.set_name(Some(name));
      node.time = Some(time);
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
