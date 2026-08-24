#[cfg(test)]
mod tests {
  use crate::topology_order::*;
  use pretty_assertions::assert_eq;

  #[test]
  fn topology_order_descendant_count_sorts_children_ascending() -> Result<(), Report> {
    let mut graph = fixture_tree()?;
    let original = child_names(&graph, "root")?;
    TopologyOrderSpec::default().apply(&mut graph)?;
    let ordered = &graph;

    let actual = child_names(ordered, "root")?;

    assert_eq!(vec!["A", "BC", "DEF"], actual);
    assert_eq!(vec!["DEF", "A", "BC"], original);

    Ok(())
  }

  #[test]
  fn topology_order_descendant_count_reverse_sorts_children_descending() -> Result<(), Report> {
    let mut graph = fixture_tree()?;
    TopologyOrderSpec::descendant_count(true).apply(&mut graph)?;
    let ordered = &graph;

    let actual = child_names(ordered, "root")?;

    assert_eq!(vec!["DEF", "BC", "A"], actual);

    Ok(())
  }

  #[test]
  fn topology_order_keep_preserves_outbound_order() -> Result<(), Report> {
    let mut graph = fixture_tree()?;
    TopologyOrderSpec::keep().apply(&mut graph)?;
    let ordered = &graph;

    let actual = child_names(ordered, "root")?;

    assert_eq!(vec!["DEF", "A", "BC"], actual);

    Ok(())
  }

  #[test]
  fn topology_order_dag_counts_shared_descendant_once_per_child() -> Result<(), Report> {
    let mut graph = Graph::<TestNode, TestEdge, ()>::new();
    let root = graph.add_node(TestNode::new("root"));
    let left = graph.add_node(TestNode::new("left"));
    let right = graph.add_node(TestNode::new("right"));
    let shared = graph.add_node(TestNode::new("shared"));
    let right_only = graph.add_node(TestNode::new("right_only"));

    graph.add_edge(root, right, TestEdge::new())?;
    graph.add_edge(root, left, TestEdge::new())?;
    graph.add_edge(left, shared, TestEdge::new())?;
    graph.add_edge(right, shared, TestEdge::new())?;
    graph.add_edge(right, right_only, TestEdge::new())?;
    graph.build()?;

    TopologyOrderSpec::default().apply(&mut graph)?;
    let ordered = &graph;

    let actual = child_names(ordered, "root")?;

    assert_eq!(vec!["left", "right"], actual);

    Ok(())
  }

  #[test]
  fn topology_order_target_order_uses_requested_tip_order() -> Result<(), Report> {
    let mut graph = fixture_tree()?;
    let spec = TopologyOrderSpec {
      preset: TopologyOrderPreset::TargetOrder,
      target_order: vec!["D", "E", "F", "B", "C", "A"]
        .into_iter()
        .map(str::to_owned)
        .collect(),
      target_aggregate: TopologyOrderTargetAggregate::Mean,
    };
    spec.apply(&mut graph)?;
    let ordered = &graph;

    let actual = child_names(ordered, "root")?;

    assert_eq!(vec!["DEF", "BC", "A"], actual);

    Ok(())
  }

  #[test]
  fn topology_order_rejects_cycles() -> Result<(), Report> {
    let mut graph = Graph::<TestNode, TestEdge, ()>::new();
    let a = graph.add_node(TestNode::new("A"));
    let b = graph.add_node(TestNode::new("B"));
    let c = graph.add_node(TestNode::new("C"));

    graph.add_edge(a, b, TestEdge::new())?;
    graph.add_edge(b, c, TestEdge::new())?;
    graph.add_edge(c, a, TestEdge::new())?;
    graph.build()?;

    let err = TopologyOrderSpec::default().apply(&mut graph).unwrap_err();

    assert!(err.to_string().contains("directed cycle"));

    Ok(())
  }

  #[test]
  fn topology_order_height_sorts_by_subtree_depth() -> Result<(), Report> {
    let mut graph = fixture_deep_tree()?;
    let spec = TopologyOrderSpec {
      preset: TopologyOrderPreset::Height,
      ..TopologyOrderSpec::default()
    };
    spec.apply(&mut graph)?;
    let ordered = &graph;

    let actual = child_names(ordered, "root")?;

    // A is leaf (height 0), shallow has height 1, deep has height 2
    assert_eq!(vec!["A", "shallow", "deep"], actual);

    Ok(())
  }

  #[test]
  fn topology_order_height_reverse_sorts_deepest_first() -> Result<(), Report> {
    let mut graph = fixture_deep_tree()?;
    let spec = TopologyOrderSpec {
      preset: TopologyOrderPreset::HeightReverse,
      ..TopologyOrderSpec::default()
    };
    spec.apply(&mut graph)?;
    let ordered = &graph;

    let actual = child_names(ordered, "root")?;

    assert_eq!(vec!["deep", "shallow", "A"], actual);

    Ok(())
  }

  #[test]
  fn topology_order_divergence_sorts_by_total_branch_length() -> Result<(), Report> {
    // root -> [short(0.1) -> [A(0.1), B(0.1)], long(0.5) -> [C(0.2)], D(0.3)]
    // short: max divergence = 0.1 + 0.1 = 0.2
    // long:  max divergence = 0.5 + 0.2 = 0.7
    // D:     leaf, divergence = 0.0
    let mut graph = fixture_branch_length_tree()?;
    let spec = TopologyOrderSpec {
      preset: TopologyOrderPreset::Divergence,
      ..TopologyOrderSpec::default()
    };
    spec.apply(&mut graph)?;
    let ordered = &graph;

    let actual = child_names(ordered, "root")?;

    assert_eq!(vec!["D", "short", "long"], actual);

    Ok(())
  }

  #[test]
  fn topology_order_divergence_reverse_sorts_longest_first() -> Result<(), Report> {
    let mut graph = fixture_branch_length_tree()?;
    let spec = TopologyOrderSpec {
      preset: TopologyOrderPreset::DivergenceReverse,
      ..TopologyOrderSpec::default()
    };
    spec.apply(&mut graph)?;
    let ordered = &graph;

    let actual = child_names(ordered, "root")?;

    assert_eq!(vec!["long", "short", "D"], actual);

    Ok(())
  }

  #[test]
  fn topology_order_label_sorts_alphabetically() -> Result<(), Report> {
    let mut graph = fixture_tree()?;
    let spec = TopologyOrderSpec {
      preset: TopologyOrderPreset::Label,
      ..TopologyOrderSpec::default()
    };
    spec.apply(&mut graph)?;
    let ordered = &graph;

    let actual = child_names(ordered, "root")?;

    // A has label "A", BC has min label "B", DEF has min label "D"
    assert_eq!(vec!["A", "BC", "DEF"], actual);

    Ok(())
  }

  #[test]
  fn topology_order_label_reverse_sorts_descending() -> Result<(), Report> {
    let mut graph = fixture_tree()?;
    let spec = TopologyOrderSpec {
      preset: TopologyOrderPreset::LabelReverse,
      ..TopologyOrderSpec::default()
    };
    spec.apply(&mut graph)?;
    let ordered = &graph;

    let actual = child_names(ordered, "root")?;

    assert_eq!(vec!["DEF", "BC", "A"], actual);

    Ok(())
  }

  #[test]
  fn topology_order_target_order_median_uses_median_position() -> Result<(), Report> {
    let mut graph = fixture_tree()?;
    // Target order: D=0, E=1, F=2, B=3, C=4, A=5
    // DEF median of [0,1,2] = 1, BC median of [3,4] = 3.5, A = 5
    let spec = TopologyOrderSpec {
      preset: TopologyOrderPreset::TargetOrder,
      target_order: vec!["D", "E", "F", "B", "C", "A"]
        .into_iter()
        .map(str::to_owned)
        .collect(),
      target_aggregate: TopologyOrderTargetAggregate::Median,
    };
    spec.apply(&mut graph)?;
    let ordered = &graph;

    let actual = child_names(ordered, "root")?;

    assert_eq!(vec!["DEF", "BC", "A"], actual);

    Ok(())
  }

  #[test]
  fn topology_order_propagates_through_nested_levels() -> Result<(), Report> {
    let mut graph = fixture_deep_tree()?;
    TopologyOrderSpec::default().apply(&mut graph)?;
    let ordered = &graph;

    let deep_children = child_names(ordered, "deep")?;
    // mid (2 leaves) vs D (1 leaf): D first
    assert_eq!(vec!["D", "mid"], deep_children);

    let mid_children = child_names(ordered, "mid")?;
    assert_eq!(vec!["E", "F"], mid_children);

    Ok(())
  }

  #[test]
  fn topology_order_target_order_reverse_inverts_order() -> Result<(), Report> {
    let mut graph = fixture_tree()?;
    let spec = TopologyOrderSpec {
      preset: TopologyOrderPreset::TargetOrderReverse,
      target_order: vec!["D", "E", "F", "B", "C", "A"]
        .into_iter()
        .map(str::to_owned)
        .collect(),
      target_aggregate: TopologyOrderTargetAggregate::Mean,
    };
    spec.apply(&mut graph)?;
    let ordered = &graph;

    let actual = child_names(ordered, "root")?;

    assert_eq!(vec!["A", "BC", "DEF"], actual);

    Ok(())
  }

  #[test]
  fn topology_order_target_order_rejects_empty() {
    let mut graph = fixture_tree().unwrap();
    let spec = TopologyOrderSpec {
      preset: TopologyOrderPreset::TargetOrder,
      target_order: vec![],
      target_aggregate: TopologyOrderTargetAggregate::Mean,
    };
    let err = spec.apply(&mut graph).unwrap_err();
    assert!(err.to_string().contains("non-empty target order"));
  }

  #[test]
  fn topology_order_target_order_rejects_duplicate_ranking_labels() {
    let mut graph = fixture_tree().unwrap();
    let spec = TopologyOrderSpec {
      preset: TopologyOrderPreset::TargetOrder,
      target_order: vec!["A", "B", "B", "C", "D", "E", "F"]
        .into_iter()
        .map(str::to_owned)
        .collect(),
      target_aggregate: TopologyOrderTargetAggregate::Mean,
    };

    let error = spec.apply(&mut graph).unwrap_err();

    assert!(error.to_string().contains("duplicate leaf label 'B'"));
  }

  #[test]
  fn topology_order_target_order_rejects_duplicate_final_leaf_labels() -> Result<(), Report> {
    let mut graph = fixture_tree()?;
    graph
      .get_node(find_node(&graph, "C")?)
      .expect("fixture node C must exist")
      .write_arc()
      .payload()
      .write_arc()
      .0 = "B".to_owned();
    let spec = TopologyOrderSpec {
      preset: TopologyOrderPreset::TargetOrder,
      target_order: vec!["A", "B", "D", "E", "F"].into_iter().map(str::to_owned).collect(),
      target_aggregate: TopologyOrderTargetAggregate::Mean,
    };

    let error = spec.apply(&mut graph).unwrap_err();

    assert!(error.to_string().contains("final leaf label 'B' is duplicated"));
    Ok(())
  }

  #[test]
  fn topology_order_target_order_ignores_absent_ranking_labels() -> Result<(), Report> {
    let mut graph = fixture_tree()?;
    let spec = TopologyOrderSpec {
      preset: TopologyOrderPreset::TargetOrder,
      target_order: vec!["removed", "D", "E", "F", "B", "C", "A"]
        .into_iter()
        .map(str::to_owned)
        .collect(),
      target_aggregate: TopologyOrderTargetAggregate::Mean,
    };

    spec.apply(&mut graph)?;

    assert_eq!(vec!["DEF", "BC", "A"], child_names(&graph, "root")?);
    Ok(())
  }

  #[test]
  fn topology_order_is_idempotent_and_preserves_graph_data_identity() -> Result<(), Report> {
    let graph = fixture_tree()?;
    let mut graph = graph.map_data(NonCloneData);
    let data = std::ptr::from_ref(graph.data());

    TopologyOrderSpec::default().apply(&mut graph)?;
    let first = child_names_with_data(&graph, "root")?;
    TopologyOrderSpec::default().apply(&mut graph)?;
    let second = child_names_with_data(&graph, "root")?;

    assert_eq!(first, second);
    assert!(std::ptr::eq(data, graph.data()));
    Ok(())
  }

  /// root -> [deep -> [D, mid -> [E, F]], shallow -> [B, C], A]
  fn fixture_deep_tree() -> Result<Graph<TestNode, TestEdge, ()>, Report> {
    let mut graph = Graph::<TestNode, TestEdge, ()>::new();
    let root = graph.add_node(TestNode::new("root"));
    let deep = graph.add_node(TestNode::new("deep"));
    let tip_a = graph.add_node(TestNode::new("A"));
    let shallow = graph.add_node(TestNode::new("shallow"));
    let mid = graph.add_node(TestNode::new("mid"));
    let tip_d = graph.add_node(TestNode::new("D"));
    let tip_e = graph.add_node(TestNode::new("E"));
    let tip_f = graph.add_node(TestNode::new("F"));
    let tip_b = graph.add_node(TestNode::new("B"));
    let tip_c = graph.add_node(TestNode::new("C"));

    graph.add_edge(root, deep, TestEdge::new())?;
    graph.add_edge(root, tip_a, TestEdge::new())?;
    graph.add_edge(root, shallow, TestEdge::new())?;
    graph.add_edge(deep, tip_d, TestEdge::new())?;
    graph.add_edge(deep, mid, TestEdge::new())?;
    graph.add_edge(mid, tip_e, TestEdge::new())?;
    graph.add_edge(mid, tip_f, TestEdge::new())?;
    graph.add_edge(shallow, tip_b, TestEdge::new())?;
    graph.add_edge(shallow, tip_c, TestEdge::new())?;
    graph.build()?;

    Ok(graph)
  }

  /// root -> [short(0.1) -> [A(0.1), B(0.1)], long(0.5) -> [C(0.2)], D(0.3)]
  fn fixture_branch_length_tree() -> Result<Graph<TestNode, TestEdge, ()>, Report> {
    let mut graph = Graph::<TestNode, TestEdge, ()>::new();
    let root = graph.add_node(TestNode::new("root"));
    let short = graph.add_node(TestNode::new("short"));
    let long = graph.add_node(TestNode::new("long"));
    let tip_a = graph.add_node(TestNode::new("A"));
    let tip_b = graph.add_node(TestNode::new("B"));
    let tip_c = graph.add_node(TestNode::new("C"));
    let tip_d = graph.add_node(TestNode::new("D"));

    graph.add_edge(root, short, TestEdge::with_length(0.1))?;
    graph.add_edge(root, long, TestEdge::with_length(0.5))?;
    graph.add_edge(root, tip_d, TestEdge::with_length(0.3))?;
    graph.add_edge(short, tip_a, TestEdge::with_length(0.1))?;
    graph.add_edge(short, tip_b, TestEdge::with_length(0.1))?;
    graph.add_edge(long, tip_c, TestEdge::with_length(0.2))?;
    graph.build()?;

    Ok(graph)
  }

  fn fixture_tree() -> Result<Graph<TestNode, TestEdge, ()>, Report> {
    let mut graph = Graph::<TestNode, TestEdge, ()>::new();
    let root = graph.add_node(TestNode::new("root"));
    let def = graph.add_node(TestNode::new("DEF"));
    let tip_a = graph.add_node(TestNode::new("A"));
    let bc = graph.add_node(TestNode::new("BC"));
    let tip_d = graph.add_node(TestNode::new("D"));
    let tip_e = graph.add_node(TestNode::new("E"));
    let tip_f = graph.add_node(TestNode::new("F"));
    let tip_b = graph.add_node(TestNode::new("B"));
    let tip_c = graph.add_node(TestNode::new("C"));

    graph.add_edge(root, def, TestEdge::new())?;
    graph.add_edge(root, tip_a, TestEdge::new())?;
    graph.add_edge(root, bc, TestEdge::new())?;
    graph.add_edge(def, tip_d, TestEdge::new())?;
    graph.add_edge(def, tip_e, TestEdge::new())?;
    graph.add_edge(def, tip_f, TestEdge::new())?;
    graph.add_edge(bc, tip_b, TestEdge::new())?;
    graph.add_edge(bc, tip_c, TestEdge::new())?;
    graph.build()?;

    Ok(graph)
  }

  fn child_names(graph: &Graph<TestNode, TestEdge, ()>, parent_name: &str) -> Result<Vec<String>, Report> {
    let parent_key = find_node(graph, parent_name)?;
    let parent = graph
      .get_node(parent_key)
      .ok_or_else(|| make_report!("Node {parent_key} not found"))?;
    Ok(
      graph
        .children_of(&parent.read_arc())
        .into_iter()
        .map(|(node, _)| node.read_arc().payload().read_arc().0.clone())
        .collect_vec(),
    )
  }

  fn find_node(graph: &Graph<TestNode, TestEdge, ()>, name: &str) -> Result<GraphNodeKey, Report> {
    graph
      .find_node(|node| node.0 == name)
      .ok_or_else(|| make_report!("Node '{name}' not found"))
  }

  #[derive(Debug, Eq, PartialEq)]
  struct TestNode(String);

  impl TestNode {
    fn new(name: &str) -> Self {
      Self(name.to_owned())
    }
  }

  impl GraphNode for TestNode {}

  impl Named for TestNode {
    fn name(&self) -> Option<impl AsRef<str>> {
      Some(&self.0)
    }

    fn set_name(&mut self, name: Option<impl AsRef<str>>) {
      self.0 = name.map(|name| name.as_ref().to_owned()).unwrap_or_default();
    }
  }

  #[derive(Debug, PartialEq)]
  struct TestEdge {
    branch_length: Option<f64>,
  }

  impl TestEdge {
    fn new() -> Self {
      Self { branch_length: None }
    }

    fn with_length(len: f64) -> Self {
      Self {
        branch_length: Some(len),
      }
    }
  }

  impl GraphEdge for TestEdge {}

  impl HasBranchLength for TestEdge {
    fn branch_length(&self) -> Option<f64> {
      self.branch_length
    }

    fn set_branch_length(&mut self, branch_length: Option<f64>) {
      self.branch_length = branch_length;
    }
  }

  #[derive(Debug)]
  struct NonCloneData;

  fn child_names_with_data<D: Send + Sync>(
    graph: &Graph<TestNode, TestEdge, D>,
    parent_name: &str,
  ) -> Result<Vec<String>, Report> {
    let parent_key = graph
      .find_node(|node| node.0 == parent_name)
      .ok_or_else(|| make_report!("Node '{parent_name}' not found"))?;
    let parent = graph
      .get_node(parent_key)
      .ok_or_else(|| make_report!("Node {parent_key} not found"))?;
    Ok(
      graph
        .children_of(&parent.read_arc())
        .into_iter()
        .map(|(node, _)| node.read_arc().payload().read_arc().0.clone())
        .collect_vec(),
    )
  }
}
