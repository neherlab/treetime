#[cfg(test)]
mod tests {
  use crate::topology_order::*;
  use pretty_assertions::assert_eq;

  #[test]
  fn topology_order_descendant_count_sorts_children_ascending() -> Result<(), Report> {
    let (mut graph, names) = fixture_tree()?;
    let original = child_names(&graph, &names, "root")?;
    let __bl = edge_branch_lengths(&graph);
    TopologyOrderSpec::default().apply(&mut graph, &names, &__bl)?;
    let ordered = &graph;

    let actual = child_names(ordered, &names, "root")?;

    assert_eq!(vec!["A", "BC", "DEF"], actual);
    assert_eq!(vec!["DEF", "A", "BC"], original);

    Ok(())
  }

  #[test]
  fn topology_order_descendant_count_reverse_sorts_children_descending() -> Result<(), Report> {
    let (mut graph, names) = fixture_tree()?;
    let __bl = edge_branch_lengths(&graph);
    TopologyOrderSpec::descendant_count(true).apply(&mut graph, &names, &__bl)?;
    let ordered = &graph;

    let actual = child_names(ordered, &names, "root")?;

    assert_eq!(vec!["DEF", "BC", "A"], actual);

    Ok(())
  }

  #[test]
  fn topology_order_keep_preserves_outbound_order() -> Result<(), Report> {
    let (mut graph, names) = fixture_tree()?;
    let __bl = edge_branch_lengths(&graph);
    TopologyOrderSpec::keep().apply(&mut graph, &names, &__bl)?;
    let ordered = &graph;

    let actual = child_names(ordered, &names, "root")?;

    assert_eq!(vec!["DEF", "A", "BC"], actual);

    Ok(())
  }

  #[test]
  fn topology_order_dag_counts_shared_descendant_once_per_child() -> Result<(), Report> {
    let mut graph = Graph::<()>::new();
    let root = graph.add_node();
    let left = graph.add_node();
    let right = graph.add_node();
    let shared = graph.add_node();
    let right_only = graph.add_node();

    graph.add_edge(root, right)?;
    graph.add_edge(root, left)?;
    graph.add_edge(left, shared)?;
    graph.add_edge(right, shared)?;
    graph.add_edge(right, right_only)?;
    graph.build()?;

    let names = make_names(vec![
      (root, "root"),
      (left, "left"),
      (right, "right"),
      (shared, "shared"),
      (right_only, "right_only"),
    ]);
    let __bl = edge_branch_lengths(&graph);
    TopologyOrderSpec::default().apply(&mut graph, &names, &__bl)?;
    let ordered = &graph;

    let actual = child_names(ordered, &names, "root")?;

    assert_eq!(vec!["left", "right"], actual);

    Ok(())
  }

  #[test]
  fn topology_order_target_order_uses_requested_tip_order() -> Result<(), Report> {
    let (mut graph, names) = fixture_tree()?;
    let spec = TopologyOrderSpec {
      preset: TopologyOrderPreset::TargetOrder,
      target_order: vec!["D", "E", "F", "B", "C", "A"]
        .into_iter()
        .map(str::to_owned)
        .collect(),
      target_aggregate: TopologyOrderTargetAggregate::Mean,
    };
    let __bl = edge_branch_lengths(&graph);
    spec.apply(&mut graph, &names, &__bl)?;
    let ordered = &graph;

    let actual = child_names(ordered, &names, "root")?;

    assert_eq!(vec!["DEF", "BC", "A"], actual);

    Ok(())
  }

  #[test]
  fn topology_order_rejects_cycles() -> Result<(), Report> {
    let mut graph = Graph::<()>::new();
    let a = graph.add_node();
    let b = graph.add_node();
    let c = graph.add_node();

    graph.add_edge(a, b)?;
    graph.add_edge(b, c)?;
    graph.add_edge(c, a)?;
    graph.build()?;

    let names = make_names(vec![(a, "A"), (b, "B"), (c, "C")]);
    let __bl = edge_branch_lengths(&graph);
    let err = TopologyOrderSpec::default()
      .apply(&mut graph, &names, &__bl)
      .unwrap_err();

    assert!(err.to_string().contains("directed cycle"));

    Ok(())
  }

  #[test]
  fn topology_order_height_sorts_by_subtree_depth() -> Result<(), Report> {
    let (mut graph, names) = fixture_deep_tree()?;
    let spec = TopologyOrderSpec {
      preset: TopologyOrderPreset::Height,
      ..TopologyOrderSpec::default()
    };
    let __bl = edge_branch_lengths(&graph);
    spec.apply(&mut graph, &names, &__bl)?;
    let ordered = &graph;

    let actual = child_names(ordered, &names, "root")?;

    // A is leaf (height 0), shallow has height 1, deep has height 2
    assert_eq!(vec!["A", "shallow", "deep"], actual);

    Ok(())
  }

  #[test]
  fn topology_order_height_reverse_sorts_deepest_first() -> Result<(), Report> {
    let (mut graph, names) = fixture_deep_tree()?;
    let spec = TopologyOrderSpec {
      preset: TopologyOrderPreset::HeightReverse,
      ..TopologyOrderSpec::default()
    };
    let __bl = edge_branch_lengths(&graph);
    spec.apply(&mut graph, &names, &__bl)?;
    let ordered = &graph;

    let actual = child_names(ordered, &names, "root")?;

    assert_eq!(vec!["deep", "shallow", "A"], actual);

    Ok(())
  }

  #[test]
  fn topology_order_divergence_sorts_by_total_branch_length() -> Result<(), Report> {
    // root -> [short(0.1) -> [A(0.1), B(0.1)], long(0.5) -> [C(0.2)], D(0.3)]
    // short: max divergence = 0.1 + 0.1 = 0.2
    // long:  max divergence = 0.5 + 0.2 = 0.7
    // D:     leaf, divergence = 0.0
    let (mut graph, names, __bl) = fixture_branch_length_tree()?;
    let spec = TopologyOrderSpec {
      preset: TopologyOrderPreset::Divergence,
      ..TopologyOrderSpec::default()
    };

    spec.apply(&mut graph, &names, &__bl)?;
    let ordered = &graph;

    let actual = child_names(ordered, &names, "root")?;

    assert_eq!(vec!["D", "short", "long"], actual);

    Ok(())
  }

  #[test]
  fn topology_order_divergence_reverse_sorts_longest_first() -> Result<(), Report> {
    let (mut graph, names, __bl) = fixture_branch_length_tree()?;
    let spec = TopologyOrderSpec {
      preset: TopologyOrderPreset::DivergenceReverse,
      ..TopologyOrderSpec::default()
    };

    spec.apply(&mut graph, &names, &__bl)?;
    let ordered = &graph;

    let actual = child_names(ordered, &names, "root")?;

    assert_eq!(vec!["long", "short", "D"], actual);

    Ok(())
  }

  #[test]
  fn topology_order_label_sorts_alphabetically() -> Result<(), Report> {
    let (mut graph, names) = fixture_tree()?;
    let spec = TopologyOrderSpec {
      preset: TopologyOrderPreset::Label,
      ..TopologyOrderSpec::default()
    };
    let __bl = edge_branch_lengths(&graph);
    spec.apply(&mut graph, &names, &__bl)?;
    let ordered = &graph;

    let actual = child_names(ordered, &names, "root")?;

    // A has label "A", BC has min label "B", DEF has min label "D"
    assert_eq!(vec!["A", "BC", "DEF"], actual);

    Ok(())
  }

  #[test]
  fn topology_order_label_reverse_sorts_descending() -> Result<(), Report> {
    let (mut graph, names) = fixture_tree()?;
    let spec = TopologyOrderSpec {
      preset: TopologyOrderPreset::LabelReverse,
      ..TopologyOrderSpec::default()
    };
    let __bl = edge_branch_lengths(&graph);
    spec.apply(&mut graph, &names, &__bl)?;
    let ordered = &graph;

    let actual = child_names(ordered, &names, "root")?;

    assert_eq!(vec!["DEF", "BC", "A"], actual);

    Ok(())
  }

  #[test]
  fn topology_order_target_order_median_uses_median_position() -> Result<(), Report> {
    let (mut graph, names) = fixture_tree()?;
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
    let __bl = edge_branch_lengths(&graph);
    spec.apply(&mut graph, &names, &__bl)?;
    let ordered = &graph;

    let actual = child_names(ordered, &names, "root")?;

    assert_eq!(vec!["DEF", "BC", "A"], actual);

    Ok(())
  }

  #[test]
  fn topology_order_propagates_through_nested_levels() -> Result<(), Report> {
    let (mut graph, names) = fixture_deep_tree()?;
    let __bl = edge_branch_lengths(&graph);
    TopologyOrderSpec::default().apply(&mut graph, &names, &__bl)?;
    let ordered = &graph;

    let deep_children = child_names(ordered, &names, "deep")?;
    // mid (2 leaves) vs D (1 leaf): D first
    assert_eq!(vec!["D", "mid"], deep_children);

    let mid_children = child_names(ordered, &names, "mid")?;
    assert_eq!(vec!["E", "F"], mid_children);

    Ok(())
  }

  #[test]
  fn topology_order_target_order_reverse_inverts_order() -> Result<(), Report> {
    let (mut graph, names) = fixture_tree()?;
    let spec = TopologyOrderSpec {
      preset: TopologyOrderPreset::TargetOrderReverse,
      target_order: vec!["D", "E", "F", "B", "C", "A"]
        .into_iter()
        .map(str::to_owned)
        .collect(),
      target_aggregate: TopologyOrderTargetAggregate::Mean,
    };
    let __bl = edge_branch_lengths(&graph);
    spec.apply(&mut graph, &names, &__bl)?;
    let ordered = &graph;

    let actual = child_names(ordered, &names, "root")?;

    assert_eq!(vec!["A", "BC", "DEF"], actual);

    Ok(())
  }

  #[test]
  fn topology_order_target_order_rejects_empty() {
    let (mut graph, names) = fixture_tree().unwrap();
    let spec = TopologyOrderSpec {
      preset: TopologyOrderPreset::TargetOrder,
      target_order: vec![],
      target_aggregate: TopologyOrderTargetAggregate::Mean,
    };
    let __bl = edge_branch_lengths(&graph);
    let err = spec.apply(&mut graph, &names, &__bl).unwrap_err();
    assert!(err.to_string().contains("non-empty target order"));
  }

  #[test]
  fn topology_order_target_order_rejects_duplicate_ranking_labels() {
    let (mut graph, names) = fixture_tree().unwrap();
    let spec = TopologyOrderSpec {
      preset: TopologyOrderPreset::TargetOrder,
      target_order: vec!["A", "B", "B", "C", "D", "E", "F"]
        .into_iter()
        .map(str::to_owned)
        .collect(),
      target_aggregate: TopologyOrderTargetAggregate::Mean,
    };
    let __bl = edge_branch_lengths(&graph);
    let error = spec.apply(&mut graph, &names, &__bl).unwrap_err();

    assert!(error.to_string().contains("duplicate leaf label 'B'"));
  }

  #[test]
  fn topology_order_target_order_rejects_duplicate_final_leaf_labels() -> Result<(), Report> {
    let (mut graph, mut names) = fixture_tree()?;
    let node_c = find_node(&names, "C")?;
    names.insert(node_c, Some("B".to_owned()));
    let spec = TopologyOrderSpec {
      preset: TopologyOrderPreset::TargetOrder,
      target_order: vec!["A", "B", "D", "E", "F"].into_iter().map(str::to_owned).collect(),
      target_aggregate: TopologyOrderTargetAggregate::Mean,
    };

    let __bl = edge_branch_lengths(&graph);
    let error = spec.apply(&mut graph, &names, &__bl).unwrap_err();

    assert!(error.to_string().contains("final leaf label 'B' is duplicated"));
    Ok(())
  }

  #[test]
  fn topology_order_target_order_ignores_absent_ranking_labels() -> Result<(), Report> {
    let (mut graph, names) = fixture_tree()?;
    let spec = TopologyOrderSpec {
      preset: TopologyOrderPreset::TargetOrder,
      target_order: vec!["removed", "D", "E", "F", "B", "C", "A"]
        .into_iter()
        .map(str::to_owned)
        .collect(),
      target_aggregate: TopologyOrderTargetAggregate::Mean,
    };
    let __bl = edge_branch_lengths(&graph);
    spec.apply(&mut graph, &names, &__bl)?;

    assert_eq!(vec!["DEF", "BC", "A"], child_names(&graph, &names, "root")?);
    Ok(())
  }

  #[test]
  fn topology_order_is_idempotent_and_preserves_graph_data_identity() -> Result<(), Report> {
    let (graph, names) = fixture_tree()?;
    let mut graph = graph.map_data(NonCloneData);
    let data = std::ptr::from_ref(graph.data());

    let __bl = edge_branch_lengths(&graph);
    TopologyOrderSpec::default().apply(&mut graph, &names, &__bl)?;
    let first = child_names_with_data(&graph, &names, "root")?;
    let __bl = edge_branch_lengths(&graph);
    TopologyOrderSpec::default().apply(&mut graph, &names, &__bl)?;
    let second = child_names_with_data(&graph, &names, "root")?;

    assert_eq!(first, second);
    assert!(std::ptr::eq(data, graph.data()));
    Ok(())
  }

  /// root -> [deep -> [D, mid -> [E, F]], shallow -> [B, C], A]
  fn fixture_deep_tree() -> Result<(Graph<()>, BTreeMap<GraphNodeKey, Option<String>>), Report> {
    let mut graph = Graph::<()>::new();
    let root = graph.add_node();
    let deep = graph.add_node();
    let tip_a = graph.add_node();
    let shallow = graph.add_node();
    let mid = graph.add_node();
    let tip_d = graph.add_node();
    let tip_e = graph.add_node();
    let tip_f = graph.add_node();
    let tip_b = graph.add_node();
    let tip_c = graph.add_node();

    graph.add_edge(root, deep)?;
    graph.add_edge(root, tip_a)?;
    graph.add_edge(root, shallow)?;
    graph.add_edge(deep, tip_d)?;
    graph.add_edge(deep, mid)?;
    graph.add_edge(mid, tip_e)?;
    graph.add_edge(mid, tip_f)?;
    graph.add_edge(shallow, tip_b)?;
    graph.add_edge(shallow, tip_c)?;
    graph.build()?;

    let names = make_names(vec![
      (root, "root"),
      (deep, "deep"),
      (tip_a, "A"),
      (shallow, "shallow"),
      (mid, "mid"),
      (tip_d, "D"),
      (tip_e, "E"),
      (tip_f, "F"),
      (tip_b, "B"),
      (tip_c, "C"),
    ]);
    Ok((graph, names))
  }

  /// root -> [short(0.1) -> [A(0.1), B(0.1)], long(0.5) -> [C(0.2)], D(0.3)]
  #[allow(clippy::type_complexity)]
  fn fixture_branch_length_tree() -> Result<
    (
      Graph<()>,
      BTreeMap<GraphNodeKey, Option<String>>,
      BTreeMap<GraphEdgeKey, Option<f64>>,
    ),
    Report,
  > {
    let mut graph = Graph::<()>::new();
    let root = graph.add_node();
    let short = graph.add_node();
    let long = graph.add_node();
    let tip_a = graph.add_node();
    let tip_b = graph.add_node();
    let tip_c = graph.add_node();
    let tip_d = graph.add_node();

    let branch_lengths = make_branch_lengths(vec![
      (graph.add_edge(root, short)?, 0.1),
      (graph.add_edge(root, long)?, 0.5),
      (graph.add_edge(root, tip_d)?, 0.3),
      (graph.add_edge(short, tip_a)?, 0.1),
      (graph.add_edge(short, tip_b)?, 0.1),
      (graph.add_edge(long, tip_c)?, 0.2),
    ]);
    graph.build()?;

    let names = make_names(vec![
      (root, "root"),
      (short, "short"),
      (long, "long"),
      (tip_a, "A"),
      (tip_b, "B"),
      (tip_c, "C"),
      (tip_d, "D"),
    ]);
    Ok((graph, names, branch_lengths))
  }

  fn fixture_tree() -> Result<(Graph<()>, BTreeMap<GraphNodeKey, Option<String>>), Report> {
    let mut graph = Graph::<()>::new();
    let root = graph.add_node();
    let def = graph.add_node();
    let tip_a = graph.add_node();
    let bc = graph.add_node();
    let tip_d = graph.add_node();
    let tip_e = graph.add_node();
    let tip_f = graph.add_node();
    let tip_b = graph.add_node();
    let tip_c = graph.add_node();

    graph.add_edge(root, def)?;
    graph.add_edge(root, tip_a)?;
    graph.add_edge(root, bc)?;
    graph.add_edge(def, tip_d)?;
    graph.add_edge(def, tip_e)?;
    graph.add_edge(def, tip_f)?;
    graph.add_edge(bc, tip_b)?;
    graph.add_edge(bc, tip_c)?;
    graph.build()?;

    let names = make_names(vec![
      (root, "root"),
      (def, "DEF"),
      (tip_a, "A"),
      (bc, "BC"),
      (tip_d, "D"),
      (tip_e, "E"),
      (tip_f, "F"),
      (tip_b, "B"),
      (tip_c, "C"),
    ]);
    Ok((graph, names))
  }

  fn child_names(
    graph: &Graph<()>,
    names: &BTreeMap<GraphNodeKey, Option<String>>,
    parent_name: &str,
  ) -> Result<Vec<String>, Report> {
    child_names_with_data(graph, names, parent_name)
  }

  fn child_names_with_data<D: Send + Sync>(
    graph: &Graph<D>,
    names: &BTreeMap<GraphNodeKey, Option<String>>,
    parent_name: &str,
  ) -> Result<Vec<String>, Report> {
    let parent_key = find_node(names, parent_name)?;
    let parent = graph
      .get_node(parent_key)
      .ok_or_else(|| make_report!("Node {parent_key} not found"))?;
    Ok(
      graph
        .children_of(&parent.read_arc())
        .into_iter()
        .map(|(node, _)| node_name(names, node.read_arc().key()))
        .collect_vec(),
    )
  }

  fn find_node(names: &BTreeMap<GraphNodeKey, Option<String>>, name: &str) -> Result<GraphNodeKey, Report> {
    names
      .iter()
      .find_map(|(key, value)| (value.as_deref() == Some(name)).then_some(*key))
      .ok_or_else(|| make_report!("Node '{name}' not found"))
  }

  fn node_name(names: &BTreeMap<GraphNodeKey, Option<String>>, key: GraphNodeKey) -> String {
    names[&key].clone().unwrap_or_default()
  }

  /// Build the node-name value map the production parse threads to `TopologyOrderSpec::apply`.
  fn make_names(pairs: Vec<(GraphNodeKey, &str)>) -> BTreeMap<GraphNodeKey, Option<String>> {
    pairs
      .into_iter()
      .map(|(key, name)| (key, Some(name.to_owned())))
      .collect()
  }

  /// Build the edge branch-length value map the production parse threads to `TopologyOrderSpec::apply`.
  fn make_branch_lengths(pairs: Vec<(GraphEdgeKey, f64)>) -> BTreeMap<GraphEdgeKey, Option<f64>> {
    pairs.into_iter().map(|(key, len)| (key, Some(len))).collect()
  }

  /// The edge branch-length value map for a graph whose edges carry no length (all `None`).
  fn edge_branch_lengths<D: Send + Sync>(graph: &Graph<D>) -> BTreeMap<GraphEdgeKey, Option<f64>> {
    graph
      .get_edges()
      .iter()
      .map(|edge| (edge.read_arc().key(), None))
      .collect()
  }

  #[derive(Debug)]
  struct NonCloneData;
}
