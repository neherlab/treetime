#[cfg(test)]
mod tests {
  use crate::graph::Graph;
  use crate::node::GraphNodeKey;
  use crate::pass::{GraphMapOutputs, GraphPass, GraphPassNodeOutput};
  use eyre::Report;
  use maplit::btreemap;
  use parking_lot::Mutex;
  use pretty_assertions::assert_eq;
  use std::collections::{BTreeMap, BTreeSet};
  use treetime_utils::{assert_error, make_report, o};

  use self::helpers::{
    child_order_by_parent, edge_values_by_child_name, fixture_chain, fixture_ordering, fixture_tree,
    own_value_pass_values, run_backward_sum, run_forward_sum, values_by_name,
  };

  #[test]
  fn test_pass_map_backward_collects_returned_subtree_sums() -> Result<(), Report> {
    let (graph, names) = fixture_tree()?;
    let outputs = run_backward_sum(&graph, &names, 4)?;

    let actual_nodes = values_by_name(&names, &outputs.nodes);
    let expected_nodes = btreemap! {
      o!("A") => 1,
      o!("B") => 2,
      o!("C") => 3,
      o!("AB") => 13,
      o!("root") => 116,
    };
    assert_eq!(expected_nodes, actual_nodes);

    let actual_edges = edge_values_by_child_name(&graph, &names, &outputs.edges)?;
    let expected_edges = btreemap! {
      o!("A") => 1,
      o!("B") => 2,
      o!("C") => 3,
      o!("AB") => 13,
    };
    assert_eq!(expected_edges, actual_edges);

    Ok(())
  }

  #[test]
  fn test_pass_map_backward_is_thread_count_independent() -> Result<(), Report> {
    let (graph, names) = fixture_tree()?;

    let single = run_backward_sum(&graph, &names, 1)?;
    let multi = run_backward_sum(&graph, &names, 4)?;

    assert_eq!(single.nodes, multi.nodes);
    assert_eq!(single.edges, multi.edges);
    Ok(())
  }

  #[test]
  fn test_pass_map_backward_without_edge_messages_collects_nodes_only() -> Result<(), Report> {
    let (graph, names) = fixture_tree()?;
    let (nodes, edges) = own_value_pass_values(&graph, &names);
    let pass = GraphPass::new(&graph)?;

    let outputs: GraphMapOutputs<usize, ()> = pass.map_backward(
      &nodes,
      &edges,
      |_| Ok(0),
      |context| {
        let children_sum = context.children.iter().map(|child| *child.node).sum::<usize>();
        Ok(GraphPassNodeOutput {
          node: *context.input + children_sum,
          parent_message: None,
        })
      },
    )?;

    assert!(outputs.edges.is_empty());
    let actual_nodes = values_by_name(&names, &outputs.nodes);
    let expected_nodes = btreemap! {
      o!("A") => 1,
      o!("B") => 2,
      o!("C") => 3,
      o!("AB") => 13,
      o!("root") => 116,
    };
    assert_eq!(expected_nodes, actual_nodes);
    Ok(())
  }

  #[test]
  fn test_pass_map_forward_accumulates_root_to_leaf() -> Result<(), Report> {
    let (graph, names) = fixture_tree()?;
    let outputs = run_forward_sum(&graph, &names, 4)?;

    let actual_nodes = values_by_name(&names, &outputs.nodes);
    let expected_nodes = btreemap! {
      o!("root") => 100,
      o!("AB") => 110,
      o!("A") => 111,
      o!("B") => 112,
      o!("C") => 103,
    };
    assert_eq!(expected_nodes, actual_nodes);

    let actual_edges = edge_values_by_child_name(&graph, &names, &outputs.edges)?;
    let expected_edges = btreemap! {
      o!("AB") => 110,
      o!("A") => 111,
      o!("B") => 112,
      o!("C") => 103,
    };
    assert_eq!(expected_edges, actual_edges);

    Ok(())
  }

  #[test]
  fn test_pass_map_forward_is_thread_count_independent() -> Result<(), Report> {
    let (graph, names) = fixture_tree()?;

    let single = run_forward_sum(&graph, &names, 1)?;
    let multi = run_forward_sum(&graph, &names, 4)?;

    assert_eq!(single.nodes, multi.nodes);
    assert_eq!(single.edges, multi.edges);
    Ok(())
  }

  #[test]
  fn test_pass_map_backward_children_arrive_in_children_of_order() -> Result<(), Report> {
    let (graph, parent_key) = fixture_ordering()?;
    let (nodes, edges) = pass_zeros(&graph);

    let seen: Mutex<BTreeMap<GraphNodeKey, Vec<GraphNodeKey>>> = Mutex::new(BTreeMap::new());
    let _outputs: GraphMapOutputs<usize, ()> = GraphPass::new(&graph)?.map_backward(
      &nodes,
      &edges,
      |_| Ok(0),
      |context| {
        let order = context.children.iter().map(|child| child.node_key).collect::<Vec<_>>();
        seen.lock().insert(context.key, order);
        Ok(GraphPassNodeOutput {
          node: 0,
          parent_message: None,
        })
      },
    )?;

    let actual = seen.lock().get(&parent_key).cloned().expect("Parent must be visited");
    let expected = child_order_by_parent(&graph, parent_key);
    assert_eq!(expected, actual);
    assert_ne!(expected, {
      let mut sorted = expected.clone();
      sorted.sort_unstable();
      sorted
    });
    Ok(())
  }

  #[test]
  fn test_pass_map_backward_chain_propagates_leaf_to_root() -> Result<(), Report> {
    let (graph, keys) = fixture_chain()?;
    let nodes = keys.iter().copied().zip([100, 10, 1]).collect::<BTreeMap<_, _>>();
    let edges = zero_edges(&graph);

    let outputs: GraphMapOutputs<usize, usize> = GraphPass::new(&graph)?.map_backward(
      &nodes,
      &edges,
      |_| Ok(0),
      |context| {
        let children_sum = context.children.iter().map(|child| *child.node).sum::<usize>();
        let node = *context.input + children_sum;
        let parent_message = (!context.is_root).then_some(node);
        Ok(GraphPassNodeOutput { node, parent_message })
      },
    )?;

    let expected = btreemap! { keys[0] => 111, keys[1] => 11, keys[2] => 1 };
    assert_eq!(expected, outputs.nodes);
    Ok(())
  }

  #[test]
  fn test_pass_map_backward_single_root_only_node() -> Result<(), Report> {
    let mut graph = Graph::new();
    let root = graph.add_node();
    graph.build()?;
    let nodes = btreemap! { root => 7 };
    let edges = zero_edges(&graph);

    let outputs: GraphMapOutputs<usize, usize> = GraphPass::new(&graph)?.map_backward(
      &nodes,
      &edges,
      |_| Ok(0),
      |context| {
        assert!(context.is_root && context.is_leaf && context.children.is_empty() && context.parent_edge.is_none());
        Ok(GraphPassNodeOutput {
          node: *context.input,
          parent_message: None,
        })
      },
    )?;

    assert_eq!(btreemap! { root => 7 }, outputs.nodes);
    assert!(outputs.edges.is_empty());
    Ok(())
  }

  #[test]
  fn test_pass_map_backward_handles_deleted_key_gaps() -> Result<(), Report> {
    let mut graph = Graph::new();
    let root = graph.add_node();
    let child = graph.add_node();
    let removed = graph.add_node();
    graph.add_edge(root, child)?;
    graph.remove_node(removed)?;
    graph.build()?;

    assert!(
      removed.as_usize() > child.as_usize(),
      "the removed key sits between live keys"
    );
    let nodes = btreemap! { root => 100, child => 1 };
    let edges = zero_edges(&graph);

    let outputs: GraphMapOutputs<usize, usize> = GraphPass::new(&graph)?.map_backward(
      &nodes,
      &edges,
      |_| Ok(0),
      |context| {
        let children_sum = context.children.iter().map(|child| *child.node).sum::<usize>();
        let node = *context.input + children_sum;
        let parent_message = (!context.is_root).then_some(node);
        Ok(GraphPassNodeOutput { node, parent_message })
      },
    )?;

    assert_eq!(btreemap! { root => 101, child => 1 }, outputs.nodes);
    Ok(())
  }

  #[test]
  fn test_pass_map_identity_preserves_all_inputs() -> Result<(), Report> {
    let (graph, _names) = fixture_tree()?;
    let (nodes, edges) = key_indices(&graph);

    let outputs: GraphMapOutputs<usize, usize> = GraphPass::new(&graph)?.map_backward(
      &nodes,
      &edges,
      |_| Ok(0),
      |context| {
        let parent_message = context.parent_edge.map(|(_, edge)| *edge);
        Ok(GraphPassNodeOutput {
          node: *context.input,
          parent_message,
        })
      },
    )?;

    assert_eq!(nodes, outputs.nodes);
    assert_eq!(edges, outputs.edges);
    Ok(())
  }

  #[test]
  fn test_pass_map_backward_retry_after_failure_succeeds() -> Result<(), Report> {
    let (graph, names) = fixture_tree()?;
    let (nodes, edges) = own_value_pass_values(&graph, &names);
    let pass = GraphPass::new(&graph)?;

    let failed: Result<GraphMapOutputs<usize, usize>, Report> = pass.map_backward(
      &nodes,
      &edges,
      |_| Ok(0),
      |_| Err(make_report!("injected pass failure")),
    );
    assert_error!(failed, "injected pass failure");

    let outputs = pass.map_backward(
      &nodes,
      &edges,
      |_| Ok(0),
      |context| {
        let children_sum = context.children.iter().map(|child| *child.node).sum::<usize>();
        let node = *context.input + children_sum;
        let parent_message = (!context.is_root).then_some(node);
        Ok(GraphPassNodeOutput { node, parent_message })
      },
    )?;
    let actual_nodes = values_by_name(&names, &outputs.nodes);
    assert_eq!(&116, &actual_nodes[&o!("root")]);
    Ok(())
  }

  #[test]
  fn test_pass_map_backward_failing_child_blocks_ancestors() -> Result<(), Report> {
    let (graph, names) = fixture_tree()?;
    let (nodes, edges) = own_value_pass_values(&graph, &names);
    let key_by_name = names
      .iter()
      .map(|(key, name)| (name.clone(), *key))
      .collect::<BTreeMap<_, _>>();

    let visited: Mutex<BTreeSet<GraphNodeKey>> = Mutex::new(BTreeSet::new());
    let failed: Result<GraphMapOutputs<usize, usize>, Report> = GraphPass::new(&graph)?.map_backward(
      &nodes,
      &edges,
      |_| Ok(0),
      |context| {
        visited.lock().insert(context.key);
        if context.key == key_by_name[&o!("B")] {
          return Err(make_report!("injected child failure"));
        }
        Ok(GraphPassNodeOutput {
          node: *context.input,
          parent_message: (!context.is_root).then_some(0),
        })
      },
    );
    assert_error!(failed, "injected child failure");

    let visited = visited.lock();
    assert!(
      !visited.contains(&key_by_name[&o!("AB")]),
      "the failing child's parent must not run"
    );
    assert!(
      !visited.contains(&key_by_name[&o!("root")]),
      "an ancestor of the failing child must not run"
    );
    Ok(())
  }

  fn pass_zeros(
    graph: &Graph,
  ) -> (
    BTreeMap<GraphNodeKey, usize>,
    BTreeMap<crate::edge::GraphEdgeKey, usize>,
  ) {
    (
      graph.get_nodes().map(|node| (node.key(), 0)).collect(),
      zero_edges(graph),
    )
  }

  fn zero_edges(graph: &Graph) -> BTreeMap<crate::edge::GraphEdgeKey, usize> {
    graph.get_edges().map(|edge| (edge.key(), 0)).collect()
  }

  fn key_indices(
    graph: &Graph,
  ) -> (
    BTreeMap<GraphNodeKey, usize>,
    BTreeMap<crate::edge::GraphEdgeKey, usize>,
  ) {
    let nodes = graph
      .get_nodes()
      .map(|node| {
        let key = node.key();
        (key, key.as_usize())
      })
      .collect();
    let edges = graph
      .get_edges()
      .map(|edge| {
        let key = edge.key();
        (key, key.as_usize())
      })
      .collect();
    (nodes, edges)
  }

  mod helpers {
    use crate::edge::GraphEdgeKey;
    use crate::graph::Graph;
    use crate::node::GraphNodeKey;
    use crate::pass::{GraphMapOutputs, GraphPass, GraphPassNodeOutput};
    use eyre::Report;
    use maplit::btreemap;
    use rayon::ThreadPoolBuilder;
    use std::collections::BTreeMap;
    use treetime_utils::o;

    pub(super) type Names = BTreeMap<GraphNodeKey, String>;

    pub(super) fn fixture_tree() -> Result<(Graph, Names), Report> {
      let mut graph = Graph::new();
      let root = graph.add_node();
      let ab = graph.add_node();
      let tip_a = graph.add_node();
      let tip_b = graph.add_node();
      let tip_c = graph.add_node();

      graph.add_edge(root, ab)?;
      graph.add_edge(root, tip_c)?;
      graph.add_edge(ab, tip_a)?;
      graph.add_edge(ab, tip_b)?;
      graph.build()?;

      let names = btreemap! {
        root => o!("root"),
        ab => o!("AB"),
        tip_a => o!("A"),
        tip_b => o!("B"),
        tip_c => o!("C"),
      };

      Ok((graph, names))
    }

    pub(super) fn fixture_chain() -> Result<(Graph, Vec<GraphNodeKey>), Report> {
      let mut graph = Graph::new();
      let root = graph.add_node();
      let mid = graph.add_node();
      let leaf = graph.add_node();
      graph.add_edge(root, mid)?;
      graph.add_edge(mid, leaf)?;
      graph.build()?;
      Ok((graph, vec![root, mid, leaf]))
    }

    pub(super) fn fixture_ordering() -> Result<(Graph, GraphNodeKey), Report> {
      let mut graph = Graph::new();
      let parent = graph.add_node();
      let first_child = graph.add_node();
      let second_child = graph.add_node();
      graph.add_edge(parent, second_child)?;
      graph.add_edge(parent, first_child)?;
      graph.build()?;
      Ok((graph, parent))
    }

    pub(super) fn child_order_by_parent(graph: &Graph, parent: GraphNodeKey) -> Vec<GraphNodeKey> {
      let node = graph.get_node(parent).expect("Parent must exist");
      graph.children_of(node).map(|(child, _)| child.key()).collect()
    }

    pub(super) fn own_value_pass_values(
      graph: &Graph,
      names: &Names,
    ) -> (BTreeMap<GraphNodeKey, usize>, BTreeMap<GraphEdgeKey, usize>) {
      let by_name = btreemap! {
        o!("A") => 1,
        o!("B") => 2,
        o!("C") => 3,
        o!("AB") => 10,
        o!("root") => 100,
      };
      let nodes = graph
        .get_nodes()
        .map(|node| {
          let key = node.key();
          (key, by_name[&names[&key]])
        })
        .collect();
      let edges = graph.get_edges().map(|edge| (edge.key(), 0)).collect();
      (nodes, edges)
    }

    pub(super) fn run_backward_sum(
      graph: &Graph,
      names: &Names,
      threads: usize,
    ) -> Result<GraphMapOutputs<usize, usize>, Report> {
      let (nodes, edges) = own_value_pass_values(graph, names);
      let pass = GraphPass::new(graph)?;
      let pool = ThreadPoolBuilder::new().num_threads(threads).build()?;
      pool.install(|| {
        pass.map_backward(
          &nodes,
          &edges,
          |_| Ok(0),
          |context| {
            let children_sum = context.children.iter().map(|child| *child.node).sum::<usize>();
            let node = *context.input + children_sum;
            let parent_message = (!context.is_root).then_some(node);
            Ok(GraphPassNodeOutput { node, parent_message })
          },
        )
      })
    }

    pub(super) fn run_forward_sum(
      graph: &Graph,
      names: &Names,
      threads: usize,
    ) -> Result<GraphMapOutputs<usize, usize>, Report> {
      let (nodes, edges) = own_value_pass_values(graph, names);
      let pass = GraphPass::new(graph)?;
      let pool = ThreadPoolBuilder::new().num_threads(threads).build()?;
      pool.install(|| {
        pass.map_forward(
          &nodes,
          &edges,
          |_| Ok(0),
          |context| {
            let parent_sum = context.parent.copied().unwrap_or(0);
            let node = *context.input + parent_sum;
            let parent_message = (!context.is_root).then_some(node);
            Ok(GraphPassNodeOutput { node, parent_message })
          },
        )
      })
    }

    pub(super) fn edge_values_by_child_name(
      graph: &Graph,
      names: &Names,
      values: &BTreeMap<GraphEdgeKey, usize>,
    ) -> Result<BTreeMap<String, usize>, Report> {
      values
        .iter()
        .map(|(edge_key, value)| {
          let child_key = graph.get_target_node_key(*edge_key)?;
          Ok((names[&child_key].clone(), *value))
        })
        .collect()
    }

    pub(super) fn values_by_name(names: &Names, values: &BTreeMap<GraphNodeKey, usize>) -> BTreeMap<String, usize> {
      values.iter().map(|(key, value)| (names[key].clone(), *value)).collect()
    }
  }
}
