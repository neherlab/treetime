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
    // Tree `((A,B)AB,C)root` with distinct own-values per node:
    //   A=1, B=2, C=3, AB=10, root=100.
    // Each node returns NodeOut = own-value + sum of children NodeOut, and sends its NodeOut up as
    // the message on its own parent edge. Derived by hand from the tree:
    //   A=1, B=2, C=3 (leaves), AB=10+1+2=13, root=100+13+3=116 (= total of all own-values).
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
    // Publication must be race-free: identical outputs under a 1-thread and a 4-thread rayon pool.
    let (graph, names) = fixture_tree()?;

    let single = run_backward_sum(&graph, &names, 1)?;
    let multi = run_backward_sum(&graph, &names, 4)?;

    assert_eq!(single.nodes, multi.nodes);
    assert_eq!(single.edges, multi.edges);
    Ok(())
  }

  #[test]
  fn test_pass_map_backward_without_edge_messages_collects_nodes_only() -> Result<(), Report> {
    // A node-only backward map: every visitor returns `parent_message: None` (here `EdgeOut = ()`), so
    // no node emits an upward edge message. The map must complete without panic, collect every node
    // output, and leave the per-edge map empty because no message travelled any edge.
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
    // Root-to-leaf prefix sums over `((A,B)AB,C)root` with own-values A=1, B=2, C=3, AB=10, root=100:
    //   root=100, AB=110, A=111, B=112, C=103.
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
    // The value engine must hand each visitor its children in the graph's canonical `children_of`
    // (outbound-edge) order, which the ordered reductions depend on. The fixture is built so that
    // outbound order differs from node-key order, so a pass that sorted children by key would fail.
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
    // The two children are not in ascending key order, so the check is not vacuous.
    assert_ne!(expected, {
      let mut sorted = expected.clone();
      sorted.sort_unstable();
      sorted
    });
    Ok(())
  }

  #[test]
  fn test_pass_map_backward_chain_propagates_leaf_to_root() -> Result<(), Report> {
    // A linear chain `root -> mid -> leaf`. Backward subtree sums with own-values leaf=1, mid=10,
    // root=100 give leaf=1, mid=11, root=111.
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
    // A single node that is both root and leaf: no children, no parent edge, no messages.
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
    // A removed node leaves a gap in the graph's node-key space (keys are never reused). The pass keys
    // its topology by actual node key, so a non-contiguous key set must still map correctly.
    let mut graph = Graph::new();
    let root = graph.add_node();
    let child = graph.add_node();
    let removed = graph.add_node(); // takes a key in the middle
    graph.add_edge(root, child)?;
    graph.remove_node(removed)?; // leaves a gap at the removed node's key
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
    // A map whose visitor returns its input unchanged and its parent edge unchanged reproduces the
    // input maps exactly (a value round-trip, replacing the in-place engine's restore round-trip).
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
  fn test_pass_map_backward_failure_leaves_inputs_unchanged_and_retries() -> Result<(), Report> {
    // A failing visitor returns the error, publishes no partial result, and leaves the borrowed input
    // maps untouched, so the caller can retry from the same inputs and succeed.
    let (graph, names) = fixture_tree()?;
    let (nodes, edges) = own_value_pass_values(&graph, &names);
    let nodes_before = nodes.clone();
    let edges_before = edges.clone();
    let pass = GraphPass::new(&graph)?;

    let failed: Result<GraphMapOutputs<usize, usize>, Report> = pass.map_backward(
      &nodes,
      &edges,
      |_| Ok(0),
      |_| Err(make_report!("injected pass failure")),
    );
    assert_error!(failed, "injected pass failure");

    // Borrowed inputs are never mutated by a pass, so a failure leaves them exactly as supplied.
    assert_eq!(nodes_before, nodes);
    assert_eq!(edges_before, edges);

    // Retry from the same inputs succeeds and produces the expected subtree sums.
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
  fn test_pass_map_backward_failing_child_blocks_ancestors_not_sibling() -> Result<(), Report> {
    // Tree `((A,B)AB,C)root`. Leaf B fails. Its parent AB depends on B, so AB never becomes ready and
    // is never visited; root depends on AB, so it is never visited either. The ready sibling A and the
    // independent leaf C may run. The error is returned and the inputs are left unchanged.
    let (graph, names) = fixture_tree()?;
    let (nodes, edges) = own_value_pass_values(&graph, &names);
    let nodes_before = nodes.clone();
    let edges_before = edges.clone();
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
    assert_eq!(nodes_before, nodes);
    assert_eq!(edges_before, edges);
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

    /// Node names threaded as a value map keyed by node key.
    pub type Names = BTreeMap<GraphNodeKey, String>;

    /// `((A,B)AB,C)root`. Returns the graph and the node-name value map, keyed by node key in
    /// creation order (`root`, `AB`, `A`, `B`, `C`).
    pub fn fixture_tree() -> Result<(Graph, Names), Report> {
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

    /// Linear chain `root -> mid -> leaf`. Returns the graph and keys `[root, mid, leaf]`.
    pub fn fixture_chain() -> Result<(Graph, Vec<GraphNodeKey>), Report> {
      let mut graph = Graph::new();
      let root = graph.add_node();
      let mid = graph.add_node();
      let leaf = graph.add_node();
      graph.add_edge(root, mid)?;
      graph.add_edge(mid, leaf)?;
      graph.build()?;
      Ok((graph, vec![root, mid, leaf]))
    }

    /// A single parent with two children whose outbound-edge order is the reverse of their key order:
    /// the parent's second-added child has the lower key. Returns the graph and the parent key.
    pub fn fixture_ordering() -> Result<(Graph, GraphNodeKey), Report> {
      let mut graph = Graph::new();
      let parent = graph.add_node();
      let first_child = graph.add_node();
      let second_child = graph.add_node();
      // Add the higher-key child's edge first, then the lower-key child's, so `children_of` order
      // (outbound order: [second_child, first_child]) differs from ascending key order.
      graph.add_edge(parent, second_child)?;
      graph.add_edge(parent, first_child)?;
      graph.build()?;
      Ok((graph, parent))
    }

    /// The child node keys of `parent` in the graph's `children_of` (outbound-edge) order.
    pub fn child_order_by_parent(graph: &Graph, parent: GraphNodeKey) -> Vec<GraphNodeKey> {
      let node = graph.get_node(parent).expect("Parent must exist");
      graph.children_of(node).map(|(child, _)| child.key()).collect()
    }

    /// Pass inputs with a distinct own-value per node (by name) and zero edge inputs.
    pub fn own_value_pass_values(
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

    /// Run the value-returning backward map on a pool of `threads` workers, computing each node's
    /// subtree sum and sending it up as the parent-edge message.
    pub fn run_backward_sum(
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

    /// Run the value-returning forward map on a pool of `threads` workers, computing each node's
    /// root-to-leaf prefix sum and sending it down as the parent-edge message.
    pub fn run_forward_sum(
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

    /// Map per-edge outputs to the name of the child node the edge points to.
    pub fn edge_values_by_child_name(
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

    pub fn values_by_name(names: &Names, values: &BTreeMap<GraphNodeKey, usize>) -> BTreeMap<String, usize> {
      values.iter().map(|(key, value)| (names[key].clone(), *value)).collect()
    }
  }
}
