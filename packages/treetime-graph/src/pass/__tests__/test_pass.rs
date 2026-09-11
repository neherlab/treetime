#[cfg(test)]
mod tests {
  use crate::pass::{GraphMapOutputs, GraphPass, GraphPassNodeOutput};
  use eyre::Report;
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use treetime_utils::{assert_error, make_report, o};

  use self::helpers::{
    edge_values_by_child_name, fixture_tree, key_payloads, own_value_pass_values, pass_values, run_backward_sum,
    run_forward_sum, values_by_name,
  };

  #[test]
  fn test_pass_backward_visits_children_before_parent() -> Result<(), Report> {
    let graph = fixture_tree()?;
    let (mut nodes, mut edges) = pass_values(&graph);
    let mut pass = GraphPass::new(&graph, &mut nodes, &mut edges, |_| Ok(0))?;

    pass.try_for_each_backward(|dependencies, slot| {
      let graph_node = graph.get_node(slot.key).expect("Indexed node must exist");
      slot.node = graph
        .children_of(&graph_node.read_arc())
        .iter()
        .map(|(child, _)| dependencies.node(child.read_arc().key()))
        .sum::<usize>()
        + 1;
      Ok(())
    })?;
    let (nodes, _) = pass.into_maps()?;
    let actual = values_by_name(&graph, &nodes);
    let expected = btreemap! {
      o!("A") => 1,
      o!("AB") => 3,
      o!("B") => 1,
      o!("C") => 1,
      o!("root") => 5,
    };

    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  fn test_pass_map_backward_collects_returned_subtree_sums() -> Result<(), Report> {
    // Tree `((A,B)AB,C)root` with distinct own-values per node:
    //   A=1, B=2, C=3, AB=10, root=100.
    // Each node returns NodeOut = own-value + sum of children NodeOut, and sends its NodeOut up as
    // the message on its own parent edge. Derived by hand from the tree:
    //   A=1, B=2, C=3 (leaves), AB=10+1+2=13, root=100+13+3=116 (= total of all own-values).
    let graph = fixture_tree()?;
    let outputs = run_backward_sum(&graph, 4)?;

    let actual_nodes = values_by_name(&graph, &outputs.nodes);
    let expected_nodes = btreemap! {
      o!("A") => 1,
      o!("B") => 2,
      o!("C") => 3,
      o!("AB") => 13,
      o!("root") => 116,
    };
    assert_eq!(expected_nodes, actual_nodes);

    // The root's total (116) and the interior AB subtree sum (13) are the load-bearing checks.
    assert_eq!(&116, &actual_nodes[&o!("root")]);
    assert_eq!(&13, &actual_nodes[&o!("AB")]);

    // Each edge carries the child's upward message, equal to that child's NodeOut. The root has no
    // parent edge, so it contributes no message.
    let actual_edges = edge_values_by_child_name(&graph, &outputs.edges)?;
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
    let graph = fixture_tree()?;

    let single = run_backward_sum(&graph, 1)?;
    let multi = run_backward_sum(&graph, 4)?;

    assert_eq!(single.nodes, multi.nodes);
    assert_eq!(single.edges, multi.edges);
    Ok(())
  }

  #[test]
  fn test_pass_map_backward_without_edge_messages_collects_nodes_only() -> Result<(), Report> {
    // A node-only backward map: every visitor returns `parent_message: None` (here `EdgeOut = ()`), so
    // no node emits an upward edge message. The map must complete without panic, collect every node
    // output, and leave the per-edge map empty because no message travelled any edge.
    let graph = fixture_tree()?;
    let (mut nodes, mut edges) = own_value_pass_values(&graph);
    let pass = GraphPass::new(&graph, &mut nodes, &mut edges, |_| Ok(0))?;

    let outputs: GraphMapOutputs<usize, ()> = pass.try_map_backward(|context| {
      let children_sum = context.children.iter().map(|child| *child.node).sum::<usize>();
      Ok(GraphPassNodeOutput {
        node: context.input + children_sum,
        parent_message: None,
      })
    })?;

    assert!(outputs.edges.is_empty());
    let actual_nodes = values_by_name(&graph, &outputs.nodes);
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
    // Tree `((A,B)AB,C)root` with distinct own-values per node:
    //   A=1, B=2, C=3, AB=10, root=100.
    // Each node returns NodeOut = own-value + parent's NodeOut (0 at the root), and sends its own
    // NodeOut down as the message on its own parent edge. These are root-to-leaf prefix sums,
    // derived by hand from the tree:
    //   root=100, AB=100+10=110, A=110+1=111, B=110+2=112, C=100+3=103.
    let graph = fixture_tree()?;
    let outputs = run_forward_sum(&graph, 4)?;

    let actual_nodes = values_by_name(&graph, &outputs.nodes);
    let expected_nodes = btreemap! {
      o!("root") => 100,
      o!("AB") => 110,
      o!("A") => 111,
      o!("B") => 112,
      o!("C") => 103,
    };
    assert_eq!(expected_nodes, actual_nodes);

    // Each edge carries the child's downward message, equal to that child's NodeOut. The root has no
    // parent edge, so it contributes no message.
    let actual_edges = edge_values_by_child_name(&graph, &outputs.edges)?;
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
    // Publication must be race-free: identical outputs under a 1-thread and a 4-thread rayon pool.
    let graph = fixture_tree()?;

    let single = run_forward_sum(&graph, 1)?;
    let multi = run_forward_sum(&graph, 4)?;

    assert_eq!(single.nodes, multi.nodes);
    assert_eq!(single.edges, multi.edges);
    Ok(())
  }

  #[test]
  fn test_pass_forward_visits_parent_before_children() -> Result<(), Report> {
    let graph = fixture_tree()?;
    let (mut nodes, mut edges) = pass_values(&graph);
    let mut pass = GraphPass::new(&graph, &mut nodes, &mut edges, |_| Ok(0))?;

    pass.try_for_each_forward(|dependencies, slot| {
      slot.node = slot.parent_key.map_or(0, |parent| dependencies.node(parent) + 1);
      Ok(())
    })?;
    let (nodes, _) = pass.into_maps()?;
    let actual = values_by_name(&graph, &nodes);
    let expected = btreemap! {
      o!("A") => 2,
      o!("AB") => 1,
      o!("B") => 2,
      o!("C") => 1,
      o!("root") => 0,
    };

    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  fn test_pass_roundtrip_preserves_all_values() -> Result<(), Report> {
    let graph = fixture_tree()?;
    let (mut nodes, mut edges) = key_payloads(&graph);
    let expected_nodes = nodes.clone();
    let expected_edges = edges.clone();

    let pass = GraphPass::new(&graph, &mut nodes, &mut edges, |_| {
      unreachable!("all graph nodes are present")
    })?;
    let (actual_nodes, actual_edges) = pass.into_maps()?;

    assert_eq!(expected_nodes, actual_nodes);
    assert_eq!(expected_edges, actual_edges);
    Ok(())
  }

  #[test]
  fn test_pass_error_restores_all_values() -> Result<(), Report> {
    let graph = fixture_tree()?;
    let (mut nodes, mut edges) = key_payloads(&graph);
    let expected_nodes = nodes.clone();
    let expected_edges = edges.clone();
    let mut pass = GraphPass::new(&graph, &mut nodes, &mut edges, |_| {
      unreachable!("all graph nodes are present")
    })?;

    let result = pass.try_for_each_backward(|_, _| Err(make_report!("injected pass failure")));
    assert_error!(result, "injected pass failure");
    let (actual_nodes, actual_edges) = pass.into_maps()?;

    assert_eq!(expected_nodes, actual_nodes);
    assert_eq!(expected_edges, actual_edges);
    Ok(())
  }

  mod helpers {
    use crate::edge::{GraphEdge, GraphEdgeKey};
    use crate::graph::Graph;
    use crate::node::{GraphNode, GraphNodeKey};
    use crate::pass::{GraphMapOutputs, GraphPass, GraphPassNodeOutput};
    use eyre::Report;
    use maplit::btreemap;
    use rayon::ThreadPoolBuilder;
    use std::collections::BTreeMap;
    use treetime_utils::o;

    /// `((A,B)AB,C)root` with branch lengths on every edge.
    pub fn fixture_tree() -> Result<Graph<TestNode, TestEdge, ()>, Report> {
      let mut graph = Graph::<TestNode, TestEdge, ()>::new();
      let root = graph.add_node(TestNode::new("root"));
      let ab = graph.add_node(TestNode::new("AB"));
      let tip_a = graph.add_node(TestNode::new("A"));
      let tip_b = graph.add_node(TestNode::new("B"));
      let tip_c = graph.add_node(TestNode::new("C"));

      graph.add_edge(root, ab, TestEdge::with_length(3.0))?;
      graph.add_edge(root, tip_c, TestEdge::with_length(4.0))?;
      graph.add_edge(ab, tip_a, TestEdge::with_length(1.0))?;
      graph.add_edge(ab, tip_b, TestEdge::with_length(2.0))?;
      graph.build()?;

      Ok(graph)
    }

    /// Zero-initialized pass payloads for every node and edge.
    pub fn pass_values(
      graph: &Graph<TestNode, TestEdge, ()>,
    ) -> (BTreeMap<GraphNodeKey, usize>, BTreeMap<GraphEdgeKey, usize>) {
      let nodes = graph
        .get_nodes()
        .iter()
        .map(|node| (node.read_arc().key(), 0))
        .collect();
      let edges = graph
        .get_edges()
        .iter()
        .map(|edge| (edge.read_arc().key(), 0))
        .collect();
      (nodes, edges)
    }

    /// Pass payloads with a distinct own-value per node (by name) and zero edge inputs.
    pub fn own_value_pass_values(
      graph: &Graph<TestNode, TestEdge, ()>,
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
        .iter()
        .map(|node| {
          let node = node.read_arc();
          (node.key(), by_name[&node.payload().read_arc().0])
        })
        .collect();
      let edges = graph
        .get_edges()
        .iter()
        .map(|edge| (edge.read_arc().key(), 0))
        .collect();
      (nodes, edges)
    }

    /// Run the value-returning backward map on a pool of `threads` workers, computing each node's
    /// subtree sum and sending it up as the parent-edge message.
    pub fn run_backward_sum(
      graph: &Graph<TestNode, TestEdge, ()>,
      threads: usize,
    ) -> Result<GraphMapOutputs<usize, usize>, Report> {
      let (mut nodes, mut edges) = own_value_pass_values(graph);
      let pass = GraphPass::new(graph, &mut nodes, &mut edges, |_| Ok(0))?;
      let pool = ThreadPoolBuilder::new().num_threads(threads).build()?;
      pool.install(|| {
        pass.try_map_backward(|context| {
          let children_sum = context.children.iter().map(|child| *child.node).sum::<usize>();
          let node = context.input + children_sum;
          let parent_message = (!context.is_root).then_some(node);
          Ok(GraphPassNodeOutput { node, parent_message })
        })
      })
    }

    /// Run the value-returning forward map on a pool of `threads` workers, computing each node's
    /// root-to-leaf prefix sum and sending it down as the parent-edge message.
    pub fn run_forward_sum(
      graph: &Graph<TestNode, TestEdge, ()>,
      threads: usize,
    ) -> Result<GraphMapOutputs<usize, usize>, Report> {
      let (mut nodes, mut edges) = own_value_pass_values(graph);
      let pass = GraphPass::new(graph, &mut nodes, &mut edges, |_| Ok(0))?;
      let pool = ThreadPoolBuilder::new().num_threads(threads).build()?;
      pool.install(|| {
        pass.try_map_forward(|context| {
          let parent_sum = context.parent.copied().unwrap_or(0);
          let node = context.input + parent_sum;
          let parent_message = (!context.is_root).then_some(node);
          Ok(GraphPassNodeOutput { node, parent_message })
        })
      })
    }

    /// Map per-edge outputs to the name of the child node the edge points to.
    pub fn edge_values_by_child_name(
      graph: &Graph<TestNode, TestEdge, ()>,
      values: &BTreeMap<GraphEdgeKey, usize>,
    ) -> Result<BTreeMap<String, usize>, Report> {
      values
        .iter()
        .map(|(edge_key, value)| {
          let child_key = graph.get_target_node_key(*edge_key)?;
          let child = graph.get_node(child_key).expect("Indexed child node must exist");
          let name = child.read_arc().payload().read_arc().0.clone();
          Ok((name, *value))
        })
        .collect()
    }

    /// Pass payloads keyed and valued by the underlying key index, for round-trip checks.
    pub fn key_payloads(
      graph: &Graph<TestNode, TestEdge, ()>,
    ) -> (BTreeMap<GraphNodeKey, usize>, BTreeMap<GraphEdgeKey, usize>) {
      let nodes = graph
        .get_nodes()
        .iter()
        .map(|node| {
          let key = node.read_arc().key();
          (key, key.as_usize())
        })
        .collect();
      let edges = graph
        .get_edges()
        .iter()
        .map(|edge| {
          let key = edge.read_arc().key();
          (key, key.as_usize())
        })
        .collect();
      (nodes, edges)
    }

    pub fn values_by_name(
      graph: &Graph<TestNode, TestEdge, ()>,
      values: &BTreeMap<GraphNodeKey, usize>,
    ) -> BTreeMap<String, usize> {
      graph
        .get_nodes()
        .iter()
        .map(|node| {
          let node = node.read_arc();
          (node.payload().read_arc().0.clone(), values[&node.key()])
        })
        .collect()
    }

    #[derive(Debug, Default, Eq, PartialEq)]
    pub struct TestNode(pub String);

    impl TestNode {
      pub fn new(name: &str) -> Self {
        Self(name.to_owned())
      }
    }

    impl GraphNode for TestNode {}

    #[derive(Debug, Default, PartialEq)]
    pub struct TestEdge {
      pub branch_length: Option<f64>,
    }

    impl TestEdge {
      pub fn with_length(len: f64) -> Self {
        Self {
          branch_length: Some(len),
        }
      }
    }

    impl GraphEdge for TestEdge {}
  }
}
