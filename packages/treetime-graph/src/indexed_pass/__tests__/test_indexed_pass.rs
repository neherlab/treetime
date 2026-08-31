#[cfg(test)]
mod tests {
  use crate::indexed_pass::{IndexedPass, with_indexed_graph_payloads};
  use eyre::Report;
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use treetime_utils::{assert_error, make_report, o};

  use self::helpers::{edge_lengths, fixture_tree, key_payloads, node_names, pass_values, values_by_name};

  #[test]
  fn test_indexed_pass_backward_visits_children_before_parent() -> Result<(), Report> {
    let graph = fixture_tree()?;
    let (mut nodes, mut edges) = pass_values(&graph);
    let mut pass = IndexedPass::new(&graph, &mut nodes, &mut edges, |_| Ok(0))?;

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
  fn test_indexed_pass_forward_visits_parent_before_children() -> Result<(), Report> {
    let graph = fixture_tree()?;
    let (mut nodes, mut edges) = pass_values(&graph);
    let mut pass = IndexedPass::new(&graph, &mut nodes, &mut edges, |_| Ok(0))?;

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
  fn test_indexed_pass_roundtrip_preserves_all_values() -> Result<(), Report> {
    let graph = fixture_tree()?;
    let (mut nodes, mut edges) = key_payloads(&graph);
    let expected_nodes = nodes.clone();
    let expected_edges = edges.clone();

    let pass = IndexedPass::new(&graph, &mut nodes, &mut edges, |_| {
      unreachable!("all graph nodes are present")
    })?;
    let (actual_nodes, actual_edges) = pass.into_maps()?;

    assert_eq!(expected_nodes, actual_nodes);
    assert_eq!(expected_edges, actual_edges);
    Ok(())
  }

  #[test]
  fn test_indexed_pass_error_restores_all_values() -> Result<(), Report> {
    let graph = fixture_tree()?;
    let (mut nodes, mut edges) = key_payloads(&graph);
    let expected_nodes = nodes.clone();
    let expected_edges = edges.clone();
    let mut pass = IndexedPass::new(&graph, &mut nodes, &mut edges, |_| {
      unreachable!("all graph nodes are present")
    })?;

    let result = pass.try_for_each_backward(|_, _| Err(make_report!("injected pass failure")));
    assert_error!(result, "injected pass failure");
    let (actual_nodes, actual_edges) = pass.into_maps()?;

    assert_eq!(expected_nodes, actual_nodes);
    assert_eq!(expected_edges, actual_edges);
    Ok(())
  }

  #[test]
  fn test_indexed_graph_payloads_error_restores_graph() -> Result<(), Report> {
    let graph = fixture_tree()?;
    let expected_nodes = node_names(&graph);
    let expected_edges = edge_lengths(&graph);

    let result = with_indexed_graph_payloads(&graph, |_| Err::<(), _>(make_report!("injected graph pass failure")));
    assert_error!(result, "injected graph pass failure");
    let actual_nodes = node_names(&graph);
    let actual_edges = edge_lengths(&graph);

    assert_eq!(expected_nodes, actual_nodes);
    assert_eq!(expected_edges, actual_edges);
    Ok(())
  }

  mod helpers {
    use crate::edge::{GraphEdge, GraphEdgeKey};
    use crate::graph::Graph;
    use crate::node::{GraphNode, GraphNodeKey};
    use eyre::Report;
    use std::collections::BTreeMap;

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

    pub fn node_names(graph: &Graph<TestNode, TestEdge, ()>) -> BTreeMap<GraphNodeKey, String> {
      graph
        .get_nodes()
        .iter()
        .map(|node| {
          let node = node.read_arc();
          (node.key(), node.payload().read_arc().0.clone())
        })
        .collect()
    }

    pub fn edge_lengths(graph: &Graph<TestNode, TestEdge, ()>) -> BTreeMap<GraphEdgeKey, Option<f64>> {
      graph
        .get_edges()
        .iter()
        .map(|edge| {
          let edge = edge.read_arc();
          (edge.key(), edge.payload().read_arc().branch_length)
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
