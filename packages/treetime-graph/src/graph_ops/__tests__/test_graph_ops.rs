#[cfg(test)]
mod tests {
  use crate::edge::{GraphEdge, GraphEdgeKey};
  use crate::graph::Graph;
  use crate::node::{GraphNode, GraphNodeKey};
  use eyre::Report;
  use pretty_assertions::assert_eq;

  #[derive(Debug, PartialEq)]
  struct TestNode(&'static str);
  impl GraphNode for TestNode {}

  /// Carries a value so tests can tell a preserved payload from a rebuilt one.
  #[derive(Debug, PartialEq)]
  struct TestEdge(f64);
  impl GraphEdge for TestEdge {}

  type TestGraph = Graph<TestNode, TestEdge, ()>;

  /// `root -> {a, b}`, `a -> c`. Returns the graph and the keys in that order.
  fn fixture() -> Result<(TestGraph, [GraphNodeKey; 4], GraphEdgeKey), Report> {
    let mut graph = TestGraph::new();
    let root = graph.add_node(TestNode("root"));
    let a = graph.add_node(TestNode("a"));
    let b = graph.add_node(TestNode("b"));
    let c = graph.add_node(TestNode("c"));
    graph.add_edge(root, a, TestEdge(1.0))?;
    graph.add_edge(root, b, TestEdge(2.0))?;
    let a_to_c = graph.add_edge(a, c, TestEdge(3.0))?;
    graph.build()?;
    Ok((graph, [root, a, b, c], a_to_c))
  }

  fn outbound(graph: &TestGraph, node_key: GraphNodeKey) -> Vec<GraphEdgeKey> {
    graph
      .get_node(node_key)
      .expect("node exists")
      .read_arc()
      .outbound()
      .to_vec()
  }

  fn inbound(graph: &TestGraph, node_key: GraphNodeKey) -> Vec<GraphEdgeKey> {
    graph
      .get_node(node_key)
      .expect("node exists")
      .read_arc()
      .inbound()
      .to_vec()
  }

  #[test]
  fn test_reparent_edge_moves_the_edge_and_keeps_key_and_payload() -> Result<(), Report> {
    let (mut graph, [root, a, b, c], a_to_c) = fixture()?;

    graph.reparent_edge(a_to_c, b)?;

    let edge = graph.get_edge(a_to_c).expect("edge survives reparenting");
    assert_eq!(edge.read_arc().source(), b);
    assert_eq!(edge.read_arc().target(), c);
    assert_eq!(
      *edge.read_arc().payload().read_arc(),
      TestEdge(3.0),
      "payload must survive"
    );

    assert!(!outbound(&graph, a).contains(&a_to_c), "old source must drop the edge");
    assert!(outbound(&graph, b).contains(&a_to_c), "new source must gain the edge");
    assert_eq!(
      inbound(&graph, c),
      vec![a_to_c],
      "the target's inbound list is untouched"
    );
    assert_eq!(outbound(&graph, root).len(), 2, "unrelated nodes are untouched");
    Ok(())
  }

  #[test]
  fn test_reparent_edge_to_the_current_source_is_a_noop() -> Result<(), Report> {
    let (mut graph, [_, a, _, _], a_to_c) = fixture()?;
    let before = outbound(&graph, a);

    graph.reparent_edge(a_to_c, a)?;

    assert_eq!(outbound(&graph, a), before, "a no-op must not duplicate the edge key");
    Ok(())
  }

  #[test]
  fn test_reparent_edge_rejects_making_the_target_its_own_source() -> Result<(), Report> {
    let (mut graph, [_, _, _, c], a_to_c) = fixture()?;
    assert!(graph.reparent_edge(a_to_c, c).is_err());
    Ok(())
  }

  #[test]
  fn test_reparent_edge_rejects_a_duplicate_connection() -> Result<(), Report> {
    let (mut graph, [root, _, _, c], a_to_c) = fixture()?;
    // `root` already reaches `c` directly, so moving `a -> c` under `root` would create a
    // second `root -> c` edge.
    graph.add_edge(root, c, TestEdge(4.0))?;

    assert!(graph.reparent_edge(a_to_c, root).is_err());

    let edge = graph
      .get_edge(a_to_c)
      .expect("a rejected reparent must leave the edge in place");
    assert_eq!(
      edge.read_arc().source(),
      fixture()?.1[1],
      "the edge must not have moved"
    );
    Ok(())
  }

  #[test]
  fn test_reparent_edge_rejects_an_unknown_edge() -> Result<(), Report> {
    let (mut graph, [root, _, _, _], _) = fixture()?;
    assert!(graph.reparent_edge(GraphEdgeKey::invalid(), root).is_err());
    Ok(())
  }
}
