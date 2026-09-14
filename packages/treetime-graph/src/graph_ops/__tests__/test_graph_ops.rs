#[cfg(test)]
mod tests {
  use crate::edge::GraphEdgeKey;
  use crate::graph::Graph;
  use crate::node::GraphNodeKey;
  use eyre::Report;
  use pretty_assertions::assert_eq;

  type TestGraph = Graph;

  /// `root -> {a, b}`, `a -> c`. Returns the graph and the keys in that order.
  fn fixture() -> Result<(TestGraph, [GraphNodeKey; 4], GraphEdgeKey), Report> {
    let mut graph = TestGraph::new();
    let root = graph.add_node();
    let a = graph.add_node();
    let b = graph.add_node();
    let c = graph.add_node();
    graph.add_edge(root, a)?;
    graph.add_edge(root, b)?;
    let a_to_c = graph.add_edge(a, c)?;
    graph.build()?;
    Ok((graph, [root, a, b, c], a_to_c))
  }

  fn outbound(graph: &TestGraph, node_key: GraphNodeKey) -> Vec<GraphEdgeKey> {
    graph.get_node(node_key).expect("node exists").outbound().to_vec()
  }

  fn inbound(graph: &TestGraph, node_key: GraphNodeKey) -> Vec<GraphEdgeKey> {
    graph.get_node(node_key).expect("node exists").inbound().to_vec()
  }

  #[test]
  fn test_reparent_edge_moves_the_edge_and_keeps_key() -> Result<(), Report> {
    let (mut graph, [root, a, b, c], a_to_c) = fixture()?;

    graph.reparent_edge(a_to_c, b)?;

    let edge = graph.get_edge(a_to_c).expect("edge survives reparenting");
    assert_eq!(edge.key(), a_to_c, "the edge keeps its key (not rebuilt)");
    assert_eq!(edge.source(), b);
    assert_eq!(edge.target(), c);

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
    graph.add_edge(root, c)?;

    assert!(graph.reparent_edge(a_to_c, root).is_err());

    let edge = graph
      .get_edge(a_to_c)
      .expect("a rejected reparent must leave the edge in place");
    assert_eq!(edge.source(), fixture()?.1[1], "the edge must not have moved");
    Ok(())
  }

  #[test]
  fn test_reparent_edge_rejects_an_unknown_edge() -> Result<(), Report> {
    let (mut graph, [root, _, _, _], _) = fixture()?;
    assert!(graph.reparent_edge(GraphEdgeKey::invalid(), root).is_err());
    Ok(())
  }
}
