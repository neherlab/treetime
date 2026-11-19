#[cfg(test)]
mod tests {
  use self::helpers::{fixture, inbound, outbound};
  use crate::edge::GraphEdgeKey;
  use eyre::Report;
  use pretty_assertions::assert_eq;
  use treetime_utils::assert_error;

  #[test]
  fn test_reparent_edge_moves_the_edge_and_keeps_key() -> Result<(), Report> {
    let (mut graph, [root, a, b, c], a_to_c) = fixture()?;
    let root_outbound = outbound(&graph, root);

    graph.reparent_edge(a_to_c, b)?;

    let edge = graph.get_edge(a_to_c).expect("the edge keeps its key");
    assert_eq!(
      ((b, c), vec![], vec![a_to_c], vec![a_to_c], root_outbound),
      (
        (edge.source(), edge.target()),
        outbound(&graph, a),
        outbound(&graph, b),
        inbound(&graph, c),
        outbound(&graph, root),
      )
    );
    Ok(())
  }

  #[test]
  fn test_reparent_edge_to_the_current_source_is_a_noop() -> Result<(), Report> {
    let (mut graph, [_, a, _, _], a_to_c) = fixture()?;
    let before = outbound(&graph, a);

    graph.reparent_edge(a_to_c, a)?;

    assert_eq!(before, outbound(&graph, a), "a no-op must not duplicate the edge key");
    Ok(())
  }

  #[test]
  fn test_reparent_edge_rejects_making_the_target_its_own_source() -> Result<(), Report> {
    let (mut graph, [_, _, _, c], a_to_c) = fixture()?;
    assert_error!(
      graph.reparent_edge(a_to_c, c),
      format!("When reparenting edge {a_to_c} to {c}: Attempted to connect node {c} to itself.")
    );
    Ok(())
  }

  #[test]
  fn test_reparent_edge_rejects_a_duplicate_connection() -> Result<(), Report> {
    let (mut graph, [root, a, _, c], a_to_c) = fixture()?;
    graph.add_edge(root, c)?;

    assert_error!(
      graph.reparent_edge(a_to_c, root),
      format!("When reparenting edge {a_to_c} to {root}: Nodes {root} and {c} are already connected.")
    );

    let edge = graph
      .get_edge(a_to_c)
      .expect("a rejected reparent must leave the edge in place");
    assert_eq!(a, edge.source(), "the edge must not have moved");
    Ok(())
  }

  #[test]
  fn test_reparent_edge_rejects_an_unknown_edge() -> Result<(), Report> {
    let (mut graph, [root, _, _, _], _) = fixture()?;
    let unknown = GraphEdgeKey(usize::MAX);
    assert_error!(
      graph.reparent_edge(unknown, root),
      format!(
        "When reparenting edge {unknown}: edge not found. This is an internal error. Please report it to developers."
      )
    );
    Ok(())
  }

  mod helpers {
    use crate::edge::GraphEdgeKey;
    use crate::graph::Graph;
    use crate::node::GraphNodeKey;
    use eyre::Report;

    pub(super) fn fixture() -> Result<(Graph, [GraphNodeKey; 4], GraphEdgeKey), Report> {
      let mut graph = Graph::new();
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

    pub(super) fn outbound(graph: &Graph, node_key: GraphNodeKey) -> Vec<GraphEdgeKey> {
      graph.get_node(node_key).expect("node exists").outbound().to_vec()
    }

    pub(super) fn inbound(graph: &Graph, node_key: GraphNodeKey) -> Vec<GraphEdgeKey> {
      graph.get_node(node_key).expect("node exists").inbound().to_vec()
    }
  }
}
