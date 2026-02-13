#[cfg(test)]
mod tests {
  use eyre::Report;
  use pretty_assertions::assert_eq;

  use crate::graph::graph_tests::tests::{TestEdge, TestNode};
  use treetime_graph::edge::{GraphEdgeKey, invert_edge};
  use treetime_io::nwk::nwk_read_str;

  #[test]
  fn edge_inverts() -> Result<(), Report> {
    let mut graph =
      nwk_read_str::<TestNode, TestEdge, ()>("((((h:0.7)e:0.6)d:0.4)b:0.,((g:0.5)c:0.2,(i:0.8)f:0.3)a:0.1)r1;")?;

    let edge = graph.get_edge(GraphEdgeKey(3)).unwrap();
    let input_source = edge.read().source();
    let input_target = edge.read().target();

    invert_edge(&mut graph, &edge);

    let edge = graph.get_edge(GraphEdgeKey(3)).unwrap();
    let output_source = edge.read().source();
    let output_target = edge.read().target();

    assert_eq!(input_source, output_target);
    assert_eq!(input_target, output_source);

    Ok(())
  }
}
