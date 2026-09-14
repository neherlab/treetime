#[cfg(test)]
mod tests {
  use eyre::Report;
  use pretty_assertions::assert_eq;

  use treetime_graph::edge::{GraphEdgeKey, invert_edge};
  use treetime_io::nwk::nwk_read_str;

  #[test]
  fn test_edge_inverts() -> Result<(), Report> {
    let mut graph = nwk_read_str("((((h:0.7)e:0.6)d:0.4)b:0.,((g:0.5)c:0.2,(i:0.8)f:0.3)a:0.1)r1;")?.graph;

    let edge_key = GraphEdgeKey(3);
    let (input_source, input_target) = {
      let edge = graph.get_edge(edge_key).unwrap();
      (edge.source(), edge.target())
    };

    invert_edge(&mut graph, edge_key);

    let (output_source, output_target) = {
      let edge = graph.get_edge(edge_key).unwrap();
      (edge.source(), edge.target())
    };

    assert_eq!(input_source, output_target);
    assert_eq!(input_target, output_source);

    Ok(())
  }
}
