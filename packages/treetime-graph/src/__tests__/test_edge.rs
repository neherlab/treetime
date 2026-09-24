#[cfg(test)]
mod tests {
  use super::super::test_graph_support::tests::NamedGraph;
  use crate::edge::invert_edge;
  use eyre::Report;
  use pretty_assertions::assert_eq;

  #[test]
  fn test_edge_invert_swaps_endpoints_and_adjacency() -> Result<(), Report> {
    let mut tree = NamedGraph::new(&["root", "A", "B"], &[("root", "A"), ("root", "B")])?;
    let inverted = tree.edge("root", "A")?;
    let untouched = tree.edge("root", "B")?;

    invert_edge(&mut tree.graph, inverted);

    let edge = tree.graph.get_edge(inverted).expect("edge survives");
    let root = tree.graph.get_node(tree.key("root")).expect("node exists");
    let leaf = tree.graph.get_node(tree.key("A")).expect("node exists");
    assert_eq!(
      (("A", "root"), vec![untouched], vec![inverted], vec![inverted], vec![]),
      (
        (tree.name(edge.source()), tree.name(edge.target())),
        root.outbound().to_vec(),
        root.inbound().to_vec(),
        leaf.outbound().to_vec(),
        leaf.inbound().to_vec(),
      )
    );
    Ok(())
  }
}
