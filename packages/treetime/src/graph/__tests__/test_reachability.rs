#[cfg(test)]
mod tests {
  use eyre::Report;
  use pretty_assertions::assert_eq;
  use treetime_graph::reachability::exists_forward_path_between;

  use self::helpers::graph_with_edges;

  #[test]
  fn test_reachability_follows_directed_chain() -> Result<(), Report> {
    let (graph, nodes) = graph_with_edges(3, &[(0, 1), (1, 2)])?;

    let root_to_leaf = exists_forward_path_between(&graph, &nodes[0], &nodes[2]);
    let root_to_internal = exists_forward_path_between(&graph, &nodes[0], &nodes[1]);
    let internal_to_leaf = exists_forward_path_between(&graph, &nodes[1], &nodes[2]);
    let leaf_to_root = exists_forward_path_between(&graph, &nodes[2], &nodes[0]);

    // Forward reachability follows edge directions; the reverse direction is not reachable.
    assert_eq!(
      (true, true, true, false),
      (root_to_leaf, root_to_internal, internal_to_leaf, leaf_to_root)
    );

    Ok(())
  }

  mod helpers {
    use eyre::Report;
    use treetime_graph::graph::{Graph, SafeNode};

    use crate::test_utils::{TestEdge, TestNode};

    pub fn graph_with_edges(
      node_count: usize,
      edges: &[(usize, usize)],
    ) -> Result<(Graph<TestNode, TestEdge, ()>, Vec<SafeNode<TestNode>>), Report> {
      let mut graph = Graph::new();
      let keys = (0..node_count)
        .map(|index| graph.add_node(TestNode(Some(format!("node_{index}")))))
        .collect::<Vec<_>>();
      edges
        .iter()
        .try_for_each(|(source, target)| graph.add_edge(keys[*source], keys[*target], TestEdge(None)).map(|_| ()))?;
      graph.build()?;
      let nodes = keys.iter().map(|key| graph.get_node(*key).unwrap()).collect::<Vec<_>>();
      Ok((graph, nodes))
    }
  }
}
