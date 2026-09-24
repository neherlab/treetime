#[cfg(test)]
mod tests {
  use crate::reachability::exists_forward_path_between;
  use eyre::Report;
  use pretty_assertions::assert_eq;

  use self::helpers::graph_with_edges;

  #[test]
  fn test_reachability_follows_directed_chain() -> Result<(), Report> {
    let (graph, keys) = graph_with_edges(3, &[(0, 1), (1, 2)])?;

    let root_to_leaf = exists_forward_path_between(&graph, keys[0], keys[2]);
    let root_to_internal = exists_forward_path_between(&graph, keys[0], keys[1]);
    let internal_to_leaf = exists_forward_path_between(&graph, keys[1], keys[2]);
    let leaf_to_root = exists_forward_path_between(&graph, keys[2], keys[0]);

    assert_eq!(
      (true, true, true, false),
      (root_to_leaf, root_to_internal, internal_to_leaf, leaf_to_root)
    );

    Ok(())
  }

  mod helpers {
    use crate::graph::Graph;
    use crate::node::GraphNodeKey;
    use eyre::Report;

    pub(super) fn graph_with_edges(
      node_count: usize,
      edges: &[(usize, usize)],
    ) -> Result<(Graph, Vec<GraphNodeKey>), Report> {
      let mut graph = Graph::new();
      let keys = std::iter::repeat_with(|| graph.add_node())
        .take(node_count)
        .collect::<Vec<_>>();
      edges
        .iter()
        .try_for_each(|(source, target)| graph.add_edge(keys[*source], keys[*target]).map(|_| ()))?;
      graph.build()?;
      Ok((graph, keys))
    }
  }
}
