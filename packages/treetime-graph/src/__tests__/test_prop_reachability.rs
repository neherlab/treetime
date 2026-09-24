#[cfg(test)]
mod tests {
  use crate::reachability::exists_forward_path_between;
  use proptest::prelude::*;
  use proptest::proptest;

  use self::helpers::graph_chain;

  proptest! {
    #[test]
    fn test_prop_reachability_queries_are_independent(
      (node_count, start_a, finish_a, start_b, finish_b) in (2_usize..32).prop_flat_map(|node_count| {
        (Just(node_count), 0..node_count, 0..node_count, 0..node_count, 0..node_count)
      }),
    ) {
      let (graph, keys) = graph_chain(node_count).unwrap();

      let actual_a_first = exists_forward_path_between(&graph, keys[start_a], keys[finish_a]);
      let actual_b = exists_forward_path_between(&graph, keys[start_b], keys[finish_b]);
      let actual_a_second = exists_forward_path_between(&graph, keys[start_a], keys[finish_a]);

      prop_assert_eq!(start_a <= finish_a, actual_a_first);
      prop_assert_eq!(start_b <= finish_b, actual_b);
      prop_assert_eq!(actual_a_first, actual_a_second);
    }
  }

  mod helpers {
    use crate::graph::Graph;
    use crate::node::GraphNodeKey;
    use eyre::Report;
    use itertools::Itertools;

    pub(super) fn graph_chain(node_count: usize) -> Result<(Graph, Vec<GraphNodeKey>), Report> {
      let mut graph = Graph::new();
      let keys = std::iter::repeat_with(|| graph.add_node())
        .take(node_count)
        .collect_vec();
      keys
        .iter()
        .tuple_windows()
        .try_for_each(|(source, target)| graph.add_edge(*source, *target).map(|_| ()))?;
      graph.build()?;
      Ok((graph, keys))
    }
  }
}
