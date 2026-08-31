#[cfg(test)]
mod tests {
  use proptest::prelude::*;
  use proptest::proptest;
  use treetime_graph::reachability::exists_forward_path_between;

  use self::helpers::graph_chain;

  proptest! {
    #[test]
    fn test_prop_graph_path_queries_are_independent(
      (node_count, start_a, finish_a, start_b, finish_b) in (2_usize..32).prop_flat_map(|node_count| {
        (Just(node_count), 0..node_count, 0..node_count, 0..node_count, 0..node_count)
      }),
    ) {
      let (graph, nodes) = graph_chain(node_count).unwrap();

      let actual_a_first = exists_forward_path_between(&graph, &nodes[start_a], &nodes[finish_a]);
      let actual_b = exists_forward_path_between(&graph, &nodes[start_b], &nodes[finish_b]);
      let actual_a_second = exists_forward_path_between(&graph, &nodes[start_a], &nodes[finish_a]);

      // A directed chain reaches exactly the nodes at or after the starting index.
      prop_assert_eq!(start_a <= finish_a, actual_a_first);
      prop_assert_eq!(start_b <= finish_b, actual_b);
      prop_assert_eq!(actual_a_first, actual_a_second);
    }
  }

  mod helpers {
    use eyre::Report;
    use itertools::Itertools;
    use treetime_graph::graph::{Graph, SafeNode};

    use crate::test_utils::{TestEdge, TestNode};

    pub fn graph_chain(node_count: usize) -> Result<(Graph<TestNode, TestEdge, ()>, Vec<SafeNode<TestNode>>), Report> {
      let mut graph = Graph::new();
      let keys = (0..node_count)
        .map(|index| graph.add_node(TestNode(Some(format!("node_{index}")))))
        .collect_vec();
      keys
        .iter()
        .tuple_windows()
        .try_for_each(|(source, target)| graph.add_edge(*source, *target, TestEdge(None)).map(|_| ()))?;
      graph.build()?;
      let nodes = keys.iter().map(|key| graph.get_node(*key).unwrap()).collect_vec();
      Ok((graph, nodes))
    }
  }
}
