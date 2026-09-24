#[cfg(test)]
mod tests {
  use proptest::prelude::*;

  use self::generators::gen_tree_with_edge;
  use self::helpers::{adjacency_from_edges, adjacency_from_nodes, build_tree, children_by_node};

  proptest! {
    #[test]
    fn test_prop_graph_ops_collapse_edge_preserves_adjacency((parents, collapsed) in gen_tree_with_edge()) {
      let mut graph = build_tree(&parents).unwrap();
      let collapsed_key = graph.get_edges().nth(collapsed).unwrap().key();

      graph.collapse_edge(collapsed_key).unwrap();

      prop_assert_eq!(adjacency_from_edges(&graph), adjacency_from_nodes(&graph));
      prop_assert_eq!((parents.len(), parents.len() - 1), (graph.get_nodes().count(), graph.get_edges().count()));
    }

    #[test]
    fn test_prop_graph_ops_collapse_edge_moves_target_children_to_source((parents, collapsed) in gen_tree_with_edge()) {
      let mut graph = build_tree(&parents).unwrap();
      let edge = graph.get_edges().nth(collapsed).unwrap().clone();
      let mut expected = children_by_node(&graph);
      let target_children = expected.remove(&edge.target()).unwrap();
      let source_children = expected.get_mut(&edge.source()).unwrap();
      source_children.retain(|child| *child != edge.target());
      source_children.extend(target_children);

      graph.collapse_edge(edge.key()).unwrap();

      prop_assert_eq!(expected, children_by_node(&graph));
    }
  }

  mod generators {
    use proptest::prelude::*;

    pub(super) fn gen_tree_with_edge() -> impl Strategy<Value = (Vec<usize>, usize)> {
      (2_usize..24).prop_flat_map(|node_count| {
        let parents = (1..node_count).map(|child| 0..child).collect::<Vec<_>>();
        (parents, 0..node_count - 1)
      })
    }
  }

  mod helpers {
    use crate::edge::GraphEdgeKey;
    use crate::graph::Graph;
    use crate::node::GraphNodeKey;
    use eyre::Report;
    use itertools::Itertools;
    use std::collections::BTreeMap;

    type Adjacency = BTreeMap<GraphNodeKey, (Vec<GraphEdgeKey>, Vec<GraphEdgeKey>)>;

    pub(super) fn build_tree(parents: &[usize]) -> Result<Graph, Report> {
      let mut graph = Graph::new();
      let keys = std::iter::repeat_with(|| graph.add_node())
        .take(parents.len() + 1)
        .collect_vec();
      for (child, parent) in parents.iter().enumerate() {
        graph.add_edge(keys[*parent], keys[child + 1])?;
      }
      graph.build()?;
      Ok(graph)
    }

    pub(super) fn adjacency_from_edges(graph: &Graph) -> Adjacency {
      let mut adjacency: Adjacency = graph.get_nodes().map(|node| (node.key(), (vec![], vec![]))).collect();
      for edge in graph.get_edges() {
        adjacency.entry(edge.source()).or_default().0.push(edge.key());
        adjacency.entry(edge.target()).or_default().1.push(edge.key());
      }
      adjacency
    }

    pub(super) fn adjacency_from_nodes(graph: &Graph) -> Adjacency {
      graph
        .get_nodes()
        .map(|node| {
          let outbound = node.outbound().iter().copied().sorted().collect_vec();
          let inbound = node.inbound().iter().copied().sorted().collect_vec();
          (node.key(), (outbound, inbound))
        })
        .collect()
    }

    pub(super) fn children_by_node(graph: &Graph) -> BTreeMap<GraphNodeKey, Vec<GraphNodeKey>> {
      graph
        .get_nodes()
        .map(|node| {
          let children = graph.children_keys_of(node).map(|(child, _)| child).collect_vec();
          (node.key(), children)
        })
        .collect()
    }
  }
}
