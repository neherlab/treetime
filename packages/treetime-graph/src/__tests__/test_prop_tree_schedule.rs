#[cfg(test)]
mod tests {
  use super::super::test_graph_support::tests::helpers::graph_from_parents;
  use crate::node::Node;
  use crate::schedule::GraphSchedule;
  use itertools::Itertools;
  use proptest::prelude::*;

  proptest! {
    #[test]
    fn test_prop_tree_schedule_follows_graph_edges_with_parents_first(
      parents in (1_usize..40).prop_flat_map(|count| {
        (0..count).map(|index| 0..index.max(1)).collect_vec()
      }),
    ) {
      let graph = graph_from_parents(&parents).unwrap();

      let tree = GraphSchedule::new(&graph).unwrap().into_tree().unwrap();

      prop_assert_eq!(vec![tree.root()], graph.get_roots().map(Node::key).collect_vec());
      prop_assert_eq!(tree.schedule().preorder(), tree.preorder());
      let position = |target| tree.preorder().iter().position(|&key| key == target);
      for node in graph.get_nodes() {
        let key = node.key();
        prop_assert_eq!(graph.node_parent(key).unwrap(), tree.parent(key));
        prop_assert_eq!(graph.children_keys_of(node).collect_vec(), tree.children(key).to_vec());
        if let Some((parent, _)) = tree.parent(key) {
          prop_assert!(position(parent) < position(key));
        }
      }
    }
  }
}
