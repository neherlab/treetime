#[cfg(test)]
mod tests {
  use self::helpers::{former_postorder, former_preorder, positions};
  use super::super::test_graph_support::tests::generators::gen_dag;
  use super::super::test_graph_support::tests::helpers::{graph_from_dag, graph_from_parents};
  use crate::node::Node;
  use crate::schedule::GraphSchedule;
  use itertools::Itertools;
  use proptest::prelude::*;

  proptest! {
    #[test]
    fn test_prop_schedule_orders_list_each_node_once_after_its_parents_or_children(dag in gen_dag(40)) {
      let graph = graph_from_dag(&dag).unwrap();

      let schedule = GraphSchedule::new(&graph).unwrap();

      let keys = graph.get_nodes().map(Node::key).collect_vec();
      prop_assert_eq!(&keys, &schedule.preorder().iter().copied().sorted().collect_vec());
      prop_assert_eq!(&keys, &schedule.postorder().iter().copied().sorted().collect_vec());
      let (pre, post) = (positions(schedule.preorder()), positions(schedule.postorder()));
      for &key in &keys {
        for &(parent, _) in schedule.parents(key) {
          prop_assert!(pre[&parent] < pre[&key], "node {key} precedes its parent {parent} in the preorder");
        }
        for &(child, _) in schedule.children(key) {
          prop_assert!(post[&child] < post[&key], "node {key} precedes its child {child} in the postorder");
        }
      }
    }

    #[test]
    fn test_prop_schedule_adjacency_roots_leaves_and_slots_follow_graph_edges(dag in gen_dag(40)) {
      let graph = graph_from_dag(&dag).unwrap();

      let schedule = GraphSchedule::new(&graph).unwrap();

      for node in graph.get_nodes() {
        let key = node.key();
        prop_assert_eq!(graph.parents_keys_of(node).collect_vec(), schedule.parents(key).to_vec());
        prop_assert_eq!(graph.children_keys_of(node).collect_vec(), schedule.children(key).to_vec());
        for (slot, &(_, edge)) in schedule.parents(key).iter().enumerate() {
          prop_assert_eq!(slot, schedule.parent_slot(edge));
        }
        for (slot, &(_, edge)) in schedule.children(key).iter().enumerate() {
          prop_assert_eq!(slot, schedule.child_slot(edge));
        }
      }
      prop_assert_eq!(graph.get_roots().map(Node::key).collect_vec(), schedule.roots().to_vec());
      prop_assert_eq!(
        graph.get_leaves().map(Node::key).sorted().collect_vec(),
        schedule.leaves().iter().copied().sorted().collect_vec()
      );
    }

    #[test]
    fn test_prop_schedule_tree_orders_match_graph_traversals(
      parents in (1_usize..40).prop_flat_map(|count| {
        (0..count).map(|index| 0..index.max(1)).collect_vec()
      }),
    ) {
      let graph = graph_from_parents(&parents).unwrap();

      let schedule = GraphSchedule::new(&graph).unwrap();

      prop_assert_eq!(former_preorder(&graph), schedule.preorder().to_vec());
      prop_assert_eq!(former_postorder(&graph), schedule.postorder().to_vec());
      schedule.require_tree().unwrap();
    }
  }

  mod helpers {
    use crate::graph::Graph;
    use crate::node::GraphNodeKey;
    use std::collections::{BTreeMap, BTreeSet};

    pub(super) fn former_preorder(graph: &Graph) -> Vec<GraphNodeKey> {
      let mut order = vec![];
      let mut stack = vec![single_root(graph)];
      while let Some(key) = stack.pop() {
        order.push(key);
        stack.extend(graph.children_keys_of(&graph[key]).rev().map(|(child, _)| child));
      }
      order
    }

    pub(super) fn former_postorder(graph: &Graph) -> Vec<GraphNodeKey> {
      let mut order = vec![];
      let mut stack = vec![single_root(graph)];
      let mut visited = BTreeSet::new();
      while let Some(key) = stack.pop() {
        if visited.insert(key) {
          stack.push(key);
          let unvisited = graph
            .children_keys_of(&graph[key])
            .rev()
            .map(|(child, _)| child)
            .filter(|child| !visited.contains(child));
          stack.extend(unvisited);
        } else {
          order.push(key);
        }
      }
      order
    }

    fn single_root(graph: &Graph) -> GraphNodeKey {
      let &[root] = graph.roots.as_slice() else {
        panic!("the former traversals require exactly one root");
      };
      root
    }

    pub(super) fn positions(order: &[GraphNodeKey]) -> BTreeMap<GraphNodeKey, usize> {
      order
        .iter()
        .enumerate()
        .map(|(position, &key)| (key, position))
        .collect()
    }
  }
}
