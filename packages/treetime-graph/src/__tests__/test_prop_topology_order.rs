#[cfg(test)]
mod tests {
  use self::helpers::{former_postorder, spec_cases};
  use super::super::test_graph_support::tests::generators::gen_dag;
  use super::super::test_graph_support::tests::helpers::graph_from_dag;
  use crate::key_map::EdgeMap;
  use crate::key_map::NodeMap;
  use crate::node_names::NodeNames;
  use crate::schedule::GraphSchedule;
  use itertools::Itertools;
  use proptest::collection::vec;
  use proptest::prelude::*;

  proptest! {
    #[test]
    fn test_prop_topology_order_on_schedule_postorder_equals_order_on_former_postorder(
      dag in gen_dag(40),
      lengths in vec(0.0..1.0_f64, 120),
      positions in Just((0..40_usize).collect_vec()).prop_shuffle(),
    ) {
      let graph = graph_from_dag(&dag).unwrap();
      let input = graph.node_keys().map(|key| (key, Some(format!("n{key}")))).collect();
      let names = NodeNames::from_input(input, &graph).unwrap().names;
      let branch_lengths: EdgeMap<_> = graph
        .get_edges()
        .map(|edge| (edge.key(), Some(lengths[edge.key().as_usize()])))
        .collect();
      let target_order: NodeMap<_> = graph.node_keys().map(|key| (key, positions[key.as_usize()])).collect();
      let schedule = GraphSchedule::new(&graph).unwrap();
      let former = former_postorder(&graph);

      for spec in spec_cases(&target_order) {
        let expected = spec.sorted_order(&graph, &former, &names, &branch_lengths).unwrap();
        let actual = spec.sorted_order(&graph, schedule.postorder(), &names, &branch_lengths).unwrap();
        prop_assert_eq!(expected, actual, "preset {:?}, aggregate {:?}", spec.preset, spec.target_aggregate);
      }
    }
  }

  mod helpers {
    use crate::graph::Graph;
    use crate::key_map::NodeMap;
    use crate::node::GraphNodeKey;
    use crate::topology_order::{TopologyOrderPreset, TopologyOrderSpec, TopologyOrderTargetAggregate};
    use itertools::Itertools;
    use std::collections::{BTreeMap, VecDeque};

    pub(super) fn spec_cases(target_order: &NodeMap<usize>) -> Vec<TopologyOrderSpec> {
      let presets = [
        TopologyOrderPreset::Keep,
        TopologyOrderPreset::DescendantCount,
        TopologyOrderPreset::DescendantCountReverse,
        TopologyOrderPreset::Height,
        TopologyOrderPreset::HeightReverse,
        TopologyOrderPreset::Divergence,
        TopologyOrderPreset::DivergenceReverse,
        TopologyOrderPreset::Label,
        TopologyOrderPreset::LabelReverse,
        TopologyOrderPreset::TargetOrder,
        TopologyOrderPreset::TargetOrderReverse,
      ];
      let aggregates = [TopologyOrderTargetAggregate::Mean, TopologyOrderTargetAggregate::Median];
      presets
        .into_iter()
        .cartesian_product(aggregates)
        .map(|(preset, target_aggregate)| TopologyOrderSpec {
          preset,
          target_order: target_order.clone(),
          target_aggregate,
        })
        .collect_vec()
    }

    pub(super) fn former_postorder(graph: &Graph) -> Vec<GraphNodeKey> {
      let mut remaining_inbound: BTreeMap<_, _> = graph
        .get_nodes()
        .map(|node| (node.key(), node.inbound().len()))
        .collect();
      let mut queue: VecDeque<_> = remaining_inbound
        .iter()
        .filter_map(|(&key, &count)| (count == 0).then_some(key))
        .collect();
      let mut ordered = vec![];
      while let Some(key) = queue.pop_front() {
        ordered.push(key);
        for (child, _) in graph.children_keys_of(&graph[key]) {
          let count = remaining_inbound.get_mut(&child).unwrap();
          *count -= 1;
          if *count == 0 {
            queue.push_back(child);
          }
        }
      }
      ordered.reverse();
      ordered
    }
  }
}
