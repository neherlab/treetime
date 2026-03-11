use crate::gtr::gtr::GTR;
use crate::representation::discrete_states::DiscreteStates;
use crate::representation::partition::traits::HasLogLh;
use crate::representation::payload::discrete::{DiscreteEdgeData, DiscreteNodeData};
use ndarray::Array1;
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::node::GraphNodeKey;

#[derive(Clone, Debug)]
pub struct PartitionDiscrete {
  pub index: usize,
  pub gtr: GTR,
  pub states: DiscreteStates,
  pub nodes: BTreeMap<GraphNodeKey, DiscreteNodeData>,
  pub edges: BTreeMap<GraphEdgeKey, DiscreteEdgeData>,
}

impl PartitionDiscrete {
  pub fn new(index: usize, gtr: GTR, states: DiscreteStates) -> Self {
    Self {
      index,
      gtr,
      states,
      nodes: BTreeMap::new(),
      edges: BTreeMap::new(),
    }
  }

  pub fn n_states(&self) -> usize {
    self.states.len()
  }

  pub fn get_reconstructed_trait(&self, node_key: GraphNodeKey) -> Option<String> {
    let node = self.nodes.get(&node_key)?;
    let argmax = argmax_first_1d(&node.profile)?;
    Some(self.states.get_name(argmax).to_owned())
  }

  pub fn get_confidence(&self, node_key: GraphNodeKey) -> Option<&Array1<f64>> {
    self.nodes.get(&node_key).map(|n| &n.profile)
  }
}

impl HasLogLh for PartitionDiscrete {
  fn get_log_lh(&self, node_key: GraphNodeKey) -> f64 {
    self.nodes.get(&node_key).map_or(0.0, |node| node.log_lh)
  }
}

fn argmax_first_1d(arr: &Array1<f64>) -> Option<usize> {
  if arr.is_empty() {
    return None;
  }
  let mut max_idx = 0;
  let mut max_val = arr[0];
  for (i, &v) in arr.iter().enumerate().skip(1) {
    if v > max_val {
      max_val = v;
      max_idx = i;
    }
  }
  Some(max_idx)
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::gtr::gtr::GTRParams;
  use approx::assert_abs_diff_eq;
  use ndarray::array;
  use pretty_assertions::assert_eq;

  fn make_test_partition() -> PartitionDiscrete {
    let states = DiscreteStates::from_values(["USA", "Germany", "China"].iter().copied(), "?");
    let gtr = GTR::new(GTRParams {
      n_states: 3,
      mu: 1.0,
      W: None,
      pi: array![1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0],
    })
    .unwrap();
    PartitionDiscrete::new(0, gtr, states)
  }

  #[test]
  fn test_new_partition() {
    let partition = make_test_partition();

    assert_eq!(partition.index, 0);
    assert_eq!(partition.n_states(), 3);
    assert!(partition.nodes.is_empty());
    assert!(partition.edges.is_empty());
  }

  #[test]
  fn test_get_reconstructed_trait() {
    let mut partition = make_test_partition();
    let node_key = GraphNodeKey(0);

    partition.nodes.insert(
      node_key,
      DiscreteNodeData {
        observed: Some(1),
        profile: array![0.1, 0.7, 0.2],
        log_lh: 0.0,
      },
    );

    let trait_name = partition.get_reconstructed_trait(node_key);
    assert_eq!(trait_name, Some("Germany".to_owned()));
  }

  #[test]
  fn test_get_confidence() {
    let mut partition = make_test_partition();
    let node_key = GraphNodeKey(0);
    partition.nodes.insert(
      node_key,
      DiscreteNodeData {
        observed: Some(1),
        profile: array![0.1, 0.7, 0.2],
        log_lh: 0.0,
      },
    );

    let confidence = partition.get_confidence(node_key).unwrap();
    assert_abs_diff_eq!(confidence[0], 0.1, epsilon = 1e-10);
    assert_abs_diff_eq!(confidence[1], 0.7, epsilon = 1e-10);
    assert_abs_diff_eq!(confidence[2], 0.2, epsilon = 1e-10);
  }

  #[test]
  fn test_get_log_lh() {
    let mut partition = make_test_partition();
    let node_key = GraphNodeKey(0);

    partition.nodes.insert(
      node_key,
      DiscreteNodeData {
        observed: None,
        profile: array![0.25, 0.25, 0.25],
        log_lh: -2.5,
      },
    );

    assert_abs_diff_eq!(partition.get_log_lh(node_key), -2.5, epsilon = 1e-10);
  }

  #[test]
  fn test_argmax_first_1d() {
    assert_eq!(argmax_first_1d(&array![0.1, 0.5, 0.4]), Some(1));
    assert_eq!(argmax_first_1d(&array![0.5, 0.5, 0.0]), Some(0)); // tie: first wins
    assert_eq!(argmax_first_1d(&Array1::zeros(0)), None);
  }
}
