use crate::gtr::gtr::GTR;
use crate::representation::discrete_states::DiscreteStates;
use crate::representation::partition::traits::HasLogLh;
use crate::representation::payload::discrete::{DiscreteEdgeData, DiscreteNodeData};
use eyre::Report;
use ndarray::Array1;
use std::collections::BTreeMap;
use treetime_graph::edge::{EdgeOptimizeOps, GraphEdgeKey};
use treetime_graph::graph::Graph;
use treetime_graph::graph_traverse::{GraphNodeBackward, GraphNodeForward};
use treetime_graph::node::{GraphNode, GraphNodeKey};

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

  pub fn process_node_backward<N, E>(&mut self, node: &GraphNodeBackward<N, E, ()>) -> Result<(), Report>
  where
    N: GraphNode,
    E: EdgeOptimizeOps,
  {
    let n_states = self.n_states();

    let msg_to_parent = if node.is_leaf {
      self.nodes[&node.key].profile.clone()
    } else {
      // Combine child messages via element-wise product in log space
      let child_msgs: Vec<Array1<f64>> = node
        .child_keys
        .iter()
        .map(|(_, edge_key)| self.edges[edge_key].msg_from_child.clone())
        .collect();

      let mut log_profile = child_msgs[0].mapv(f64::ln);
      for msg in &child_msgs[1..] {
        log_profile += &msg.mapv(f64::ln);
      }

      // Normalize from log space
      let (profile, delta_ll) = normalize_from_log_1d(&log_profile);

      let log_lh: f64 = node
        .child_keys
        .iter()
        .map(|(_, edge_key)| self.edges[edge_key].log_lh)
        .sum();

      // Store internal node data
      self.nodes.insert(
        node.key,
        DiscreteNodeData {
          observed: None,
          profile: profile.clone(),
          log_lh: log_lh + delta_ll,
        },
      );

      profile
    };

    if node.is_root {
      // Root: weight by equilibrium frequencies
      let mut profile = &msg_to_parent * &self.gtr.pi;
      let delta_ll = normalize_inplace_1d(&mut profile);

      let node_data = self.nodes.get_mut(&node.key).unwrap();
      node_data.profile = profile;
      node_data.log_lh += delta_ll;
    } else {
      // Non-root: propagate message to parent via transition matrix
      let edge_key = node.parent_edge_keys[0];
      let branch_length = node.parent_edges[0].branch_length().unwrap_or(0.0);
      let exp_qt = self.gtr.expQt(branch_length);

      // msg_from_child[j] = sum_i msg_to_parent[i] * P(i->j|t)
      // P is column-stochastic: P[i,j] = P(j->i), so we need P^T
      let msg_from_child = exp_qt.t().dot(&msg_to_parent);

      let log_lh = self.nodes.get(&node.key).map_or(0.0, |n| n.log_lh);

      self.edges.insert(
        edge_key,
        DiscreteEdgeData {
          msg_to_parent,
          msg_from_child,
          msg_to_child: Array1::zeros(n_states),
          log_lh,
        },
      );
    }

    Ok(())
  }

  pub fn process_node_forward<N, E>(
    &mut self,
    graph: &Graph<N, E, ()>,
    node: &GraphNodeForward<N, E, ()>,
  ) -> Result<(), Report>
  where
    N: GraphNode,
    E: EdgeOptimizeOps,
  {
    if !node.is_root {
      // Combine messages from parent
      let mut msgs_to_combine: Vec<Array1<f64>> = vec![];
      let mut log_lh = 0.0;

      for (_, edge_key) in &node.parent_keys {
        let edge = &self.edges[edge_key];
        let edge_payload = graph.get_edge(*edge_key).unwrap().read_arc().payload().read_arc();
        let branch_length = edge_payload.branch_length().unwrap_or(0.0);
        let exp_qt = self.gtr.expQt(branch_length);

        msgs_to_combine.push(edge.msg_to_parent.clone());
        log_lh += edge.log_lh;

        // msg_to_child propagated through transition matrix
        let propagated = exp_qt.dot(&edge.msg_to_child);
        msgs_to_combine.push(propagated);
      }

      // Element-wise product and normalize
      let mut profile = msgs_to_combine[0].clone();
      for msg in &msgs_to_combine[1..] {
        profile *= msg;
      }
      let delta_ll = normalize_inplace_1d(&mut profile);
      log_lh += delta_ll;

      let node_data = self.nodes.get_mut(&node.key).unwrap();
      node_data.profile = profile;
      node_data.log_lh = log_lh;
    }

    // Compute msg_to_child for each child edge
    for child_edge_key in &node.child_edge_keys {
      let node_data = &self.nodes[&node.key];
      let edge_data = self.edges.get_mut(child_edge_key).unwrap();

      // msg_to_child = node_profile / msg_from_child (element-wise)
      let mut msg_to_child = &node_data.profile / &edge_data.msg_from_child;
      normalize_inplace_1d(&mut msg_to_child);
      edge_data.msg_to_child = msg_to_child;
    }

    Ok(())
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

fn normalize_inplace_1d(arr: &mut Array1<f64>) -> f64 {
  let sum = arr.sum();
  if sum > 0.0 {
    *arr /= sum;
  }
  sum.ln()
}

fn normalize_from_log_1d(log_arr: &Array1<f64>) -> (Array1<f64>, f64) {
  let max_val = log_arr.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b));

  let mut arr = log_arr.mapv(|v| (v - max_val).exp());
  let sum = arr.sum();
  if sum > 0.0 {
    arr /= sum;
  }

  let log_lh = max_val + sum.ln();
  (arr, log_lh)
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
