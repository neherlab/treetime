#[cfg(test)]
mod __tests__;

use crate::edge::{GraphEdge, GraphEdgeKey, HasBranchLength};
use crate::graph::{Graph, SafeNode};
use crate::node::{GraphNode, GraphNodeKey, Named};
use eyre::Report;
use itertools::Itertools;
use ordered_float::OrderedFloat;
use serde::{Deserialize, Serialize};
use std::cmp::Ordering;
use std::collections::{BTreeMap, VecDeque};
use treetime_utils::{make_error, make_report};

#[derive(Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub struct TopologyOrderSpec {
  pub preset: TopologyOrderPreset,
  pub target_order: Vec<String>,
  pub target_aggregate: TopologyOrderTargetAggregate,
}

impl Default for TopologyOrderSpec {
  fn default() -> Self {
    Self {
      preset: TopologyOrderPreset::DescendantCount,
      target_order: vec![],
      target_aggregate: TopologyOrderTargetAggregate::Mean,
    }
  }
}

impl TopologyOrderSpec {
  pub fn keep() -> Self {
    Self {
      preset: TopologyOrderPreset::Keep,
      ..Self::default()
    }
  }

  pub fn descendant_count(reverse: bool) -> Self {
    Self {
      preset: if reverse {
        TopologyOrderPreset::DescendantCountReverse
      } else {
        TopologyOrderPreset::DescendantCount
      },
      ..Self::default()
    }
  }

  /// Apply the requested logical topology order to a graph.
  ///
  /// All fallible computation and validation completes before the graph is
  /// mutated. Node and edge keys and their slot storage remain unchanged.
  pub fn apply<N, E, D>(&self, graph: &mut Graph<N, E, D>) -> Result<(), Report>
  where
    N: GraphNode + Named,
    E: GraphEdge + HasBranchLength,
    D: Sync + Send,
  {
    let order = if self.preset == TopologyOrderPreset::Keep {
      build_order_unmodified(graph)
    } else {
      let postorder = postorder_keys(graph)?;
      let reverse = self.preset.is_reverse();
      match self.preset {
        TopologyOrderPreset::Keep => unreachable!(),
        TopologyOrderPreset::DescendantCount | TopologyOrderPreset::DescendantCountReverse => {
          let keys = compute_descendant_counts(graph, &postorder);
          build_order(graph, &keys, reverse)
        },
        TopologyOrderPreset::Height | TopologyOrderPreset::HeightReverse => {
          let keys = compute_heights(graph, &postorder);
          build_order(graph, &keys, reverse)
        },
        TopologyOrderPreset::Divergence | TopologyOrderPreset::DivergenceReverse => {
          let keys = compute_divergences(graph, &postorder);
          build_order(graph, &keys, reverse)
        },
        TopologyOrderPreset::Label | TopologyOrderPreset::LabelReverse => {
          let keys = compute_labels(graph, &postorder)?;
          build_order(graph, &keys, reverse)
        },
        TopologyOrderPreset::TargetOrder | TopologyOrderPreset::TargetOrderReverse => {
          let keys = compute_target_scores(graph, &postorder, &self.target_order, self.target_aggregate)?;
          build_order(graph, &keys, reverse)
        },
      }
    }?;

    let ordered_nodes: Vec<(SafeNode<N>, Vec<GraphEdgeKey>)> = order
      .outbound_edges
      .into_iter()
      .map(|(node_key, outbound_edges)| {
        graph
          .get_node(node_key)
          .map(|node| (node, outbound_edges))
          .ok_or_else(|| make_report!("Node {node_key} disappeared while applying topology order"))
      })
      .try_collect()?;

    for (node, outbound_edges) in ordered_nodes {
      *node.write_arc().outbound_mut() = outbound_edges;
    }
    graph.roots = order.roots;
    graph.leaves = order.leaves;
    Ok(())
  }
}

#[derive(Copy, Clone, Debug, Default, Eq, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum TopologyOrderPreset {
  Keep,
  #[default]
  DescendantCount,
  DescendantCountReverse,
  Height,
  HeightReverse,
  Divergence,
  DivergenceReverse,
  Label,
  LabelReverse,
  TargetOrder,
  TargetOrderReverse,
}

impl TopologyOrderPreset {
  pub fn is_target_order(self) -> bool {
    matches!(self, Self::TargetOrder | Self::TargetOrderReverse)
  }

  fn is_reverse(self) -> bool {
    matches!(
      self,
      Self::DescendantCountReverse
        | Self::HeightReverse
        | Self::DivergenceReverse
        | Self::LabelReverse
        | Self::TargetOrderReverse
    )
  }
}

#[derive(Copy, Clone, Debug, Default, Eq, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum TopologyOrderTargetAggregate {
  #[default]
  Mean,
  Median,
}

#[derive(Debug, Eq, PartialEq)]
struct TopologyOrder {
  roots: Vec<GraphNodeKey>,
  leaves: Vec<GraphNodeKey>,
  outbound_edges: BTreeMap<GraphNodeKey, Vec<GraphEdgeKey>>,
}

fn build_order_unmodified<N, E, D>(graph: &Graph<N, E, D>) -> Result<TopologyOrder, Report>
where
  N: GraphNode,
  E: GraphEdge,
  D: Sync + Send,
{
  Ok(TopologyOrder {
    roots: graph.roots.clone(),
    leaves: graph.leaves.clone(),
    outbound_edges: graph
      .nodes
      .iter()
      .filter_map(Option::as_ref)
      .map(|node| {
        let node = node.read_arc();
        (node.key(), node.outbound().to_vec())
      })
      .collect(),
  })
}

fn build_order<N, E, D, K: Ord>(
  graph: &Graph<N, E, D>,
  keys: &BTreeMap<GraphNodeKey, K>,
  reverse: bool,
) -> Result<TopologyOrder, Report>
where
  N: GraphNode,
  E: GraphEdge,
  D: Sync + Send,
{
  let compare = |a: &GraphNodeKey, b: &GraphNodeKey| -> Ordering {
    let ord = keys[a].cmp(&keys[b]);
    if reverse { ord.reverse() } else { ord }
  };

  let sort_node_keys = |node_keys: &[GraphNodeKey]| -> Vec<GraphNodeKey> {
    node_keys
      .iter()
      .copied()
      .enumerate()
      .sorted_by(|(pos_a, a), (pos_b, b)| compare(a, b).then(pos_a.cmp(pos_b)))
      .map(|(_, key)| key)
      .collect_vec()
  };

  let sort_edge_keys = |edge_keys: &[GraphEdgeKey]| -> Vec<GraphEdgeKey> {
    edge_keys
      .iter()
      .copied()
      .enumerate()
      .sorted_by(|(pos_a, a), (pos_b, b)| {
        let target_a = graph.get_target_node_key(*a);
        let target_b = graph.get_target_node_key(*b);
        match (target_a, target_b) {
          (Ok(ta), Ok(tb)) => compare(&ta, &tb).then(pos_a.cmp(pos_b)),
          _ => pos_a.cmp(pos_b),
        }
      })
      .map(|(_, key)| key)
      .collect_vec()
  };

  Ok(TopologyOrder {
    roots: sort_node_keys(&graph.roots),
    leaves: sort_node_keys(&graph.leaves),
    outbound_edges: graph
      .nodes
      .iter()
      .filter_map(Option::as_ref)
      .map(|node| {
        let node = node.read_arc();
        (node.key(), sort_edge_keys(node.outbound()))
      })
      .collect(),
  })
}

#[derive(Clone, Debug, Eq, PartialEq)]
struct TargetScore {
  numerator: usize,
  denominator: usize,
}

impl Ord for TargetScore {
  fn cmp(&self, other: &Self) -> Ordering {
    ((self.numerator as u128) * (other.denominator as u128))
      .cmp(&((other.numerator as u128) * (self.denominator as u128)))
  }
}

impl PartialOrd for TargetScore {
  fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
    Some(self.cmp(other))
  }
}

fn compute_descendant_counts<N, E, D>(
  graph: &Graph<N, E, D>,
  postorder: &[GraphNodeKey],
) -> BTreeMap<GraphNodeKey, usize>
where
  N: GraphNode,
  E: GraphEdge,
  D: Sync + Send,
{
  let mut counts = BTreeMap::new();
  for &node_key in postorder {
    let node = graph.get_node(node_key).unwrap();
    let child_keys = graph.child_keys_of(&node.read_arc());
    let count = if child_keys.is_empty() {
      1
    } else {
      child_keys.iter().map(|ck| counts[ck]).sum()
    };
    counts.insert(node_key, count);
  }
  counts
}

fn compute_heights<N, E, D>(graph: &Graph<N, E, D>, postorder: &[GraphNodeKey]) -> BTreeMap<GraphNodeKey, usize>
where
  N: GraphNode,
  E: GraphEdge,
  D: Sync + Send,
{
  let mut heights = BTreeMap::new();
  for &node_key in postorder {
    let node = graph.get_node(node_key).unwrap();
    let child_keys = graph.child_keys_of(&node.read_arc());
    let height = child_keys.iter().map(|ck| heights[ck] + 1).max().unwrap_or(0);
    heights.insert(node_key, height);
  }
  heights
}

fn compute_divergences<N, E, D>(
  graph: &Graph<N, E, D>,
  postorder: &[GraphNodeKey],
) -> BTreeMap<GraphNodeKey, OrderedFloat<f64>>
where
  N: GraphNode,
  E: GraphEdge + HasBranchLength,
  D: Sync + Send,
{
  let mut divergences: BTreeMap<GraphNodeKey, OrderedFloat<f64>> = BTreeMap::new();
  for &node_key in postorder {
    let node = graph.get_node(node_key).unwrap();
    let node = node.read_arc();
    let divergence = graph
      .children_of(&node)
      .iter()
      .map(|(child, edge)| {
        let child_key = child.read_arc().key();
        let edge_len = edge.read_arc().payload().read_arc().branch_length().unwrap_or(0.0);
        divergences[&child_key].0 + edge_len
      })
      .reduce(f64::max)
      .unwrap_or(0.0);
    divergences.insert(node_key, OrderedFloat(divergence));
  }
  divergences
}

fn compute_labels<N, E, D>(
  graph: &Graph<N, E, D>,
  postorder: &[GraphNodeKey],
) -> Result<BTreeMap<GraphNodeKey, String>, Report>
where
  N: GraphNode + Named,
  E: GraphEdge,
  D: Sync + Send,
{
  let mut labels: BTreeMap<GraphNodeKey, String> = BTreeMap::new();
  for &node_key in postorder {
    let node = graph.get_node(node_key).unwrap();
    let node = node.read_arc();
    let child_keys = graph.child_keys_of(&node);
    let label = if child_keys.is_empty() {
      node
        .payload()
        .read_arc()
        .name()
        .map(|name| name.as_ref().to_owned())
        .ok_or_else(|| make_report!("When ordering topology by labels: leaf node {} has no name", node_key))?
    } else {
      child_keys.iter().map(|ck| labels[ck].clone()).min().unwrap_or_default()
    };
    labels.insert(node_key, label);
  }
  Ok(labels)
}

fn compute_target_scores<N, E, D>(
  graph: &Graph<N, E, D>,
  postorder: &[GraphNodeKey],
  target_order: &[String],
  aggregate: TopologyOrderTargetAggregate,
) -> Result<BTreeMap<GraphNodeKey, TargetScore>, Report>
where
  N: GraphNode + Named,
  E: GraphEdge,
  D: Sync + Send,
{
  if target_order.is_empty() {
    return make_error!("When ordering topology: target-order mode requires a non-empty target order");
  }

  let position_of: BTreeMap<&str, usize> = target_order
    .iter()
    .enumerate()
    .map(|(i, label)| (label.as_str(), i))
    .collect();

  if position_of.len() != target_order.len() {
    let duplicate = target_order
      .iter()
      .duplicates()
      .next()
      .expect("target-order length mismatch guarantees a duplicate");
    return make_error!("When ordering topology: target order contains duplicate leaf label '{duplicate}'");
  }

  validate_target_order(graph, &position_of)?;

  match aggregate {
    TopologyOrderTargetAggregate::Mean => compute_target_scores_mean(graph, postorder, &position_of),
    TopologyOrderTargetAggregate::Median => compute_target_scores_median(graph, postorder, &position_of),
  }
}

fn compute_target_scores_mean<N, E, D>(
  graph: &Graph<N, E, D>,
  postorder: &[GraphNodeKey],
  position_of: &BTreeMap<&str, usize>,
) -> Result<BTreeMap<GraphNodeKey, TargetScore>, Report>
where
  N: GraphNode + Named,
  E: GraphEdge,
  D: Sync + Send,
{
  let mut scores: BTreeMap<GraphNodeKey, TargetScore> = BTreeMap::new();
  for &node_key in postorder {
    let node = graph.get_node(node_key).unwrap();
    let node = node.read_arc();
    let child_keys = graph.child_keys_of(&node);
    let score = if child_keys.is_empty() {
      let name = node
        .payload()
        .read_arc()
        .name()
        .map(|n| n.as_ref().to_owned())
        .ok_or_else(|| make_report!("When ordering topology by target order: leaf node {node_key} has no name"))?;
      let pos = *position_of.get(name.as_str()).ok_or_else(|| {
        make_report!("When ordering topology by target order: leaf '{name}' is absent from target order")
      })?;
      TargetScore {
        numerator: pos,
        denominator: 1,
      }
    } else {
      TargetScore {
        numerator: child_keys.iter().map(|ck| scores[ck].numerator).sum(),
        denominator: child_keys.iter().map(|ck| scores[ck].denominator).sum(),
      }
    };
    scores.insert(node_key, score);
  }
  Ok(scores)
}

fn compute_target_scores_median<N, E, D>(
  graph: &Graph<N, E, D>,
  postorder: &[GraphNodeKey],
  position_of: &BTreeMap<&str, usize>,
) -> Result<BTreeMap<GraphNodeKey, TargetScore>, Report>
where
  N: GraphNode + Named,
  E: GraphEdge,
  D: Sync + Send,
{
  let mut positions: BTreeMap<GraphNodeKey, Vec<usize>> = BTreeMap::new();
  let mut scores = BTreeMap::new();
  for &node_key in postorder {
    let node = graph.get_node(node_key).unwrap();
    let node = node.read_arc();
    let child_keys = graph.child_keys_of(&node);
    let pos = if child_keys.is_empty() {
      let name = node
        .payload()
        .read_arc()
        .name()
        .map(|n| n.as_ref().to_owned())
        .ok_or_else(|| make_report!("When ordering topology by target order: leaf node {node_key} has no name"))?;
      let p = *position_of.get(name.as_str()).ok_or_else(|| {
        make_report!("When ordering topology by target order: leaf '{name}' is absent from target order")
      })?;
      vec![p]
    } else {
      child_keys
        .iter()
        .flat_map(|ck| positions[ck].iter().copied())
        .sorted_unstable()
        .collect_vec()
    };
    let score = median_score(&pos);
    scores.insert(node_key, score);
    positions.insert(node_key, pos);
  }
  Ok(scores)
}

fn median_score(sorted_positions: &[usize]) -> TargetScore {
  let n = sorted_positions.len();
  let midpoint = n / 2;
  if n.is_multiple_of(2) {
    TargetScore {
      numerator: sorted_positions[midpoint - 1] + sorted_positions[midpoint],
      denominator: 2,
    }
  } else {
    TargetScore {
      numerator: sorted_positions[midpoint],
      denominator: 1,
    }
  }
}

fn validate_target_order<N, E, D>(graph: &Graph<N, E, D>, position_of: &BTreeMap<&str, usize>) -> Result<(), Report>
where
  N: GraphNode + Named,
  E: GraphEdge,
  D: Sync + Send,
{
  let mut final_labels = BTreeMap::new();
  for leaf in graph.get_leaves() {
    let leaf = leaf.read_arc();
    let label = leaf
      .payload()
      .read_arc()
      .name()
      .map(|name| name.as_ref().to_owned())
      .ok_or_else(|| make_report!("When validating target order: leaf node {} has no name", leaf.key()))?;
    if !position_of.contains_key(label.as_str()) {
      return make_error!("When validating target order: leaf '{label}' is absent from target order");
    }
    if let Some(previous_key) = final_labels.insert(label.clone(), leaf.key()) {
      return make_error!(
        "When validating target order: final leaf label '{label}' is duplicated by nodes {previous_key} and {}",
        leaf.key()
      );
    }
  }
  Ok(())
}

fn postorder_keys<N, E, D>(graph: &Graph<N, E, D>) -> Result<Vec<GraphNodeKey>, Report>
where
  N: GraphNode,
  E: GraphEdge,
  D: Sync + Send,
{
  let node_keys = graph
    .nodes
    .iter()
    .filter_map(|node| node.as_ref().map(|node| node.read_arc().key()))
    .collect_vec();

  let mut remaining_inbound = node_keys
    .iter()
    .map(|node_key| (*node_key, graph.degree_in(*node_key).unwrap_or(0)))
    .collect::<BTreeMap<_, _>>();

  let mut queue = remaining_inbound
    .iter()
    .filter_map(|(node_key, count)| (*count == 0).then_some(*node_key))
    .collect::<VecDeque<_>>();

  let mut ordered = Vec::with_capacity(node_keys.len());
  while let Some(node_key) = queue.pop_front() {
    ordered.push(node_key);
    let node = graph
      .get_node(node_key)
      .ok_or_else(|| make_report!("When computing topology order: Node {node_key} not found"))?;
    for child_key in graph.child_keys_of(&node.read_arc()) {
      let count = remaining_inbound
        .get_mut(&child_key)
        .ok_or_else(|| make_report!("When computing topology order: Node {child_key} not found"))?;
      *count -= 1;
      if *count == 0 {
        queue.push_back(child_key);
      }
    }
  }

  if ordered.len() != node_keys.len() {
    return make_error!("When ordering topology: graph contains a directed cycle");
  }

  ordered.reverse();
  Ok(ordered)
}
