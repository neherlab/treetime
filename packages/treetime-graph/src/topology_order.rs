use crate::assign_node_names::node_name_or_key;
use crate::edge::GraphEdgeKey;
use crate::graph::Graph;
use crate::node::GraphNodeKey;
use deser::{Deserialize, Serialize};
use eyre::Report;
use itertools::Itertools;
use ordered_float::OrderedFloat;
use std::cmp::Ordering;
use std::collections::{BTreeMap, VecDeque};
use treetime_utils::{make_error, make_report};

#[derive(Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub struct TopologyOrderSpec {
  pub preset: TopologyOrderPreset,
  pub target_order: BTreeMap<GraphNodeKey, usize>,
  pub target_aggregate: TopologyOrderTargetAggregate,
}

impl Default for TopologyOrderSpec {
  fn default() -> Self {
    Self {
      preset: TopologyOrderPreset::DescendantCount,
      target_order: BTreeMap::new(),
      target_aggregate: TopologyOrderTargetAggregate::Mean,
    }
  }
}

impl TopologyOrderSpec {
  #[allow(
    clippy::expect_used,
    reason = "expect on a value an upstream invariant guarantees is present"
  )]
  pub fn apply(
    &self,
    graph: &mut Graph,
    names: &BTreeMap<GraphNodeKey, Option<String>>,
    branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  ) -> Result<(), Report> {
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
          let keys = compute_divergences(graph, &postorder, branch_lengths);
          build_order(graph, &keys, reverse)
        },
        TopologyOrderPreset::Label | TopologyOrderPreset::LabelReverse => {
          let keys = compute_labels(graph, &postorder, names)?;
          build_order(graph, &keys, reverse)
        },
        TopologyOrderPreset::TargetOrder | TopologyOrderPreset::TargetOrderReverse => {
          let keys = compute_target_scores(graph, &postorder, names, &self.target_order, self.target_aggregate)?;
          build_order(graph, &keys, reverse)
        },
      }
    }?;

    for &node_key in order.outbound_edges.keys() {
      if graph.get_node(node_key).is_none() {
        return make_error!("Node {node_key} disappeared while applying topology order");
      }
    }

    for (node_key, outbound_edges) in order.outbound_edges {
      *graph
        .get_node_mut(node_key)
        .expect("Node presence validated above")
        .outbound_mut() = outbound_edges;
    }
    graph.roots = order.roots;
    graph.leaves = order.leaves;
    Ok(())
  }
}

#[derive(Copy, Clone, Debug, Default, Eq, PartialEq, Serialize, Deserialize)]
#[deser(rename_all = "kebab-case")]
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

fn build_order_unmodified(graph: &Graph) -> Result<TopologyOrder, Report> {
  Ok(TopologyOrder {
    roots: graph.roots.clone(),
    leaves: graph.leaves.clone(),
    outbound_edges: graph
      .get_nodes()
      .map(|node| (node.key(), node.outbound().to_vec()))
      .collect(),
  })
}

fn build_order<K: Ord>(
  graph: &Graph,
  keys: &BTreeMap<GraphNodeKey, K>,
  reverse: bool,
) -> Result<TopologyOrder, Report> {
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
      .get_nodes()
      .map(|node| (node.key(), sort_edge_keys(node.outbound())))
      .collect(),
  })
}

#[derive(Debug, Eq, PartialEq)]
struct TopologyOrder {
  roots: Vec<GraphNodeKey>,
  leaves: Vec<GraphNodeKey>,
  outbound_edges: BTreeMap<GraphNodeKey, Vec<GraphEdgeKey>>,
}

#[allow(
  clippy::unwrap_used,
  reason = "unwrap on a value an upstream invariant guarantees is present"
)]
fn compute_descendant_counts(graph: &Graph, postorder: &[GraphNodeKey]) -> BTreeMap<GraphNodeKey, usize> {
  let mut counts = BTreeMap::new();
  for &node_key in postorder {
    let node = graph.get_node(node_key).unwrap();
    let child_keys = graph.children_keys_of(node).map(|(key, _)| key).collect_vec();
    let count = if child_keys.is_empty() {
      1
    } else {
      child_keys.iter().map(|ck| counts[ck]).sum()
    };
    counts.insert(node_key, count);
  }
  counts
}

#[allow(
  clippy::unwrap_used,
  reason = "unwrap on a value an upstream invariant guarantees is present"
)]
fn compute_heights(graph: &Graph, postorder: &[GraphNodeKey]) -> BTreeMap<GraphNodeKey, usize> {
  let mut heights = BTreeMap::new();
  for &node_key in postorder {
    let node = graph.get_node(node_key).unwrap();
    let child_keys = graph.children_keys_of(node).map(|(key, _)| key).collect_vec();
    let height = child_keys.iter().map(|ck| heights[ck] + 1).max().unwrap_or(0);
    heights.insert(node_key, height);
  }
  heights
}

#[allow(
  clippy::unwrap_used,
  reason = "unwrap on a value an upstream invariant guarantees is present"
)]
fn compute_divergences(
  graph: &Graph,
  postorder: &[GraphNodeKey],
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
) -> BTreeMap<GraphNodeKey, OrderedFloat<f64>> {
  let mut divergences: BTreeMap<GraphNodeKey, OrderedFloat<f64>> = BTreeMap::new();
  for &node_key in postorder {
    let node = graph.get_node(node_key).unwrap();
    let divergence = graph
      .children_of(node)
      .map(|(child, edge)| {
        let child_key = child.key();
        let edge_len = branch_lengths[&edge.key()].unwrap_or(0.0);
        divergences[&child_key].0 + edge_len
      })
      .reduce(f64::max)
      .unwrap_or(0.0);
    divergences.insert(node_key, OrderedFloat(divergence));
  }
  divergences
}

#[allow(
  clippy::unwrap_used,
  reason = "unwrap on a value an upstream invariant guarantees is present"
)]
fn compute_labels(
  graph: &Graph,
  postorder: &[GraphNodeKey],
  names: &BTreeMap<GraphNodeKey, Option<String>>,
) -> Result<BTreeMap<GraphNodeKey, String>, Report> {
  let mut labels: BTreeMap<GraphNodeKey, String> = BTreeMap::new();
  for &node_key in postorder {
    let node = graph.get_node(node_key).unwrap();
    let child_keys = graph.children_keys_of(node).map(|(key, _)| key).collect_vec();
    let label = if child_keys.is_empty() {
      names[&node_key]
        .clone()
        .ok_or_else(|| make_report!("When ordering topology by labels: leaf node {} has no name", node_key))?
    } else {
      child_keys.iter().map(|ck| labels[ck].clone()).min().unwrap_or_default()
    };
    labels.insert(node_key, label);
  }
  Ok(labels)
}

fn compute_target_scores(
  graph: &Graph,
  postorder: &[GraphNodeKey],
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  target_order: &BTreeMap<GraphNodeKey, usize>,
  aggregate: TopologyOrderTargetAggregate,
) -> Result<BTreeMap<GraphNodeKey, TargetScore>, Report> {
  if target_order.is_empty() {
    return make_error!("When ordering topology: target-order mode requires a non-empty target order");
  }

  validate_target_order(graph, names, target_order)?;

  match aggregate {
    TopologyOrderTargetAggregate::Mean => compute_target_scores_mean(graph, postorder, target_order),
    TopologyOrderTargetAggregate::Median => compute_target_scores_median(graph, postorder, target_order),
  }
}

#[derive(Copy, Clone, Debug, Default, Eq, PartialEq, Serialize, Deserialize)]
#[deser(rename_all = "kebab-case")]
pub enum TopologyOrderTargetAggregate {
  #[default]
  Mean,
  Median,
}

#[allow(
  clippy::unwrap_used,
  reason = "unwrap on a value an upstream invariant guarantees is present"
)]
fn compute_target_scores_mean(
  graph: &Graph,
  postorder: &[GraphNodeKey],
  target_order: &BTreeMap<GraphNodeKey, usize>,
) -> Result<BTreeMap<GraphNodeKey, TargetScore>, Report> {
  let mut scores: BTreeMap<GraphNodeKey, TargetScore> = BTreeMap::new();
  for &node_key in postorder {
    let node = graph.get_node(node_key).unwrap();
    let child_keys = graph.children_keys_of(node).map(|(key, _)| key).collect_vec();
    let score = if child_keys.is_empty() {
      let pos = target_position(target_order, node_key)?;
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

#[allow(
  clippy::unwrap_used,
  reason = "unwrap on a value an upstream invariant guarantees is present"
)]
fn compute_target_scores_median(
  graph: &Graph,
  postorder: &[GraphNodeKey],
  target_order: &BTreeMap<GraphNodeKey, usize>,
) -> Result<BTreeMap<GraphNodeKey, TargetScore>, Report> {
  let mut positions: BTreeMap<GraphNodeKey, Vec<usize>> = BTreeMap::new();
  let mut scores = BTreeMap::new();
  for &node_key in postorder {
    let node = graph.get_node(node_key).unwrap();
    let child_keys = graph.children_keys_of(node).map(|(key, _)| key).collect_vec();
    let pos = if child_keys.is_empty() {
      vec![target_position(target_order, node_key)?]
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

#[allow(clippy::integer_division, reason = "integer division is the intended floor division")]
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

#[derive(Clone, Debug, Eq, PartialEq)]
struct TargetScore {
  numerator: usize,
  denominator: usize,
}

impl Ord for TargetScore {
  #[allow(
    clippy::as_conversions,
    reason = "count/index numeric cast is exact for the domain range"
  )]
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

fn validate_target_order(
  graph: &Graph,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  target_order: &BTreeMap<GraphNodeKey, usize>,
) -> Result<(), Report> {
  match graph
    .get_leaves()
    .map(|leaf| leaf.key())
    .find(|key| !target_order.contains_key(key))
  {
    Some(key) => make_error!(
      "When validating target order: leaf '{}' is absent from target order",
      node_name_or_key(key, names[&key].as_deref())
    ),
    None => Ok(()),
  }
}

fn target_position(target_order: &BTreeMap<GraphNodeKey, usize>, key: GraphNodeKey) -> Result<usize, Report> {
  target_order
    .get(&key)
    .copied()
    .ok_or_else(|| make_report!("When ordering topology by target order: leaf {key} is absent from target order"))
}

fn postorder_keys(graph: &Graph) -> Result<Vec<GraphNodeKey>, Report> {
  let node_keys = graph.node_keys().collect_vec();

  let mut remaining_inbound = node_keys
    .iter()
    .map(|node_key| Ok((*node_key, graph.degree_in(*node_key)?)))
    .collect::<Result<BTreeMap<_, _>, Report>>()?;

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
    for (child_key, _) in graph.children_keys_of(node) {
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
