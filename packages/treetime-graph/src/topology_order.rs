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

#[cfg(test)]
mod tests {
  use super::*;
  use pretty_assertions::assert_eq;

  #[test]
  fn topology_order_descendant_count_sorts_children_ascending() -> Result<(), Report> {
    let mut graph = fixture_tree()?;
    let original = child_names(&graph, "root")?;
    TopologyOrderSpec::default().apply(&mut graph)?;
    let ordered = &graph;

    let actual = child_names(ordered, "root")?;

    assert_eq!(vec!["A", "BC", "DEF"], actual);
    assert_eq!(vec!["DEF", "A", "BC"], original);

    Ok(())
  }

  #[test]
  fn topology_order_descendant_count_reverse_sorts_children_descending() -> Result<(), Report> {
    let mut graph = fixture_tree()?;
    TopologyOrderSpec::descendant_count(true).apply(&mut graph)?;
    let ordered = &graph;

    let actual = child_names(ordered, "root")?;

    assert_eq!(vec!["DEF", "BC", "A"], actual);

    Ok(())
  }

  #[test]
  fn topology_order_keep_preserves_outbound_order() -> Result<(), Report> {
    let mut graph = fixture_tree()?;
    TopologyOrderSpec::keep().apply(&mut graph)?;
    let ordered = &graph;

    let actual = child_names(ordered, "root")?;

    assert_eq!(vec!["DEF", "A", "BC"], actual);

    Ok(())
  }

  #[test]
  fn topology_order_dag_counts_shared_descendant_once_per_child() -> Result<(), Report> {
    let mut graph = Graph::<TestNode, TestEdge, ()>::new();
    let root = graph.add_node(TestNode::new("root"));
    let left = graph.add_node(TestNode::new("left"));
    let right = graph.add_node(TestNode::new("right"));
    let shared = graph.add_node(TestNode::new("shared"));
    let right_only = graph.add_node(TestNode::new("right_only"));

    graph.add_edge(root, right, TestEdge::new())?;
    graph.add_edge(root, left, TestEdge::new())?;
    graph.add_edge(left, shared, TestEdge::new())?;
    graph.add_edge(right, shared, TestEdge::new())?;
    graph.add_edge(right, right_only, TestEdge::new())?;
    graph.build()?;

    TopologyOrderSpec::default().apply(&mut graph)?;
    let ordered = &graph;

    let actual = child_names(ordered, "root")?;

    assert_eq!(vec!["left", "right"], actual);

    Ok(())
  }

  #[test]
  fn topology_order_target_order_uses_requested_tip_order() -> Result<(), Report> {
    let mut graph = fixture_tree()?;
    let spec = TopologyOrderSpec {
      preset: TopologyOrderPreset::TargetOrder,
      target_order: vec!["D", "E", "F", "B", "C", "A"]
        .into_iter()
        .map(str::to_owned)
        .collect(),
      target_aggregate: TopologyOrderTargetAggregate::Mean,
    };
    spec.apply(&mut graph)?;
    let ordered = &graph;

    let actual = child_names(ordered, "root")?;

    assert_eq!(vec!["DEF", "BC", "A"], actual);

    Ok(())
  }

  #[test]
  fn topology_order_rejects_cycles() -> Result<(), Report> {
    let mut graph = Graph::<TestNode, TestEdge, ()>::new();
    let a = graph.add_node(TestNode::new("A"));
    let b = graph.add_node(TestNode::new("B"));
    let c = graph.add_node(TestNode::new("C"));

    graph.add_edge(a, b, TestEdge::new())?;
    graph.add_edge(b, c, TestEdge::new())?;
    graph.add_edge(c, a, TestEdge::new())?;
    graph.build()?;

    let err = TopologyOrderSpec::default().apply(&mut graph).unwrap_err();

    assert!(err.to_string().contains("directed cycle"));

    Ok(())
  }

  #[test]
  fn topology_order_height_sorts_by_subtree_depth() -> Result<(), Report> {
    let mut graph = fixture_deep_tree()?;
    let spec = TopologyOrderSpec {
      preset: TopologyOrderPreset::Height,
      ..TopologyOrderSpec::default()
    };
    spec.apply(&mut graph)?;
    let ordered = &graph;

    let actual = child_names(ordered, "root")?;

    // A is leaf (height 0), shallow has height 1, deep has height 2
    assert_eq!(vec!["A", "shallow", "deep"], actual);

    Ok(())
  }

  #[test]
  fn topology_order_height_reverse_sorts_deepest_first() -> Result<(), Report> {
    let mut graph = fixture_deep_tree()?;
    let spec = TopologyOrderSpec {
      preset: TopologyOrderPreset::HeightReverse,
      ..TopologyOrderSpec::default()
    };
    spec.apply(&mut graph)?;
    let ordered = &graph;

    let actual = child_names(ordered, "root")?;

    assert_eq!(vec!["deep", "shallow", "A"], actual);

    Ok(())
  }

  #[test]
  fn topology_order_divergence_sorts_by_total_branch_length() -> Result<(), Report> {
    // root -> [short(0.1) -> [A(0.1), B(0.1)], long(0.5) -> [C(0.2)], D(0.3)]
    // short: max divergence = 0.1 + 0.1 = 0.2
    // long:  max divergence = 0.5 + 0.2 = 0.7
    // D:     leaf, divergence = 0.0
    let mut graph = fixture_branch_length_tree()?;
    let spec = TopologyOrderSpec {
      preset: TopologyOrderPreset::Divergence,
      ..TopologyOrderSpec::default()
    };
    spec.apply(&mut graph)?;
    let ordered = &graph;

    let actual = child_names(ordered, "root")?;

    assert_eq!(vec!["D", "short", "long"], actual);

    Ok(())
  }

  #[test]
  fn topology_order_divergence_reverse_sorts_longest_first() -> Result<(), Report> {
    let mut graph = fixture_branch_length_tree()?;
    let spec = TopologyOrderSpec {
      preset: TopologyOrderPreset::DivergenceReverse,
      ..TopologyOrderSpec::default()
    };
    spec.apply(&mut graph)?;
    let ordered = &graph;

    let actual = child_names(ordered, "root")?;

    assert_eq!(vec!["long", "short", "D"], actual);

    Ok(())
  }

  #[test]
  fn topology_order_label_sorts_alphabetically() -> Result<(), Report> {
    let mut graph = fixture_tree()?;
    let spec = TopologyOrderSpec {
      preset: TopologyOrderPreset::Label,
      ..TopologyOrderSpec::default()
    };
    spec.apply(&mut graph)?;
    let ordered = &graph;

    let actual = child_names(ordered, "root")?;

    // A has label "A", BC has min label "B", DEF has min label "D"
    assert_eq!(vec!["A", "BC", "DEF"], actual);

    Ok(())
  }

  #[test]
  fn topology_order_label_reverse_sorts_descending() -> Result<(), Report> {
    let mut graph = fixture_tree()?;
    let spec = TopologyOrderSpec {
      preset: TopologyOrderPreset::LabelReverse,
      ..TopologyOrderSpec::default()
    };
    spec.apply(&mut graph)?;
    let ordered = &graph;

    let actual = child_names(ordered, "root")?;

    assert_eq!(vec!["DEF", "BC", "A"], actual);

    Ok(())
  }

  #[test]
  fn topology_order_target_order_median_uses_median_position() -> Result<(), Report> {
    let mut graph = fixture_tree()?;
    // Target order: D=0, E=1, F=2, B=3, C=4, A=5
    // DEF median of [0,1,2] = 1, BC median of [3,4] = 3.5, A = 5
    let spec = TopologyOrderSpec {
      preset: TopologyOrderPreset::TargetOrder,
      target_order: vec!["D", "E", "F", "B", "C", "A"]
        .into_iter()
        .map(str::to_owned)
        .collect(),
      target_aggregate: TopologyOrderTargetAggregate::Median,
    };
    spec.apply(&mut graph)?;
    let ordered = &graph;

    let actual = child_names(ordered, "root")?;

    assert_eq!(vec!["DEF", "BC", "A"], actual);

    Ok(())
  }

  #[test]
  fn topology_order_propagates_through_nested_levels() -> Result<(), Report> {
    let mut graph = fixture_deep_tree()?;
    TopologyOrderSpec::default().apply(&mut graph)?;
    let ordered = &graph;

    let deep_children = child_names(ordered, "deep")?;
    // mid (2 leaves) vs D (1 leaf): D first
    assert_eq!(vec!["D", "mid"], deep_children);

    let mid_children = child_names(ordered, "mid")?;
    assert_eq!(vec!["E", "F"], mid_children);

    Ok(())
  }

  #[test]
  fn topology_order_target_order_reverse_inverts_order() -> Result<(), Report> {
    let mut graph = fixture_tree()?;
    let spec = TopologyOrderSpec {
      preset: TopologyOrderPreset::TargetOrderReverse,
      target_order: vec!["D", "E", "F", "B", "C", "A"]
        .into_iter()
        .map(str::to_owned)
        .collect(),
      target_aggregate: TopologyOrderTargetAggregate::Mean,
    };
    spec.apply(&mut graph)?;
    let ordered = &graph;

    let actual = child_names(ordered, "root")?;

    assert_eq!(vec!["A", "BC", "DEF"], actual);

    Ok(())
  }

  #[test]
  fn topology_order_target_order_rejects_empty() {
    let mut graph = fixture_tree().unwrap();
    let spec = TopologyOrderSpec {
      preset: TopologyOrderPreset::TargetOrder,
      target_order: vec![],
      target_aggregate: TopologyOrderTargetAggregate::Mean,
    };
    let err = spec.apply(&mut graph).unwrap_err();
    assert!(err.to_string().contains("non-empty target order"));
  }

  #[test]
  fn topology_order_target_order_rejects_duplicate_ranking_labels() {
    let mut graph = fixture_tree().unwrap();
    let spec = TopologyOrderSpec {
      preset: TopologyOrderPreset::TargetOrder,
      target_order: vec!["A", "B", "B", "C", "D", "E", "F"]
        .into_iter()
        .map(str::to_owned)
        .collect(),
      target_aggregate: TopologyOrderTargetAggregate::Mean,
    };

    let error = spec.apply(&mut graph).unwrap_err();

    assert!(error.to_string().contains("duplicate leaf label 'B'"));
  }

  #[test]
  fn topology_order_target_order_rejects_duplicate_final_leaf_labels() -> Result<(), Report> {
    let mut graph = fixture_tree()?;
    graph
      .get_node(find_node(&graph, "C")?)
      .expect("fixture node C must exist")
      .write_arc()
      .payload()
      .write_arc()
      .0 = "B".to_owned();
    let spec = TopologyOrderSpec {
      preset: TopologyOrderPreset::TargetOrder,
      target_order: vec!["A", "B", "D", "E", "F"].into_iter().map(str::to_owned).collect(),
      target_aggregate: TopologyOrderTargetAggregate::Mean,
    };

    let error = spec.apply(&mut graph).unwrap_err();

    assert!(error.to_string().contains("final leaf label 'B' is duplicated"));
    Ok(())
  }

  #[test]
  fn topology_order_target_order_ignores_absent_ranking_labels() -> Result<(), Report> {
    let mut graph = fixture_tree()?;
    let spec = TopologyOrderSpec {
      preset: TopologyOrderPreset::TargetOrder,
      target_order: vec!["removed", "D", "E", "F", "B", "C", "A"]
        .into_iter()
        .map(str::to_owned)
        .collect(),
      target_aggregate: TopologyOrderTargetAggregate::Mean,
    };

    spec.apply(&mut graph)?;

    assert_eq!(vec!["DEF", "BC", "A"], child_names(&graph, "root")?);
    Ok(())
  }

  #[test]
  fn topology_order_is_idempotent_and_preserves_graph_data_identity() -> Result<(), Report> {
    let graph = fixture_tree()?;
    let mut graph = graph.map_data(NonCloneData);
    let data = std::ptr::from_ref(graph.data());

    TopologyOrderSpec::default().apply(&mut graph)?;
    let first = child_names_with_data(&graph, "root")?;
    TopologyOrderSpec::default().apply(&mut graph)?;
    let second = child_names_with_data(&graph, "root")?;

    assert_eq!(first, second);
    assert!(std::ptr::eq(data, graph.data()));
    Ok(())
  }

  /// root -> [deep -> [D, mid -> [E, F]], shallow -> [B, C], A]
  fn fixture_deep_tree() -> Result<Graph<TestNode, TestEdge, ()>, Report> {
    let mut graph = Graph::<TestNode, TestEdge, ()>::new();
    let root = graph.add_node(TestNode::new("root"));
    let deep = graph.add_node(TestNode::new("deep"));
    let tip_a = graph.add_node(TestNode::new("A"));
    let shallow = graph.add_node(TestNode::new("shallow"));
    let mid = graph.add_node(TestNode::new("mid"));
    let tip_d = graph.add_node(TestNode::new("D"));
    let tip_e = graph.add_node(TestNode::new("E"));
    let tip_f = graph.add_node(TestNode::new("F"));
    let tip_b = graph.add_node(TestNode::new("B"));
    let tip_c = graph.add_node(TestNode::new("C"));

    graph.add_edge(root, deep, TestEdge::new())?;
    graph.add_edge(root, tip_a, TestEdge::new())?;
    graph.add_edge(root, shallow, TestEdge::new())?;
    graph.add_edge(deep, tip_d, TestEdge::new())?;
    graph.add_edge(deep, mid, TestEdge::new())?;
    graph.add_edge(mid, tip_e, TestEdge::new())?;
    graph.add_edge(mid, tip_f, TestEdge::new())?;
    graph.add_edge(shallow, tip_b, TestEdge::new())?;
    graph.add_edge(shallow, tip_c, TestEdge::new())?;
    graph.build()?;

    Ok(graph)
  }

  /// root -> [short(0.1) -> [A(0.1), B(0.1)], long(0.5) -> [C(0.2)], D(0.3)]
  fn fixture_branch_length_tree() -> Result<Graph<TestNode, TestEdge, ()>, Report> {
    let mut graph = Graph::<TestNode, TestEdge, ()>::new();
    let root = graph.add_node(TestNode::new("root"));
    let short = graph.add_node(TestNode::new("short"));
    let long = graph.add_node(TestNode::new("long"));
    let tip_a = graph.add_node(TestNode::new("A"));
    let tip_b = graph.add_node(TestNode::new("B"));
    let tip_c = graph.add_node(TestNode::new("C"));
    let tip_d = graph.add_node(TestNode::new("D"));

    graph.add_edge(root, short, TestEdge::with_length(0.1))?;
    graph.add_edge(root, long, TestEdge::with_length(0.5))?;
    graph.add_edge(root, tip_d, TestEdge::with_length(0.3))?;
    graph.add_edge(short, tip_a, TestEdge::with_length(0.1))?;
    graph.add_edge(short, tip_b, TestEdge::with_length(0.1))?;
    graph.add_edge(long, tip_c, TestEdge::with_length(0.2))?;
    graph.build()?;

    Ok(graph)
  }

  fn fixture_tree() -> Result<Graph<TestNode, TestEdge, ()>, Report> {
    let mut graph = Graph::<TestNode, TestEdge, ()>::new();
    let root = graph.add_node(TestNode::new("root"));
    let def = graph.add_node(TestNode::new("DEF"));
    let tip_a = graph.add_node(TestNode::new("A"));
    let bc = graph.add_node(TestNode::new("BC"));
    let tip_d = graph.add_node(TestNode::new("D"));
    let tip_e = graph.add_node(TestNode::new("E"));
    let tip_f = graph.add_node(TestNode::new("F"));
    let tip_b = graph.add_node(TestNode::new("B"));
    let tip_c = graph.add_node(TestNode::new("C"));

    graph.add_edge(root, def, TestEdge::new())?;
    graph.add_edge(root, tip_a, TestEdge::new())?;
    graph.add_edge(root, bc, TestEdge::new())?;
    graph.add_edge(def, tip_d, TestEdge::new())?;
    graph.add_edge(def, tip_e, TestEdge::new())?;
    graph.add_edge(def, tip_f, TestEdge::new())?;
    graph.add_edge(bc, tip_b, TestEdge::new())?;
    graph.add_edge(bc, tip_c, TestEdge::new())?;
    graph.build()?;

    Ok(graph)
  }

  fn child_names(graph: &Graph<TestNode, TestEdge, ()>, parent_name: &str) -> Result<Vec<String>, Report> {
    let parent_key = find_node(graph, parent_name)?;
    let parent = graph
      .get_node(parent_key)
      .ok_or_else(|| make_report!("Node {parent_key} not found"))?;
    Ok(
      graph
        .children_of(&parent.read_arc())
        .into_iter()
        .map(|(node, _)| node.read_arc().payload().read_arc().0.clone())
        .collect_vec(),
    )
  }

  fn find_node(graph: &Graph<TestNode, TestEdge, ()>, name: &str) -> Result<GraphNodeKey, Report> {
    graph
      .find_node(|node| node.0 == name)
      .ok_or_else(|| make_report!("Node '{name}' not found"))
  }

  #[derive(Debug, Eq, PartialEq)]
  struct TestNode(String);

  impl TestNode {
    fn new(name: &str) -> Self {
      Self(name.to_owned())
    }
  }

  impl GraphNode for TestNode {}

  impl Named for TestNode {
    fn name(&self) -> Option<impl AsRef<str>> {
      Some(&self.0)
    }

    fn set_name(&mut self, name: Option<impl AsRef<str>>) {
      self.0 = name.map(|name| name.as_ref().to_owned()).unwrap_or_default();
    }
  }

  #[derive(Debug, PartialEq)]
  struct TestEdge {
    branch_length: Option<f64>,
  }

  impl TestEdge {
    fn new() -> Self {
      Self { branch_length: None }
    }

    fn with_length(len: f64) -> Self {
      Self {
        branch_length: Some(len),
      }
    }
  }

  impl GraphEdge for TestEdge {}

  impl HasBranchLength for TestEdge {
    fn branch_length(&self) -> Option<f64> {
      self.branch_length
    }

    fn set_branch_length(&mut self, branch_length: Option<f64>) {
      self.branch_length = branch_length;
    }
  }

  #[derive(Debug)]
  struct NonCloneData;

  fn child_names_with_data<D: Send + Sync>(
    graph: &Graph<TestNode, TestEdge, D>,
    parent_name: &str,
  ) -> Result<Vec<String>, Report> {
    let parent_key = graph
      .find_node(|node| node.0 == parent_name)
      .ok_or_else(|| make_report!("Node '{parent_name}' not found"))?;
    let parent = graph
      .get_node(parent_key)
      .ok_or_else(|| make_report!("Node {parent_key} not found"))?;
    Ok(
      graph
        .children_of(&parent.read_arc())
        .into_iter()
        .map(|(node, _)| node.read_arc().payload().read_arc().0.clone())
        .collect_vec(),
    )
  }
}
