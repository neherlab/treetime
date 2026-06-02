use crate::edge::{Edge, GraphEdge, GraphEdgeKey};
use crate::graph::Graph;
use crate::node::{GraphNode, GraphNodeKey, Named, Node};
use eyre::Report;
use itertools::Itertools;
use parking_lot::RwLock;
use serde::{Deserialize, Serialize};
use std::cmp::Ordering;
use std::collections::{BTreeMap, BTreeSet, VecDeque};
use std::sync::Arc;
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

  pub fn plan<N, E, D>(&self, graph: &Graph<N, E, D>) -> Result<TopologyOrderPlan, Report>
  where
    N: GraphNode + Named,
    E: GraphEdge,
    D: Sync + Send,
  {
    let mut state = TopologyOrderState::new(graph, self)?;
    let node_order = topological_order(graph)?;

    for node_key in node_order.into_iter().rev() {
      let node = graph
        .get_node(node_key)
        .ok_or_else(|| make_report!("When ordering topology: Node {node_key} not found"))?;
      state.compute_node(graph, &node.read_arc())?;
    }

    Ok(TopologyOrderPlan {
      roots: self.order_node_keys(&graph.roots, &state),
      leaves: self.order_node_keys(&graph.leaves, &state),
      outbound_edges: graph
        .nodes
        .iter()
        .filter_map(Option::as_ref)
        .map(|node| {
          let node = node.read_arc();
          (node.key(), self.order_edge_keys(graph, node.outbound(), &state))
        })
        .collect(),
    })
  }

  fn order_edge_keys<N, E, D>(
    &self,
    graph: &Graph<N, E, D>,
    edge_keys: &[GraphEdgeKey],
    state: &TopologyOrderState,
  ) -> Vec<GraphEdgeKey>
  where
    N: GraphNode,
    E: GraphEdge,
    D: Sync + Send,
  {
    if self.preset == TopologyOrderPreset::Keep {
      return edge_keys.to_vec();
    }

    edge_keys
      .iter()
      .copied()
      .enumerate()
      .sorted_by(|(lhs_pos, lhs), (rhs_pos, rhs)| {
        let lhs_target = graph.get_target_node_key(*lhs);
        let rhs_target = graph.get_target_node_key(*rhs);
        match (lhs_target, rhs_target) {
          (Ok(lhs_target), Ok(rhs_target)) => self
            .compare_nodes(lhs_target, rhs_target, state)
            .then(lhs_pos.cmp(rhs_pos)),
          _ => lhs_pos.cmp(rhs_pos),
        }
      })
      .map(|(_, edge_key)| edge_key)
      .collect_vec()
  }

  fn order_node_keys(&self, node_keys: &[GraphNodeKey], state: &TopologyOrderState) -> Vec<GraphNodeKey> {
    if self.preset == TopologyOrderPreset::Keep {
      return node_keys.to_vec();
    }

    node_keys
      .iter()
      .copied()
      .enumerate()
      .sorted_by(|(lhs_pos, lhs), (rhs_pos, rhs)| self.compare_nodes(*lhs, *rhs, state).then(lhs_pos.cmp(rhs_pos)))
      .map(|(_, node_key)| node_key)
      .collect_vec()
  }

  fn compare_nodes(&self, lhs: GraphNodeKey, rhs: GraphNodeKey, state: &TopologyOrderState) -> Ordering {
    let ordering = match self.preset {
      TopologyOrderPreset::Keep => Ordering::Equal,
      TopologyOrderPreset::DescendantCount | TopologyOrderPreset::DescendantCountReverse => {
        state.leaf_counts[&lhs].cmp(&state.leaf_counts[&rhs])
      },
      TopologyOrderPreset::Height | TopologyOrderPreset::HeightReverse => state.heights[&lhs].cmp(&state.heights[&rhs]),
      TopologyOrderPreset::Label | TopologyOrderPreset::LabelReverse => state.labels[&lhs].cmp(&state.labels[&rhs]),
      TopologyOrderPreset::TargetOrder | TopologyOrderPreset::TargetOrderReverse => {
        state.target_scores[&lhs].cmp(&state.target_scores[&rhs])
      },
    };

    if self.preset.is_reverse() {
      ordering.reverse()
    } else {
      ordering
    }
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
      Self::DescendantCountReverse | Self::HeightReverse | Self::LabelReverse | Self::TargetOrderReverse
    )
  }
}

#[derive(Copy, Clone, Debug, Default, Eq, PartialEq, Serialize, Deserialize)]
pub enum TopologyOrderTargetAggregate {
  #[default]
  Mean,
  Median,
}

#[derive(Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub struct TopologyOrderPlan {
  pub roots: Vec<GraphNodeKey>,
  pub leaves: Vec<GraphNodeKey>,
  pub outbound_edges: BTreeMap<GraphNodeKey, Vec<GraphEdgeKey>>,
}

impl TopologyOrderPlan {
  pub fn ordered_graph<N, E, D>(&self, graph: &Graph<N, E, D>) -> Result<Graph<N, E, D>, Report>
  where
    N: GraphNode,
    E: GraphEdge,
    D: Sync + Send,
  {
    let mut ordered = Graph {
      nodes: graph
        .nodes
        .iter()
        .map(|node| {
          node.as_ref().map(|node| {
            let node = node.read_arc();
            let mut cloned = Node::new(node.key(), node.payload().read_arc().clone());
            *cloned.outbound_mut() = self
              .outbound_edges
              .get(&node.key())
              .cloned()
              .unwrap_or_else(|| node.outbound().to_vec());
            *cloned.inbound_mut() = node.inbound().to_vec();
            Arc::new(RwLock::new(cloned))
          })
        })
        .collect(),
      edges: graph
        .edges
        .iter()
        .map(|edge| {
          edge.as_ref().map(|edge| {
            let edge = edge.read_arc();
            Arc::new(RwLock::new(Edge::new(
              edge.key(),
              edge.source(),
              edge.target(),
              edge.payload().read_arc().clone(),
            )))
          })
        })
        .collect(),
      roots: self.roots.clone(),
      leaves: self.leaves.clone(),
      data: graph.data(),
    };
    ordered.build()?;
    ordered.roots = self.roots.clone();
    ordered.leaves = self.leaves.clone();
    Ok(ordered)
  }
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

struct TopologyOrderState {
  leaf_sets: BTreeMap<GraphNodeKey, BTreeSet<GraphNodeKey>>,
  leaf_counts: BTreeMap<GraphNodeKey, usize>,
  heights: BTreeMap<GraphNodeKey, usize>,
  labels: BTreeMap<GraphNodeKey, String>,
  target_scores: BTreeMap<GraphNodeKey, TargetScore>,
  target_order: BTreeMap<String, usize>,
  target_aggregate: TopologyOrderTargetAggregate,
}

impl TopologyOrderState {
  fn new<N, E, D>(graph: &Graph<N, E, D>, spec: &TopologyOrderSpec) -> Result<Self, Report>
  where
    N: GraphNode + Named,
    E: GraphEdge,
    D: Sync + Send,
  {
    let target_order = spec
      .target_order
      .iter()
      .enumerate()
      .map(|(index, label)| (label.clone(), index))
      .collect();

    if spec.preset.is_target_order() && spec.target_order.is_empty() {
      return make_error!("When ordering topology: target-order mode requires a non-empty target order");
    }

    let state = Self {
      leaf_sets: BTreeMap::new(),
      leaf_counts: BTreeMap::new(),
      heights: BTreeMap::new(),
      labels: BTreeMap::new(),
      target_scores: BTreeMap::new(),
      target_order,
      target_aggregate: spec.target_aggregate,
    };

    if spec.preset.is_target_order() {
      validate_target_order(graph, &state.target_order)?;
    }

    Ok(state)
  }

  fn compute_node<N, E, D>(&mut self, graph: &Graph<N, E, D>, node: &Node<N>) -> Result<(), Report>
  where
    N: GraphNode + Named,
    E: GraphEdge,
    D: Sync + Send,
  {
    let child_keys = graph.child_keys_of(node);
    let leaf_set = if child_keys.is_empty() {
      BTreeSet::from([node.key()])
    } else {
      child_keys
        .iter()
        .flat_map(|child_key| self.leaf_sets[child_key].iter().copied())
        .collect()
    };

    let height = child_keys
      .iter()
      .map(|child_key| self.heights[child_key] + 1)
      .max()
      .unwrap_or(0);

    let label = if child_keys.is_empty() {
      node
        .payload()
        .read_arc()
        .name()
        .map(|name| name.as_ref().to_owned())
        .ok_or_else(|| make_report!("When ordering topology by labels: leaf node {} has no name", node.key()))?
    } else {
      child_keys
        .iter()
        .map(|child_key| self.labels[child_key].clone())
        .min()
        .unwrap_or_default()
    };

    let target_score = self.target_score(graph, node, &leaf_set)?;

    self.leaf_counts.insert(node.key(), leaf_set.len());
    self.leaf_sets.insert(node.key(), leaf_set);
    self.heights.insert(node.key(), height);
    self.labels.insert(node.key(), label);
    self.target_scores.insert(node.key(), target_score);

    Ok(())
  }

  fn target_score<N, E, D>(
    &self,
    graph: &Graph<N, E, D>,
    node: &Node<N>,
    leaf_set: &BTreeSet<GraphNodeKey>,
  ) -> Result<TargetScore, Report>
  where
    N: GraphNode + Named,
    E: GraphEdge,
    D: Sync + Send,
  {
    if self.target_order.is_empty() {
      return Ok(TargetScore {
        numerator: 0,
        denominator: 1,
      });
    }

    let mut positions = leaf_set
      .iter()
      .map(|leaf_key| {
        let leaf = graph
          .get_node(*leaf_key)
          .ok_or_else(|| make_report!("When ordering topology: leaf node {leaf_key} not found"))?;
        let label = leaf
          .read_arc()
          .payload()
          .read_arc()
          .name()
          .map(|name| name.as_ref().to_owned())
          .ok_or_else(|| make_report!("When ordering topology by target order: leaf node {leaf_key} has no name"))?;
        self.target_order.get(&label).copied().ok_or_else(|| {
          make_report!("When ordering topology by target order: leaf '{label}' is absent from target order")
        })
      })
      .collect::<Result<Vec<_>, _>>()?;

    positions.sort_unstable();

    match self.target_aggregate {
      TopologyOrderTargetAggregate::Mean => Ok(TargetScore {
        numerator: positions.iter().sum(),
        denominator: positions.len(),
      }),
      TopologyOrderTargetAggregate::Median => {
        let midpoint = positions.len() / 2;
        if positions.len() % 2 == 0 {
          Ok(TargetScore {
            numerator: positions[midpoint - 1] + positions[midpoint],
            denominator: 2,
          })
        } else {
          Ok(TargetScore {
            numerator: positions[midpoint],
            denominator: 1,
          })
        }
      },
    }
  }
}

fn validate_target_order<N, E, D>(graph: &Graph<N, E, D>, target_order: &BTreeMap<String, usize>) -> Result<(), Report>
where
  N: GraphNode + Named,
  E: GraphEdge,
  D: Sync + Send,
{
  for leaf in graph.get_leaves() {
    let leaf = leaf.read_arc();
    let label = leaf
      .payload()
      .read_arc()
      .name()
      .map(|name| name.as_ref().to_owned())
      .ok_or_else(|| make_report!("When validating target order: leaf node {} has no name", leaf.key()))?;
    if !target_order.contains_key(&label) {
      return make_error!("When validating target order: leaf '{label}' is absent from target order");
    }
  }
  Ok(())
}

fn topological_order<N, E, D>(graph: &Graph<N, E, D>) -> Result<Vec<GraphNodeKey>, Report>
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
      .ok_or_else(|| make_report!("When validating topology order: Node {node_key} not found"))?;
    for child_key in graph.child_keys_of(&node.read_arc()) {
      let count = remaining_inbound
        .get_mut(&child_key)
        .ok_or_else(|| make_report!("When validating topology order: Node {child_key} not found"))?;
      *count -= 1;
      if *count == 0 {
        queue.push_back(child_key);
      }
    }
  }

  if ordered.len() != node_keys.len() {
    return make_error!("When ordering topology: graph contains a directed cycle");
  }

  Ok(ordered)
}

#[cfg(test)]
mod tests {
  use super::*;
  use pretty_assertions::assert_eq;

  #[test]
  fn topology_order_descendant_count_sorts_children_ascending() -> Result<(), Report> {
    let graph = fixture_tree()?;
    let plan = TopologyOrderSpec::default().plan(&graph)?;
    let ordered = plan.ordered_graph(&graph)?;

    let actual = child_names(&ordered, "root")?;

    assert_eq!(vec!["A", "BC", "DEF"], actual);
    assert_eq!(vec!["DEF", "A", "BC"], child_names(&graph, "root")?);

    Ok(())
  }

  #[test]
  fn topology_order_descendant_count_reverse_sorts_children_descending() -> Result<(), Report> {
    let graph = fixture_tree()?;
    let plan = TopologyOrderSpec::descendant_count(true).plan(&graph)?;
    let ordered = plan.ordered_graph(&graph)?;

    let actual = child_names(&ordered, "root")?;

    assert_eq!(vec!["DEF", "BC", "A"], actual);

    Ok(())
  }

  #[test]
  fn topology_order_keep_preserves_outbound_order() -> Result<(), Report> {
    let graph = fixture_tree()?;
    let plan = TopologyOrderSpec::keep().plan(&graph)?;
    let ordered = plan.ordered_graph(&graph)?;

    let actual = child_names(&ordered, "root")?;

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

    graph.add_edge(root, right, TestEdge)?;
    graph.add_edge(root, left, TestEdge)?;
    graph.add_edge(left, shared, TestEdge)?;
    graph.add_edge(right, shared, TestEdge)?;
    graph.add_edge(right, right_only, TestEdge)?;
    graph.build()?;

    let plan = TopologyOrderSpec::default().plan(&graph)?;
    let ordered = plan.ordered_graph(&graph)?;

    let actual = child_names(&ordered, "root")?;

    assert_eq!(vec!["left", "right"], actual);

    Ok(())
  }

  #[test]
  fn topology_order_target_order_uses_requested_tip_order() -> Result<(), Report> {
    let graph = fixture_tree()?;
    let spec = TopologyOrderSpec {
      preset: TopologyOrderPreset::TargetOrder,
      target_order: vec!["D", "E", "F", "B", "C", "A"]
        .into_iter()
        .map(str::to_owned)
        .collect(),
      target_aggregate: TopologyOrderTargetAggregate::Mean,
    };
    let plan = spec.plan(&graph)?;
    let ordered = plan.ordered_graph(&graph)?;

    let actual = child_names(&ordered, "root")?;

    assert_eq!(vec!["DEF", "BC", "A"], actual);

    Ok(())
  }

  #[test]
  fn topology_order_rejects_cycles() -> Result<(), Report> {
    let mut graph = Graph::<TestNode, TestEdge, ()>::new();
    let a = graph.add_node(TestNode::new("A"));
    let b = graph.add_node(TestNode::new("B"));
    let c = graph.add_node(TestNode::new("C"));

    graph.add_edge(a, b, TestEdge)?;
    graph.add_edge(b, c, TestEdge)?;
    graph.add_edge(c, a, TestEdge)?;
    graph.build()?;

    let err = TopologyOrderSpec::default().plan(&graph).unwrap_err();

    assert!(err.to_string().contains("directed cycle"));

    Ok(())
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

    graph.add_edge(root, def, TestEdge)?;
    graph.add_edge(root, tip_a, TestEdge)?;
    graph.add_edge(root, bc, TestEdge)?;
    graph.add_edge(def, tip_d, TestEdge)?;
    graph.add_edge(def, tip_e, TestEdge)?;
    graph.add_edge(def, tip_f, TestEdge)?;
    graph.add_edge(bc, tip_b, TestEdge)?;
    graph.add_edge(bc, tip_c, TestEdge)?;
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

  #[derive(Clone, Debug, Eq, PartialEq)]
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

  #[derive(Clone, Debug, Eq, PartialEq)]
  struct TestEdge;

  impl GraphEdge for TestEdge {}
}
