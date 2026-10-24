use crate::build::mat::Mat;
use crate::build::query::Query;
use crate::build::references::Reaches;
use crate::seq::mutation::Sub;
use eyre::Report;
use std::cmp::Reverse;
use std::collections::{BTreeMap, BinaryHeap};
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::node::GraphNodeKey;
use treetime_primitives::Seq;

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Target {
  Node(GraphNodeKey),
  Edge(GraphEdgeKey),
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Placement {
  pub target: Target,
  pub cost: usize,
  pub n_best: usize,
  pub visited: usize,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct EdgeScore {
  delta: isize,
  shared: usize,
}

pub(crate) fn score_edge(subs: &[Sub], query: &Query, reference: &Seq) -> EdgeScore {
  let mut score = EdgeScore { delta: 0, shared: 0 };
  for sub in subs {
    let Some(states) = query.states(sub.pos(), reference) else {
      continue;
    };
    match (states.contains(sub.qry()), states.contains(sub.reff())) {
      (true, false) => {
        score.delta -= 1;
        score.shared += 1;
      },
      (false, true) => score.delta += 1,
      _ => {},
    }
  }
  score
}

pub(crate) fn shared_positions(subs: &[Sub], query: &Query, reference: &Seq) -> Vec<usize> {
  subs
    .iter()
    .filter(|sub| {
      query
        .states(sub.pos(), reference)
        .is_some_and(|states| states.contains(sub.qry()) && !states.contains(sub.reff()))
    })
    .map(Sub::pos)
    .collect()
}

struct Best {
  target: Option<Target>,
  cost: isize,
  n_best: usize,
  prefer: Option<Target>,
}

impl Best {
  fn new(prefer: Option<Target>) -> Self {
    Self {
      target: None,
      cost: isize::MAX,
      n_best: 0,
      prefer,
    }
  }

  fn offer(&mut self, target: Target, cost: isize) {
    let rank = |t: Target| (Some(t) != self.prefer, matches!(t, Target::Edge(_)), t);
    if cost < self.cost {
      self.cost = cost;
      self.target = Some(target);
      self.n_best = 1;
    } else if cost == self.cost {
      self.n_best += 1;
      if self.target.is_none_or(|current| rank(target) < rank(current)) {
        self.target = Some(target);
      }
    }
  }

  fn finish(self, visited: usize) -> Result<Placement, Report> {
    let target = self
      .target
      .ok_or_else(|| crate::make_internal_report!("Placement search evaluated no candidate"))?;
    Ok(Placement {
      target,
      cost: usize::try_from(self.cost)?,
      n_best: self.n_best,
      visited,
    })
  }
}

fn offer_node(mat: &Mat, best: &mut Best, node: GraphNodeKey, cost: isize) {
  if !mat.is_leaf(node) {
    best.offer(Target::Node(node), cost);
  }
}

fn offer_edge(mat: &Mat, best: &mut Best, edge: GraphEdgeKey, parent_cost: isize, score: EdgeScore) {
  if score.shared > 0 && score.shared < mat.edge_subs(edge).len() {
    best.offer(Target::Edge(edge), parent_cost - signed(score.shared));
  } else if score.shared > 0 {
    let (_, child) = mat.edge_endpoints(edge);
    if mat.is_leaf(child) {
      best.offer(Target::Edge(edge), parent_cost - signed(score.shared));
    }
  }
}

fn signed(value: usize) -> isize {
  isize::try_from(value).unwrap_or(isize::MAX)
}

pub fn place_exhaustive(
  mat: &Mat,
  query: &Query,
  prefer: Option<Target>,
  seed: Option<Placement>,
  prune: bool,
) -> Result<Placement, Report> {
  let reference = mat.reference();
  let mut best = Best::new(prefer);
  if let Some(seed) = seed {
    best.offer(seed.target, signed(seed.cost));
    best.n_best = 0;
  }
  let bound = |node: GraphNodeKey, cost: isize| cost - signed(mat.depth_below(node));
  let mut stack = vec![(mat.root(), signed(query.root_cost(reference)))];
  let mut children = vec![];
  let mut visited = 0;
  while let Some((node, cost)) = stack.pop() {
    if prune && bound(node, cost) > best.cost {
      continue;
    }
    visited += 1;
    offer_node(mat, &mut best, node, cost);
    children.clear();
    for (child, edge) in mat.children(node) {
      let score = score_edge(mat.edge_subs(edge), query, reference);
      offer_edge(mat, &mut best, edge, cost, score);
      let child_cost = cost + score.delta;
      if !prune || bound(child, child_cost) <= best.cost {
        children.push((bound(child, child_cost), child, child_cost));
      }
    }
    children.sort_unstable_by_key(|&(bound, child, _)| (Reverse(bound), child));
    stack.extend(children.iter().map(|&(_, child, child_cost)| (child, child_cost)));
  }
  best.finish(visited)
}

pub fn place_greedy(
  mat: &Mat,
  query: &Query,
  starts: &[GraphNodeKey],
  slack: usize,
  prefer: Option<Target>,
) -> Result<Placement, Report> {
  let reference = mat.reference();
  let slack = signed(slack);
  let mut best = Best::new(prefer);
  let mut costs: BTreeMap<GraphNodeKey, isize> = BTreeMap::new();
  let mut heap = BinaryHeap::new();
  for &start in starts {
    if !mat.has_node(start) || costs.contains_key(&start) {
      continue;
    }
    let cost = signed(query.cost_against(&mat.node_state(start), reference));
    costs.insert(start, cost);
    heap.push(Reverse((cost, start)));
  }
  let mut visited = 0;
  while let Some(Reverse((cost, node))) = heap.pop() {
    if best.target.is_some() && cost > best.cost + slack {
      break;
    }
    visited += 1;
    offer_node(mat, &mut best, node, cost);
    for (child, edge) in mat.children(node) {
      let score = score_edge(mat.edge_subs(edge), query, reference);
      offer_edge(mat, &mut best, edge, cost, score);
      let child_cost = cost + score.delta;
      if !costs.contains_key(&child) && child_cost <= best.cost.saturating_add(slack) {
        costs.insert(child, child_cost);
        heap.push(Reverse((child_cost, child)));
      }
    }
    if let Some((parent, edge)) = mat.parent(node) {
      let score = score_edge(mat.edge_subs(edge), query, reference);
      let parent_cost = cost - score.delta;
      offer_edge(mat, &mut best, edge, parent_cost, score);
      if !costs.contains_key(&parent) && parent_cost <= best.cost.saturating_add(slack) {
        costs.insert(parent, parent_cost);
        heap.push(Reverse((parent_cost, parent)));
      }
    }
  }
  best.finish(visited)
}

pub fn place_from(mat: &Mat, query: &Query, start: GraphNodeKey, reaches: &Reaches) -> Result<Placement, Report> {
  let reference = mat.reference();
  let mut best = Best::new(None);
  let reach = |value: Option<usize>| signed(value.unwrap_or(usize::MAX));
  let start_cost = signed(query.cost_against(&mat.node_state(start), reference));
  let mut heap = BinaryHeap::from([Reverse((isize::MIN, start, start_cost, None))]);
  let mut visited = 0;
  while let Some(Reverse((bound, node, cost, came_from))) = heap.pop() {
    if bound > best.cost {
      break;
    }
    visited += 1;
    offer_node(mat, &mut best, node, cost);
    for (child, edge) in mat.children(node) {
      if Some(child) == came_from {
        continue;
      }
      let score = score_edge(mat.edge_subs(edge), query, reference);
      offer_edge(mat, &mut best, edge, cost, score);
      let child_cost = cost + score.delta;
      let child_bound = child_cost.saturating_sub(reach(reaches.below(child)));
      if child_bound <= best.cost {
        heap.push(Reverse((child_bound, child, child_cost, Some(node))));
      }
    }
    if let Some((parent, edge)) = mat.parent(node) {
      if Some(parent) != came_from {
        let score = score_edge(mat.edge_subs(edge), query, reference);
        let parent_cost = cost - score.delta;
        offer_edge(mat, &mut best, edge, parent_cost, score);
        let beyond = reaches
          .above(node)
          .map(|above| above.saturating_sub(mat.edge_subs(edge).len()));
        let parent_bound = parent_cost.saturating_sub(reach(beyond));
        if parent_bound <= best.cost {
          heap.push(Reverse((parent_bound, parent, parent_cost, Some(node))));
        }
      }
    }
  }
  best.finish(visited)
}
