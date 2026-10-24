use crate::build::mat::Mat;
use itertools::Itertools;
use ordered_float::OrderedFloat;
use std::cmp::Reverse;
use std::collections::{BTreeMap, BTreeSet};
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::node::GraphNodeKey;
use treetime_utils::cast::ToF64Exact;

#[derive(Clone, Debug, PartialEq)]
pub enum ReferenceSelection {
  Lbi {
    max_references: usize,
    tau: Option<f64>,
    pseudo_length: f64,
  },
  LongBranches {
    cuts: usize,
  },
}

pub fn select_references(mat: &Mat, selection: &ReferenceSelection) -> Vec<GraphNodeKey> {
  match selection {
    ReferenceSelection::Lbi {
      max_references,
      tau,
      pseudo_length,
    } => {
      let tau = tau.unwrap_or_else(|| default_lbi_tau(mat));
      lbi_peaks(mat, tau, *pseudo_length, *max_references)
    },
    ReferenceSelection::LongBranches { cuts } => long_branch_centers(mat, *cuts),
  }
}

pub fn default_lbi_tau(mat: &Mat) -> f64 {
  let radius = eccentricities(mat, &BTreeSet::new())
    .values()
    .copied()
    .min()
    .unwrap_or(0);
  radius.max(1).to_f64_exact() / mat.samples().len().to_f64_exact().ln().max(1.0)
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Reaches {
  down: BTreeMap<GraphNodeKey, usize>,
  up: BTreeMap<GraphNodeKey, usize>,
}

impl Reaches {
  pub fn new(mat: &Mat, blocked: &BTreeSet<GraphEdgeKey>) -> Self {
    let preorder = mat.preorder();
    let open_children = |node: GraphNodeKey| {
      mat
        .children(node)
        .filter(|(_, edge)| !blocked.contains(edge))
        .map(|(child, edge)| (child, mat.edge_subs(edge).len()))
    };
    let mut down: BTreeMap<GraphNodeKey, usize> = BTreeMap::new();
    for &node in preorder.iter().rev() {
      let depth = open_children(node)
        .map(|(child, length)| length + down[&child])
        .max()
        .unwrap_or(0);
      down.insert(node, depth);
    }
    let mut up: BTreeMap<GraphNodeKey, usize> = BTreeMap::new();
    up.insert(mat.root(), 0);
    for &node in &preorder {
      let branches = open_children(node)
        .map(|(child, length)| (length + down[&child], child))
        .sorted_by_key(|&(reach, child)| (Reverse(reach), child))
        .collect_vec();
      for (child, edge) in mat.children(node) {
        let reach = if blocked.contains(&edge) {
          0
        } else {
          let sibling = branches
            .iter()
            .find(|&&(_, other)| other != child)
            .map_or(0, |&(reach, _)| reach);
          mat.edge_subs(edge).len() + up[&node].max(sibling)
        };
        up.insert(child, reach);
      }
    }
    Self { down, up }
  }

  pub fn below(&self, node: GraphNodeKey) -> Option<usize> {
    self.down.get(&node).copied()
  }

  pub fn above(&self, node: GraphNodeKey) -> Option<usize> {
    self.up.get(&node).copied()
  }

  pub fn eccentricities(&self) -> BTreeMap<GraphNodeKey, usize> {
    self
      .down
      .iter()
      .map(|(&node, &down)| (node, down.max(self.up.get(&node).copied().unwrap_or(0))))
      .collect()
  }
}

pub fn eccentricities(mat: &Mat, blocked: &BTreeSet<GraphEdgeKey>) -> BTreeMap<GraphNodeKey, usize> {
  Reaches::new(mat, blocked).eccentricities()
}

fn long_branch_centers(mat: &Mat, cuts: usize) -> Vec<GraphNodeKey> {
  let cut: BTreeSet<GraphEdgeKey> = mat
    .preorder()
    .into_iter()
    .filter(|&node| !mat.is_leaf(node))
    .filter_map(|node| mat.parent(node).map(|(_, edge)| edge))
    .sorted_by_key(|&edge| (Reverse(mat.edge_subs(edge).len()), edge))
    .take(cuts)
    .collect();
  let eccentricity = eccentricities(mat, &cut);
  let mut component: BTreeMap<GraphNodeKey, GraphNodeKey> = BTreeMap::new();
  for node in mat.preorder() {
    let owner = match mat.parent(node) {
      Some((parent, edge)) if !cut.contains(&edge) => component[&parent],
      _ => node,
    };
    component.insert(node, owner);
  }
  component
    .iter()
    .filter(|&(&node, _)| !mat.is_leaf(node))
    .into_group_map_by(|&(_, &owner)| owner)
    .into_values()
    .filter_map(|members| {
      members
        .into_iter()
        .map(|(&node, _)| node)
        .min_by_key(|node| (eccentricity[node], *node))
    })
    .sorted()
    .collect()
}

pub fn local_branching_index(mat: &Mat, tau: f64, pseudo_length: f64) -> BTreeMap<GraphNodeKey, f64> {
  let preorder = mat.preorder();
  let decay = |node: GraphNodeKey| -> Option<(GraphNodeKey, f64, f64)> {
    let (parent, edge) = mat.parent(node)?;
    let length = mat.edge_subs(edge).len().to_f64_exact() + pseudo_length;
    let factor = (-length / tau).exp();
    Some((parent, factor, tau * (1.0 - factor)))
  };

  let mut below: BTreeMap<GraphNodeKey, f64> = preorder.iter().map(|&node| (node, 0.0)).collect();
  let mut sent_up: BTreeMap<GraphNodeKey, f64> = BTreeMap::new();
  for &node in preorder.iter().rev() {
    if let Some((parent, factor, own)) = decay(node) {
      let message = own + factor * below[&node];
      sent_up.insert(node, message);
      *below.entry(parent).or_default() += message;
    }
  }

  let mut above: BTreeMap<GraphNodeKey, f64> = BTreeMap::new();
  above.insert(mat.root(), 0.0);
  for &node in &preorder {
    if let Some((parent, factor, own)) = decay(node) {
      let outside = above[&parent] + below[&parent] - sent_up[&node];
      above.insert(node, own + factor * outside);
    }
  }
  preorder
    .iter()
    .map(|&node| (node, below[&node] + above[&node]))
    .collect()
}

pub fn lbi_peaks(mat: &Mat, tau: f64, pseudo_length: f64, max_peaks: usize) -> Vec<GraphNodeKey> {
  let lbi = local_branching_index(mat, tau, pseudo_length);
  let is_peak = |node: GraphNodeKey| {
    let value = lbi[&node];
    let parent_lower = mat.parent(node).is_none_or(|(parent, _)| lbi[&parent] < value);
    parent_lower && mat.children(node).all(|(child, _)| lbi[&child] < value)
  };
  lbi
    .keys()
    .copied()
    .filter(|&node| !mat.is_leaf(node) && is_peak(node))
    .sorted_by_key(|&node| (Reverse(OrderedFloat(lbi[&node])), node))
    .take(max_peaks)
    .collect()
}
