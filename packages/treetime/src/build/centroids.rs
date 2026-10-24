use crate::build::mat::Mat;
use itertools::Itertools;
use ordered_float::OrderedFloat;
use std::cmp::Reverse;
use std::collections::BTreeMap;
use treetime_graph::node::GraphNodeKey;
use treetime_utils::cast::ToF64Exact;

pub fn local_branching_index(mat: &Mat, tau: f64, pseudo_length: f64) -> BTreeMap<GraphNodeKey, f64> {
  let mut preorder = vec![];
  let mut stack = vec![mat.root()];
  while let Some(node) = stack.pop() {
    preorder.push(node);
    stack.extend(mat.children(node).map(|(child, _)| child));
  }
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
