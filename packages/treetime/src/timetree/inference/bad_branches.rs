use crate::clock::date_constraints::DateConstraints;
use eyre::Report;
use std::collections::{BTreeMap, BTreeSet};
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;

pub(crate) fn bad_leaves(
  graph: &Graph,
  constraints: &DateConstraints,
  outliers: &BTreeSet<GraphNodeKey>,
) -> BTreeMap<GraphNodeKey, bool> {
  graph
    .get_leaves()
    .map(|leaf| {
      let key = leaf.key();
      (
        key,
        constraints.date_constraint(key).is_none() || outliers.contains(&key),
      )
    })
    .collect()
}

pub(crate) fn derive_bad_branches(
  graph: &Graph,
  constraints: &DateConstraints,
  leaf_bad_branches: &BTreeMap<GraphNodeKey, bool>,
) -> Result<BTreeMap<GraphNodeKey, bool>, Report> {
  let mut bad_branches = BTreeMap::new();
  graph.iter_depth_first_postorder_forward(|node| {
    let bad = if node.is_leaf {
      leaf_bad_branches[&node.key]
    } else if constraints.date_constraint(node.key).is_some() {
      false
    } else {
      node.child_keys.iter().all(|(child_key, _)| bad_branches[child_key])
    };
    bad_branches.insert(node.key, bad);
    Ok(())
  })?;
  Ok(bad_branches)
}
