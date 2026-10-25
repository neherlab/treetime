use crate::clock::date_constraints::DateConstraints;
use eyre::Report;
use std::collections::BTreeSet;
use treetime_graph::graph::Graph;
use treetime_graph::key_map::NodeMap;
use treetime_graph::node::GraphNodeKey;
use treetime_graph::schedule::GraphSchedule;

pub(crate) fn bad_leaves(
  graph: &Graph,
  constraints: &DateConstraints,
  outliers: &BTreeSet<GraphNodeKey>,
) -> NodeMap<bool> {
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
  leaf_bad_branches: &NodeMap<bool>,
) -> Result<NodeMap<bool>, Report> {
  let schedule = GraphSchedule::new(graph)?;
  let mut bad_branches = NodeMap::new();
  for &key in schedule.postorder() {
    let children = schedule.children(key);
    let bad = if children.is_empty() {
      leaf_bad_branches[&key]
    } else if constraints.date_constraint(key).is_some() {
      false
    } else {
      children.iter().all(|(child_key, _)| bad_branches[child_key])
    };
    bad_branches.insert(key, bad);
  }
  Ok(bad_branches)
}
