use crate::clock::clock_state::ClockInputs;
use crate::make_error;
use eyre::Report;
use std::collections::BTreeMap;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_io::dates_csv::DatesMap;

const MIN_GOOD_LEAVES: usize = 3;

/// Assign each node's observed date and bad-branch flag into `state`.
///
/// `bad_branch` is set bottom-up: a node is bad when it has no date and every child is bad (or it
/// is a dateless leaf). The postorder walk visits children before parents, so each child's flag is
/// already in `state` when the parent reads it.
pub fn assign_dates(
  graph: &Graph,
  dates: &DatesMap,
  inputs: &mut ClockInputs,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
) -> Result<(), Report> {
  let n_dates = dates.iter().filter(|(_, d)| d.is_some()).count();
  if n_dates == 0 {
    return make_error!("No valid date information found in {dates:#?}");
  }

  let mut n_bad_leaves = 0;
  graph.iter_depth_first_postorder_forward(|node| {
    let name = names[&node.key].clone();
    let time: Option<f64> = name
      .and_then(|name| dates.get(&name))
      .and_then(|d| d.as_ref().map(|c| c.mean()))
      .filter(|&d| d.is_finite());

    let bad_branch = time.is_none()
      && (node.is_leaf
        || node
          .child_keys
          .iter()
          .all(|(child_key, _)| inputs.node(*child_key).bad_branch));

    let node_input = inputs.node_mut(node.key);
    node_input.time = time;
    node_input.bad_branch = bad_branch;

    if node.is_leaf && bad_branch {
      n_bad_leaves += 1;
    }
    Ok(())
  })?;

  let n_leaves = graph.num_leaves();
  if n_leaves - n_bad_leaves < MIN_GOOD_LEAVES {
    return make_error!(
      "Not enough valid date constraints: there are {n_leaves} leaf nodes and {n_bad_leaves} of them have no date information"
    );
  }

  Ok(())
}
