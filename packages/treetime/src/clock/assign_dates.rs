use crate::clock::clock_state::ClockInputs;
use crate::error::input_error;
use crate::make_error;
use eyre::Report;
use std::collections::BTreeMap;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_graph::schedule::GraphSchedule;
use treetime_primitives::date::DateConstraint;

const MIN_GOOD_LEAVES: usize = 3;

pub(crate) fn assign_dates(
  graph: &Graph,
  dates: &BTreeMap<GraphNodeKey, DateConstraint>,
  inputs: &mut ClockInputs,
) -> Result<(), Report> {
  if dates.is_empty() {
    return Err(input_error(
      "No valid date information found: no node of the tree has a usable date in the dates input",
    ));
  }

  let schedule = GraphSchedule::new(graph)?;
  let mut n_bad_leaves = 0;
  for &key in schedule.postorder() {
    let time: Option<f64> = dates.get(&key).map(DateConstraint::mean).filter(|&d| d.is_finite());

    let children = schedule.children(key);
    let is_leaf = children.is_empty();
    let bad_branch = time.is_none() && children.iter().all(|&(child_key, _)| inputs.node(child_key).bad_branch);

    let node_input = inputs.node_mut(key);
    node_input.time = time;
    node_input.bad_branch = bad_branch;

    if is_leaf && bad_branch {
      n_bad_leaves += 1;
    }
  }

  let n_leaves = graph.num_leaves();
  if n_leaves - n_bad_leaves < MIN_GOOD_LEAVES {
    return make_error!(
      "Not enough valid date constraints: there are {n_leaves} leaf nodes and {n_bad_leaves} of them have no date information"
    );
  }

  Ok(())
}
