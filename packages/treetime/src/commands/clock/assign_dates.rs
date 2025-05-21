use crate::commands::clock::clock_graph::ClockGraph;
use crate::graph::node::Named;
use crate::io::dates_csv::{DateOrRange, DatesMap};
use crate::make_error;
use eyre::Report;

const MIN_GOOD_LEAVES: usize = 3;

pub fn assign_dates(graph: &ClockGraph, dates: &DatesMap) -> Result<(), Report> {
  let n_dates = dates.iter().filter(|(_, d)| d.is_some()).count();
  if n_dates == 0 {
    return make_error!("No valid date information found in {dates:#?}");
  }

  let mut n_bad_leaves = 0;
  graph.iter_depth_first_postorder_forward(|mut node| {
    let name = node.payload.name().map(|s| s.as_ref().to_owned());
    let date: Option<f64> = name
      .and_then(|name| dates.get(&name))
      .and_then(|d| d.as_ref().map(DateOrRange::mean))
      .filter(|&d| d.is_finite());

    node.payload.date = date;

    node.payload.bad_branch =
      date.is_none() && (node.is_leaf || node.children.iter().all(|(child, edge)| child.read_arc().bad_branch));

    if node.is_leaf && node.payload.bad_branch {
      n_bad_leaves += 1;
    }
  });

  // FIXME: this fails if n_leaves < MIN_GOOD_LEAVES
  let n_leaves = graph.num_leaves();
  if n_leaves - n_bad_leaves < MIN_GOOD_LEAVES {
    return make_error!(
      "Not enough valid date constraints: there are {n_leaves} leaf nodes and {n_bad_leaves} of them have no date information"
    );
  }

  Ok(())
}
