use crate::clock::date_constraints::DateConstraints;
use crate::timetree::inference::runner::{EPS, GRID_POINTS};
use crate::timetree::timetree_state::{DateEdgeState, DateNodeState, TimetreeState};
use eyre::{Report, WrapErr};
use log::{Level, debug, log_enabled, warn};
use std::collections::BTreeMap;
use std::sync::Arc;
use treetime_distribution::Distribution;
use treetime_distribution::NegLog;
use treetime_distribution::convolve_across_edge;
use treetime_distribution::distribution_division;
use treetime_distribution::distribution_multiplication;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_graph::pass::{GraphPassForwardContext, GraphPassNodeOutput};
use treetime_grid::Side;

pub(crate) fn propagate_distributions_forward(
  graph: &Graph,
  constraints: &DateConstraints,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  state: &mut TimetreeState,
) -> Result<(), Report> {
  state.map_forward(graph, |context| {
    propagate_distributions_forward_node(constraints, names, &context)
  })?;

  let contradicted = state.nodes.values().filter(|node| node.contradicted).count();
  if contradicted > 0 {
    warn!(
      "Timetree forward pass: {contradicted} node(s) carry a date that the rest of the tree gives \
       no probability at all, so their posterior came out empty and each kept the date it was \
       given, unrefined. The usual cause is a sequence whose divergence implies a date far from \
       the one it is stamped with, which the clock filter reports separately. Run with \
       `--verbosity=debug` to see which nodes and where the tree puts each of them."
    );
  }

  Ok(())
}

fn propagate_distributions_forward_node(
  constraints: &DateConstraints,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  context: &GraphPassForwardContext<'_, DateNodeState, DateEdgeState, DateNodeState>,
) -> Result<GraphPassNodeOutput<DateNodeState, DateEdgeState>, Report> {
  let mut node = context.input.clone();
  let date_constraint = constraints.date_constraints.get(&context.key).cloned().flatten();
  let edge = context.parent_edge.map(|(_, edge)| edge);
  if refine_distribution_from_parent(
    names,
    context.key,
    date_constraint.as_ref(),
    context.parent,
    edge,
    &mut node,
  )? == Refinement::ContradictedGivenDate
  {
    node.contradicted = true;
  }
  commit_node_time(
    names,
    context.key,
    date_constraint.as_ref(),
    context.parent,
    context.is_leaf,
    &mut node,
  );
  let parent_message = context.parent_edge.map(|(_, edge)| edge.clone());
  Ok(GraphPassNodeOutput { node, parent_message })
}

#[allow(
  clippy::expect_used,
  reason = "expect on a value an upstream invariant guarantees is present"
)]
fn refine_distribution_from_parent(
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  key: GraphNodeKey,
  date_constraint: Option<&Arc<Distribution<NegLog>>>,
  parent: Option<&DateNodeState>,
  edge: Option<&DateEdgeState>,
  node: &mut DateNodeState,
) -> Result<Refinement, Report> {
  let Some(parent) = parent else {
    return Ok(Refinement::Done);
  };
  if has_exact_date(date_constraint) {
    return Ok(Refinement::Done);
  }

  let edge = edge.expect("Non-root indexed node must own its parent edge");

  let (Some(parent_time_dist), Some(branch_dist)) = (&parent.time_distribution, &edge.branch_length_distribution)
  else {
    return Ok(Refinement::Done);
  };

  let Some(subtree_dist) = &node.time_distribution else {
    let dist_from_parent = convolve_across_edge(parent_time_dist, branch_dist, Side::Right, EPS, GRID_POINTS)?;
    log_refinement(names, key, parent_time_dist, &dist_from_parent);
    node.time_distribution = Some(Arc::new(dist_from_parent));
    return Ok(Refinement::Done);
  };

  let parent_except_subtree = match edge.msg_to_parent.as_deref() {
    Some(msg_to_parent) => distribution_division(parent_time_dist, msg_to_parent)?,
    None => parent_time_dist.as_ref().clone(),
  };
  let dist_from_parent = convolve_across_edge(&parent_except_subtree, branch_dist, Side::Right, EPS, GRID_POINTS)?;

  let combined = distribution_multiplication(&dist_from_parent, subtree_dist)?
    .normalize()
    .wrap_err_with(|| format!("When normalizing the time distribution of node {key}"))?;
  log_refinement(names, key, parent_time_dist, &combined);

  if combined.likely_time().is_none() && date_constraint.is_some() {
    log_kept_given_date(names, key, date_constraint, &dist_from_parent);
    return Ok(Refinement::ContradictedGivenDate);
  }
  node.time_distribution = Some(Arc::new(combined));
  Ok(Refinement::Done)
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Refinement {
  Done,
  ContradictedGivenDate,
}

fn commit_node_time(
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  key: GraphNodeKey,
  date_constraint: Option<&Arc<Distribution<NegLog>>>,
  parent: Option<&DateNodeState>,
  is_leaf: bool,
  node: &mut DateNodeState,
) {
  let parent_time = (!has_exact_date(date_constraint))
    .then(|| parent_time(parent))
    .flatten();

  let is_dateable = !is_leaf || date_constraint.is_some();
  if set_likely_time(node, parent_time).is_none() && is_dateable {
    let name = node_name(names, key);
    let name = name.as_deref().unwrap_or("<unnamed>");
    warn!(
      "Timetree forward pass: node '{name}' has an empty time distribution; no date was assigned. \
       The messages meeting at this node leave no time with any probability: the dates below it \
       and the times the rest of the tree implies have disjoint support."
    );
  }
}

fn has_exact_date(date_constraint: Option<&Arc<Distribution<NegLog>>>) -> bool {
  date_constraint.is_some_and(|dist| dist.is_point())
}

fn parent_time(parent: Option<&DateNodeState>) -> Option<f64> {
  parent?.time
}

pub(super) fn set_likely_time(node: &mut DateNodeState, parent_time: Option<f64>) -> Option<f64> {
  let time = node
    .time_distribution
    .as_ref()
    .and_then(|time_dist| time_dist.likely_time())?;

  let time = parent_time.map_or(time, |parent_time| time.max(parent_time));
  node.time = Some(time);
  Some(time)
}

fn log_refinement(
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  key: GraphNodeKey,
  parent: &Distribution<NegLog>,
  refined: &Distribution<NegLog>,
) {
  if !log_enabled!(Level::Debug) {
    return;
  }
  let name = node_name(names, key);
  let name = name.as_deref().unwrap_or("<unnamed>");
  debug!(
    "Timetree forward pass: node '{name}': parent {} -> refined {}",
    describe_grid(parent),
    describe_grid(refined)
  );
}

fn log_kept_given_date(
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  key: GraphNodeKey,
  date_constraint: Option<&Arc<Distribution<NegLog>>>,
  dist_from_parent: &Distribution<NegLog>,
) {
  if !log_enabled!(Level::Debug) {
    return;
  }
  let name = node_name(names, key);
  let name = name.as_deref().unwrap_or("<unnamed>");
  let given = date_constraint.map_or_else(|| "none".to_owned(), |constraint| describe_grid(constraint.as_ref()));
  debug!(
    "Timetree forward pass: node '{name}' keeps the date it was given, {given}: the rest of the \
     tree puts it at {}, which leaves no probability on that date",
    describe_grid(dist_from_parent)
  );
}

fn node_name(names: &BTreeMap<GraphNodeKey, Option<String>>, key: GraphNodeKey) -> Option<String> {
  names[&key].clone()
}

fn describe_grid(dist: &Distribution<NegLog>) -> String {
  let n_points = dist.t().len();
  match dist.time_bounds() {
    Some((t_min, t_max)) => format!("n={n_points}, [{t_min:.6}, {t_max:.6}]"),
    None => "empty".to_owned(),
  }
}
