use crate::clock::clock_model::{ClockLine, ClockModel};
use crate::clock::clock_state::ClockState;
use crate::partition::timetree::partition::GraphTimetree;
use crate::timetree::timetree_state::TimetreeState;
use eyre::Report;
use itertools::Itertools;
use log::warn;
use ordered_float::OrderedFloat;
use std::collections::BTreeMap;
use treetime_graph::node::GraphNodeKey;
use treetime_utils::fmt::string::truncate_right_with_ellipsis;

#[derive(Debug, Clone)]
pub struct OutlierRecord {
  pub name: String,
  pub given_date: f64,
  pub apparent_date: f64,
  pub residual: f64,
}

/// Collect outlier records from the clock state for leaves marked as outliers.
///
/// The outlier flag and divergence come from the threaded [`ClockState`] value, and the given date
/// from the date-state `given_dates` map (the same node dates the filter regressed on); the leaf name
/// stays transitional on the payload.
pub fn collect_outliers(
  graph: &GraphTimetree,
  clock_state: &ClockState,
  clock_model: &ClockModel,
  iqd: f64,
  given_dates: &BTreeMap<GraphNodeKey, Option<f64>>,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
) -> Vec<OutlierRecord> {
  graph
    .get_leaves()
    .iter()
    .filter_map(|leaf| {
      let node = leaf.read_arc();
      let state = clock_state.node(node.key());
      if !state.is_outlier {
        return None;
      }
      let name = names[&node.key()].clone()?;
      let given_date = given_dates.get(&node.key()).copied().flatten()?;
      let div = state.div;
      let apparent_date = clock_model.date(div);
      let clock_deviation = clock_model.clock_deviation(given_date, div);
      let residual = if iqd > 0.0 { clock_deviation / iqd } else { 0.0 };
      Some(OutlierRecord {
        name,
        given_date,
        apparent_date,
        residual,
      })
    })
    .sorted_by_key(|r| (OrderedFloat(r.residual.abs()), r.name.clone()))
    .rev()
    .collect_vec()
}

/// Report outlier branches that violate molecular clock.
pub fn report_bad_branches(
  graph: &GraphTimetree,
  clock_state: &ClockState,
  clock_model: &ClockModel,
  iqd: f64,
  given_dates: &BTreeMap<GraphNodeKey, Option<f64>>,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
) {
  let outliers = collect_outliers(graph, clock_state, clock_model, iqd, given_dates, names);
  if outliers.is_empty() {
    return;
  }

  warn!("Clock filter marked {} outliers:", outliers.len());
  warn!(
    "{:>20} {:>12} {:>14} {:>10}",
    "name", "given_date", "apparent_date", "residual"
  );
  for r in &outliers {
    warn!(
      "{:>20} {:>12.2} {:>14.2} {:>10.2}",
      truncate_right_with_ellipsis(&r.name, 20),
      r.given_date,
      r.apparent_date,
      r.residual
    );
  }
}

/// Convert outlier flags to bad_branch flags for backward pass exclusion.
///
/// After clock_filter_inplace marks leaves as outliers (is_outlier=true in the clock state), this
/// sets bad_branch=true on those leaves and propagates upward: an internal node
/// is bad only when all its children are bad.
///
/// The outlier flag is read from the threaded [`ClockState`] value; the bad-branch flag is written
/// into the threaded [`TimetreeState`] value, the home the coalescent and date passes read.
pub fn apply_outlier_bad_branches(
  graph: &GraphTimetree,
  clock_state: &ClockState,
  state: &mut TimetreeState,
) -> Result<(), Report> {
  for leaf in graph.get_leaves() {
    let node = leaf.read_arc();
    if clock_state.node(node.key()).is_outlier {
      state.node_mut(node.key()).bad_branch = true;
    }
  }

  propagate_bad_branches(graph, state)
}

/// Recompute internal bad-branch state from the current topology.
///
/// Each internal node's flag is the conjunction of its children's flags, read from the threaded
/// [`TimetreeState`] value and written back into it.
pub fn propagate_bad_branches(graph: &GraphTimetree, state: &mut TimetreeState) -> Result<(), Report> {
  graph.iter_depth_first_postorder_forward(|node| {
    if node.is_leaf {
      return Ok(());
    }

    let all_children_bad = node
      .child_keys
      .iter()
      .all(|(child_key, _)| state.node(*child_key).bad_branch);

    state.node_mut(node.key).bad_branch = all_children_bad;
    Ok(())
  })
}
