use crate::clock::clock_model::ClockModel;
use crate::partition::timetree::GraphTimetree;
use crate::payload::traits::ClockNode;
use itertools::Itertools;
use log::warn;
use ordered_float::OrderedFloat;
use treetime_graph::node::{Named, Outlier};
use treetime_utils::fmt::string::truncate_right_with_ellipsis;

#[derive(Debug, Clone)]
pub struct OutlierRecord {
  pub name: String,
  pub given_date: f64,
  pub apparent_date: f64,
  pub residual: f64,
}

/// Collect outlier records from the graph for nodes marked as outliers.
pub fn collect_outliers(graph: &GraphTimetree, clock_model: &ClockModel, iqd: f64) -> Vec<OutlierRecord> {
  graph
    .get_leaves()
    .iter()
    .filter_map(|leaf| {
      let node = leaf.read_arc();
      let payload_arc = node.payload();
      let payload = payload_arc.read();
      if !payload.is_outlier() {
        return None;
      }
      let name = payload.name().map(|n| n.as_ref().to_owned())?;
      let given_date = payload.likely_time()?;
      let div = payload.div();
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
pub fn report_bad_branches(graph: &GraphTimetree, clock_model: &ClockModel, iqd: f64) {
  let outliers = collect_outliers(graph, clock_model, iqd);
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
/// After clock_filter_inplace marks leaves as outliers (is_outlier=true), this
/// sets bad_branch=true on those leaves and propagates upward: an internal node
/// is bad only when all its children are bad.
pub fn apply_outlier_bad_branches(graph: &GraphTimetree) {
  // Set bad_branch on outlier leaves
  for leaf in graph.get_leaves() {
    let node = leaf.read_arc();
    let mut payload = node.payload().write_arc();
    if payload.is_outlier() {
      payload.bad_branch = true;
    }
  }

  // Propagate upward in postorder: parent is bad only when all children are bad
  graph.iter_depth_first_postorder_forward(|mut node| {
    if node.is_leaf {
      return;
    }

    let all_children_bad = node.children.iter().all(|(child, _)| child.read_arc().bad_branch);

    node.payload.bad_branch = all_children_bad;
  });
}
