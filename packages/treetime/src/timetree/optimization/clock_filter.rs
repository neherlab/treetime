use crate::clock::clock_model::{ClockLine, ClockModel};
use crate::clock::clock_state::ClockState;
use crate::timetree::timetree_state::TimetreeState;
use eyre::Report;
use itertools::Itertools;
use log::warn;
use ordered_float::OrderedFloat;
use std::collections::BTreeMap;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_utils::fmt::string::truncate_right_with_ellipsis;

#[derive(Debug, Clone)]
pub struct OutlierRecord {
  pub name: String,
  pub given_date: f64,
  pub apparent_date: f64,
  pub residual: f64,
}

pub fn collect_outliers(
  graph: &Graph,
  clock_state: &ClockState,
  clock_model: &ClockModel,
  iqd: f64,
  given_dates: &BTreeMap<GraphNodeKey, Option<f64>>,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
) -> Vec<OutlierRecord> {
  graph
    .get_leaves()
    .filter_map(|leaf| {
      let node = leaf;
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

pub fn report_bad_branches(
  graph: &Graph,
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

pub fn apply_outlier_bad_branches(
  graph: &Graph,
  clock_state: &ClockState,
  state: &mut TimetreeState,
) -> Result<(), Report> {
  for leaf in graph.get_leaves() {
    let node = leaf;
    if clock_state.node(node.key()).is_outlier {
      state.node_mut(node.key()).bad_branch = true;
    }
  }

  propagate_bad_branches(graph, state)
}

pub fn propagate_bad_branches(graph: &Graph, state: &mut TimetreeState) -> Result<(), Report> {
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
