use crate::clock::clock_model::{ClockLine, ClockModel};
use crate::progress::ProgressSink;
use crate::progress_warn;
use itertools::Itertools;
use ordered_float::OrderedFloat;
use std::collections::{BTreeMap, BTreeSet};
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_utils::fmt::string::truncate_right_with_ellipsis;

pub(crate) fn report_bad_branches(
  graph: &Graph,
  outliers: &BTreeSet<GraphNodeKey>,
  divergences: &BTreeMap<GraphNodeKey, f64>,
  clock_model: &ClockModel,
  iqd: f64,
  given_dates: &BTreeMap<GraphNodeKey, Option<f64>>,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  progress: &dyn ProgressSink,
) {
  let records = collect_outlier_records(graph, outliers, divergences, clock_model, iqd, given_dates, names);
  if records.is_empty() {
    return;
  }

  progress_warn!(progress, "Clock filter marked {} outliers:", records.len());
  progress_warn!(
    progress,
    "{:>20} {:>12} {:>14} {:>10}",
    "name",
    "given_date",
    "apparent_date",
    "residual"
  );
  for r in &records {
    progress_warn!(
      progress,
      "{:>20} {:>12.2} {:>14.2} {:>10.2}",
      truncate_right_with_ellipsis(&r.name, 20),
      r.given_date,
      r.apparent_date,
      r.residual
    );
  }
}

fn collect_outlier_records(
  graph: &Graph,
  outliers: &BTreeSet<GraphNodeKey>,
  divergences: &BTreeMap<GraphNodeKey, f64>,
  clock_model: &ClockModel,
  iqd: f64,
  given_dates: &BTreeMap<GraphNodeKey, Option<f64>>,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
) -> Vec<OutlierRecord> {
  graph
    .get_leaves()
    .filter_map(|leaf| {
      let key = leaf.key();
      if !outliers.contains(&key) {
        return None;
      }
      let name = names[&key].clone()?;
      let given_date = given_dates.get(&key).copied().flatten()?;
      let div = divergences[&key];
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

struct OutlierRecord {
  name: String,
  given_date: f64,
  apparent_date: f64,
  residual: f64,
}

pub(crate) fn mark_outlier_leaves(
  graph: &Graph,
  outliers: &BTreeSet<GraphNodeKey>,
  leaf_bad_branches: &BTreeMap<GraphNodeKey, bool>,
) -> BTreeMap<GraphNodeKey, bool> {
  graph
    .get_leaves()
    .map(|leaf| {
      let key = leaf.key();
      (key, leaf_bad_branches[&key] || outliers.contains(&key))
    })
    .collect()
}
