use crate::clock::clock_model::{ClockLine, ClockModel};
use crate::clock::clock_regression::ClockRegressionPoint;
use crate::clock::clock_state::ClockInputs;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_utils::array::serde::{false_if_missing, skip_serializing_if_false};

pub(crate) fn gather_clock_regression_results(
  graph: &Graph,
  inputs: &ClockInputs,
  divergences: &BTreeMap<GraphNodeKey, f64>,
  outliers: &BTreeSet<GraphNodeKey>,
  clock_model: &ClockModel,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
) -> Vec<ClockRegressionResult> {
  graph
    .get_nodes()
    .map(|node| {
      let key = node.key();
      let div = divergences[&key];
      let time = inputs.likely_time(key);
      ClockRegressionResult {
        name: names[&key].clone(),
        div,
        date: time,
        predicted_date: clock_model.date(div),
        clock_deviation: time.map(|time| clock_model.clock_deviation(time, div)),
        is_outlier: outliers.contains(&key),
        is_leaf: node.is_leaf(),
        date_source: None,
      }
    })
    .collect()
}

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct ClockRegressionResult {
  pub name: Option<String>,
  pub div: f64,
  pub date: Option<f64>,
  pub predicted_date: f64,
  pub clock_deviation: Option<f64>,
  #[serde(
    serialize_with = "skip_serializing_if_false",
    deserialize_with = "false_if_missing",
    default
  )]
  pub is_outlier: bool,
  #[serde(skip)]
  pub is_leaf: bool,
  #[serde(default, skip_serializing_if = "Option::is_none")]
  pub date_source: Option<ClockDateSource>,
}

/// Where the date a clock regression used for a sample came from.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Deserialize, Serialize, JsonSchema)]
#[serde(rename_all = "kebab-case")]
pub enum ClockDateSource {
  /// The sampling date given in the input.
  Input,
  /// The date the time tree inferred for a sample without an input date.
  Inferred,
  /// No date; the sample did not enter the regression.
  Missing,
}

pub(crate) fn clock_fit_regression_results(
  model: &ClockModel,
  points: &[ClockRegressionPoint],
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  date_source: impl Fn(GraphNodeKey) -> ClockDateSource,
) -> Vec<ClockRegressionResult> {
  points
    .iter()
    .map(|point| ClockRegressionResult {
      name: names[&point.key].clone(),
      div: point.div,
      date: point.date,
      predicted_date: model.date(point.div),
      clock_deviation: point.date.map(|date| model.clock_deviation(date, point.div)),
      is_outlier: point.is_outlier,
      is_leaf: true,
      date_source: Some(point.date.map_or(ClockDateSource::Missing, |_| date_source(point.key))),
    })
    .collect()
}
