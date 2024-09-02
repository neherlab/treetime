use crate::commands::clock::clock_graph::ClockGraph;
use crate::commands::clock::clock_model::ClockModel;
use crate::graph::breadth_first::GraphTraversalContinuation;
use crate::io::csv::CsvStructFileWriter;
use crossbeam_queue::ArrayQueue;
use crossbeam_skiplist::SkipMap;
use eyre::Report;
use itertools::Itertools;
use ordered_float::OrderedFloat;
use serde::{Deserialize, Serialize};
use std::path::Path;

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct ClockRegressionResult {
  pub name: Option<String>,
  pub div: f64,
  pub date: Option<f64>,
  pub predicted_date: f64,
  pub clock_deviation: Option<f64>,
  pub is_outlier: bool,
  pub is_leaf: bool,
}

/// Get results of the root-to-tip clock inference.
pub fn gather_clock_regression_results(
  graph: &ClockGraph,
  clock_model: &ClockModel,
  clock_filter_threshold: f64,
) -> Vec<ClockRegressionResult> {
  let result = ArrayQueue::new(graph.num_nodes());
  let divs = SkipMap::new();

  graph.par_iter_breadth_first_forward(|node| {
    let div = node
      .get_exactly_one_parent()
      .map(|(parent, edge)| {
        let branch_length = edge.read_arc().branch_length.unwrap_or_default();
        let parent_div = parent.read_arc().div;
        parent_div + branch_length
      })
      .unwrap_or(0.0);

    let is_leaf = node.is_leaf;

    let mut node = node.payload;
    node.div = div;

    let name = node.name.clone();
    let date = node.date;

    let is_outlier = node.is_outlier;
    let predicted_date = clock_model.date(div);
    let clock_deviation = node.date.map(|date| clock_model.clock_deviation(date, div));

    divs.insert(name.clone(), div);

    result
      .push(ClockRegressionResult {
        name,
        div,
        date,
        predicted_date,
        clock_deviation,
        is_outlier,
        is_leaf,
      })
      .expect("ArrayQueue::push() failed. Queue is full.");

    GraphTraversalContinuation::Continue
  });

  // calculate the interquartile range of the clock_deviation of all leaf nodes
  // first, generated a sorted list of clock_deviation of all leaf nodes from the result array
  let leaf_clock_deviations: Vec<f64> = result
    .into_iter()
    .filter(|result| result.is_leaf)
    .filter_map(|result| result.clock_deviation.map(OrderedFloat))
    .sorted()
    .map(|cd| cd.into_inner())
    .collect();

  // calculate the interquartile range by taking the difference between the 3/4 and 1/4 quantile
  let iqd =
    leaf_clock_deviations[3 * leaf_clock_deviations.len() / 4] - leaf_clock_deviations[leaf_clock_deviations.len() / 4];

  // mark all leaf nodes as outliers if their absolute clock_deviation is greater than clock_filter_threshold * iqd
  result.into_iter().for_each(|mut result| {
    if result.is_leaf
      && result
        .clock_deviation
        .map(|deviation| deviation.abs() > clock_filter_threshold * iqd)
        .unwrap_or(false)
    {
      result.is_outlier = true;
    }
  });

  result.into_iter().collect_vec()
}

pub fn write_clock_regression_result_csv(
  results: &[ClockRegressionResult],
  filepath: impl AsRef<Path>,
  delimiter: u8,
) -> Result<(), Report> {
  let mut rtt_writer = CsvStructFileWriter::new(filepath, delimiter)?;
  results.iter().try_for_each(|result| rtt_writer.write(result))
}
