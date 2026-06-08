use crate::clock::clock_graph::GraphClock;
use crate::clock::clock_model::{ClockLine, ClockModel};
use crossbeam_queue::ArrayQueue;
use crossbeam_skiplist::SkipMap;
use eyre::Report;
use itertools::Itertools;
use serde::{Deserialize, Serialize};
use std::path::Path;
use treetime_graph::breadth_first::GraphTraversalContinuation;
use treetime_graph::edge::HasBranchLength;
use treetime_io::csv::CsvStructFileWriter;
use treetime_utils::array::serde::skip_serializing_if_false;

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct ClockRegressionResult {
  pub name: Option<String>,
  pub div: f64,
  pub date: Option<f64>,
  pub predicted_date: f64,
  pub clock_deviation: Option<f64>,
  #[serde(serialize_with = "skip_serializing_if_false")]
  pub is_outlier: bool,
  #[serde(skip)]
  pub is_leaf: bool,
}

/// Get results of the root-to-tip clock inference.
pub fn gather_clock_regression_results(
  graph: &GraphClock,
  clock_model: &ClockModel,
) -> Result<Vec<ClockRegressionResult>, Report> {
  let result = ArrayQueue::new(graph.num_nodes());
  let divs = SkipMap::new();

  graph.par_iter_breadth_first_forward(|node| {
    let div = node.get_exactly_one_parent().map_or(0.0, |(parent, edge)| {
      let branch_length = edge.read_arc().branch_length().unwrap_or_default();
      let parent_div = parent.read_arc().div;
      parent_div + branch_length
    });

    let is_leaf = node.is_leaf;

    let mut node = node.payload;
    node.div = div;

    let name = node.name.clone();
    let time = node.time;
    let is_outlier = node.is_outlier;
    let predicted_date = clock_model.date(div);
    let clock_deviation = node.time.map(|time| clock_model.clock_deviation(time, div));

    divs.insert(name.clone(), div);

    result
      .push(ClockRegressionResult {
        name,
        div,
        date: time,
        predicted_date,
        clock_deviation,
        is_outlier,
        is_leaf,
      })
      .expect("ArrayQueue::push() failed. Queue is full.");

    Ok(GraphTraversalContinuation::Continue)
  })?;

  Ok(result.into_iter().collect_vec())
}

pub fn write_clock_regression_result_csv(
  results: &[ClockRegressionResult],
  filepath: impl AsRef<Path>,
  delimiter: u8,
) -> Result<(), Report> {
  let mut rtt_writer = CsvStructFileWriter::new(filepath, delimiter)?;
  results.iter().try_for_each(|result| rtt_writer.write(result))
}
