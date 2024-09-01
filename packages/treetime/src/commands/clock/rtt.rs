use crate::commands::clock::clock::ClockGraph;
use crate::commands::clock::clock_set::ClockModel;
use crate::io::csv::CsvStructFileWriter;
use eyre::Report;
use rayon::prelude::*;
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
}

/// Get results of the root-to-tip clock inference.
pub fn gather_clock_regression_results(graph: &ClockGraph, clock_model: &ClockModel) -> Vec<ClockRegressionResult> {
  graph
    .get_nodes()
    .par_iter()
    .map(|node| {
      let node = node.read_arc();
      let node = node.payload().read_arc();
      let name = node.name.clone();
      let date = node.date;
      let div = node.div;
      let is_outlier = node.is_outlier;
      let predicted_date = clock_model.date(div);
      let clock_deviation = node.date.map(|date| clock_model.clock_deviation(date, div));
      ClockRegressionResult {
        name,
        div,
        date,
        predicted_date,
        clock_deviation,
        is_outlier,
      }
    })
    .collect()
}

pub fn write_clock_regression_result_csv(
  results: &[ClockRegressionResult],
  filepath: impl AsRef<Path>,
  delimiter: u8,
) -> Result<(), Report> {
  let mut rtt_writer = CsvStructFileWriter::new(filepath, delimiter)?;
  results.iter().try_for_each(|result| rtt_writer.write(result))
}
