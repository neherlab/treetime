use crate::clock::clock_graph::{EdgeClock, GraphClock, NodeClock};
use crate::clock::clock_model::{ClockLine, ClockModel};
use eyre::Report;
use serde::{Deserialize, Serialize};
use std::path::Path;
use treetime_graph::edge::HasBranchLength;
use treetime_graph::pass::{GraphPassNodeOutput, with_graph_payloads_map};
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
  // Assign divergence to each node: div = parent.div + branch_length, parents before children.
  with_graph_payloads_map(graph, |pass| {
    pass.try_map_forward::<NodeClock, EdgeClock>(|context| {
      let mut node = context.input;
      let parent_message = if let Some((_, edge)) = context.parent_edge {
        let parent = context.parent.expect("Non-root node must have a parent");
        let div = parent.div + edge.branch_length().unwrap_or_default();
        node.div = div;
        Some(edge)
      } else {
        node.div = 0.0;
        None
      };
      Ok(GraphPassNodeOutput { node, parent_message })
    })
  })?;

  // One result per node, in node order.
  graph
    .get_nodes()
    .iter()
    .map(|node| {
      let node = node.read_arc();
      let is_leaf = node.is_leaf();
      let payload = node.payload();
      let payload = payload.read();
      let div = payload.div;
      let predicted_date = clock_model.date(div);
      let clock_deviation = payload.time.map(|time| clock_model.clock_deviation(time, div));
      Ok(ClockRegressionResult {
        name: payload.name.clone(),
        div,
        date: payload.time,
        predicted_date,
        clock_deviation,
        is_outlier: payload.is_outlier,
        is_leaf,
      })
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
