use crate::clock::clock_graph::GraphClock;
use crate::clock::clock_model::{ClockLine, ClockModel};
use crate::clock::clock_state::ClockState;
use eyre::Report;
use serde::{Deserialize, Serialize};
use std::path::Path;
use treetime_graph::pass::GraphPassNodeOutput;
use treetime_graph::value_maps::edge_branch_lengths;
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
  state: &mut ClockState,
  clock_model: &ClockModel,
) -> Result<Vec<ClockRegressionResult>, Report> {
  // Assign divergence to each node: div = parent.div + branch_length, parents before children.
  let branch_lengths = edge_branch_lengths(graph);
  state.map_forward(graph, |context| {
    let mut node = context.input;
    let parent_message = if let Some((edge_key, edge)) = context.parent_edge {
      let parent = context.parent.expect("Non-root node must have a parent");
      let branch_length = branch_lengths[&edge_key].unwrap_or_default();
      node.div = parent.div + branch_length;
      Some(edge)
    } else {
      node.div = 0.0;
      None
    };
    Ok(GraphPassNodeOutput { node, parent_message })
  })?;

  // One result per node, in node order.
  graph
    .get_nodes()
    .iter()
    .map(|node| {
      let node = node.read_arc();
      let is_leaf = node.is_leaf();
      let name = node.payload().read_arc().name.clone();
      let node_state = state.node(node.key());
      let div = node_state.div;
      let time = node_state.time;
      let predicted_date = clock_model.date(div);
      let clock_deviation = time.map(|time| clock_model.clock_deviation(time, div));
      Ok(ClockRegressionResult {
        name,
        div,
        date: time,
        predicted_date,
        clock_deviation,
        is_outlier: node_state.is_outlier,
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
