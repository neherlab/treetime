use crate::clock::clock_model::{ClockLine, ClockModel};
use crate::clock::clock_state::{ClockInputs, ClockState};
use eyre::Report;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::path::Path;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_graph::pass::GraphPassNodeOutput;
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
  graph: &Graph,
  inputs: &ClockInputs,
  state: &mut ClockState,
  clock_model: &ClockModel,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
) -> Result<Vec<ClockRegressionResult>, Report> {
  // Assign divergence to each node: div = parent.div + branch_length, parents before children.
  // `names` and `branch_lengths` are the post-reroot value maps threaded from the pipeline.
  state.map_forward(graph, |context| {
    let mut node = context.input.clone();
    let parent_message = if let Some((edge_key, edge)) = context.parent_edge {
      let parent = context.parent.expect("Non-root node must have a parent");
      let branch_length = branch_lengths[&edge_key].unwrap_or_default();
      node.div = parent.div + branch_length;
      Some(edge.clone())
    } else {
      node.div = 0.0;
      None
    };
    Ok(GraphPassNodeOutput { node, parent_message })
  })?;

  // One result per node, in node order.
  graph
    .get_nodes()
    .map(|node| {
      let is_leaf = node.is_leaf();
      let name = names[&node.key()].clone();
      let node_state = state.node(node.key());
      let div = node_state.div;
      let time = inputs.likely_time(node.key());
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
