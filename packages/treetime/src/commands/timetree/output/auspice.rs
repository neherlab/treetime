use crate::payload::timetree::{EdgeTimetree, NodeTimetree};
use crate::timetree::confidence::NodeConfidenceInterval;
use crate::timetree::output::auspice::build_timetree_auspice;
use eyre::{Report, WrapErr};
use log::info;
use std::path::Path;
use treetime_graph::graph::Graph;
use treetime_utils::io::json::{JsonPretty, json_write_file};

pub fn write_auspice_json(
  graph: &Graph<NodeTimetree, EdgeTimetree, ()>,
  confidence_intervals: Option<&[NodeConfidenceInterval]>,
  outdir: &Path,
) -> Result<(), Report> {
  let tree = build_timetree_auspice(graph, confidence_intervals)?;
  let filepath = outdir.join("auspice_tree.json");
  json_write_file(&filepath, &tree, JsonPretty(true)).wrap_err("Failed to write auspice JSON")?;
  info!("Wrote auspice JSON to {filepath}", filepath = filepath.display());
  Ok(())
}
