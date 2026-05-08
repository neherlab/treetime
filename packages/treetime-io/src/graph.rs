use crate::graphviz::{EdgeToGraphviz, NodeToGraphviz, graphviz_write_file};
use crate::nex::{NexWriteOptions, nex_write_file};
use crate::nwk::{EdgeToNwk, NodeToNwk, NwkWriteOptions, nwk_write_file};
use eyre::Report;
use serde::Serialize;
use std::path::Path;
use treetime_graph::edge::GraphEdge;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNode;
use treetime_utils::io::json::{JsonPretty, json_write_file};

pub fn write_graph_files<N, E, D>(outdir: impl AsRef<Path>, stem: &str, graph: &Graph<N, E, D>) -> Result<(), Report>
where
  N: GraphNode + NodeToNwk + NodeToGraphviz + Serialize,
  E: GraphEdge + EdgeToNwk + EdgeToGraphviz + Serialize,
  D: Send + Sync + Default + Serialize,
{
  let outdir = outdir.as_ref();

  nwk_write_file(outdir.join(format!("{stem}.nwk")), graph, &NwkWriteOptions::default())?;

  nex_write_file(outdir.join(format!("{stem}.nexus")), graph, &NexWriteOptions::default())?;

  json_write_file(outdir.join(format!("{stem}.graph.json")), graph, JsonPretty(true))?;

  graphviz_write_file(outdir.join(format!("{stem}.dot")), graph)?;

  Ok(())
}
