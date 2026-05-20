use crate::graphviz::{EdgeToGraphviz, NodeToGraphviz, graphviz_write_file};
use crate::nex::{NexWriteOptions, nex_write_file_with};
use crate::nwk::{CommentProviders, EdgeToNwk, NodeToNwk, NwkWriteOptions, nwk_write_file_with};
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
  write_graph_files_with(outdir, stem, graph, &CommentProviders::new())
}

pub fn write_graph_files_with<N, E, D>(
  outdir: impl AsRef<Path>,
  stem: &str,
  graph: &Graph<N, E, D>,
  providers: &CommentProviders,
) -> Result<(), Report>
where
  N: GraphNode + NodeToNwk + NodeToGraphviz + Serialize,
  E: GraphEdge + EdgeToNwk + EdgeToGraphviz + Serialize,
  D: Send + Sync + Default + Serialize,
{
  let outdir = outdir.as_ref();

  nwk_write_file_with(
    outdir.join(format!("{stem}.nwk")),
    graph,
    &NwkWriteOptions::default(),
    providers,
  )?;

  nex_write_file_with(
    outdir.join(format!("{stem}.nexus")),
    graph,
    &NexWriteOptions::default(),
    providers,
  )?;

  json_write_file(outdir.join(format!("{stem}.graph.json")), graph, JsonPretty(true))?;

  graphviz_write_file(outdir.join(format!("{stem}.dot")), graph)?;

  Ok(())
}
