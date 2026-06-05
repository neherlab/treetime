use crate::graphviz::{EdgeToGraphviz, NodeToGraphviz, graphviz_write_file};
use crate::nex::{NexWriteOptions, nex_write_file_with};
use crate::nwk::{CommentProviders, EdgeToNwk, NodeToNwk, NwkWriteOptions, nwk_write_file_with};
use eyre::Report;
use serde::Serialize;
use std::path::Path;
use treetime_graph::edge::{GraphEdge, HasBranchLength};
use treetime_graph::graph::Graph;
use treetime_graph::node::{GraphNode, Named};
use treetime_graph::topology_order::TopologyOrderSpec;
use treetime_utils::io::json::{JsonPretty, json_write_file};

#[derive(Clone, Debug, Default)]
pub struct GraphWriteOptions {
  pub topology_order: TopologyOrderSpec,
}

pub fn write_graph_files_with_options<N, E, D>(
  outdir: impl AsRef<Path>,
  stem: &str,
  graph: &Graph<N, E, D>,
  providers: &CommentProviders,
  options: &GraphWriteOptions,
) -> Result<(), Report>
where
  N: GraphNode + Named + NodeToNwk + NodeToGraphviz + Serialize,
  E: GraphEdge + HasBranchLength + EdgeToNwk + EdgeToGraphviz + Serialize,
  D: Send + Sync + Default + Serialize,
{
  let outdir = outdir.as_ref();
  let plan = options.topology_order.plan(graph)?;
  let graph = plan.ordered_graph(graph)?;

  nwk_write_file_with(
    outdir.join(format!("{stem}.nwk")),
    &graph,
    &NwkWriteOptions::default(),
    providers,
  )?;

  nex_write_file_with(
    outdir.join(format!("{stem}.nexus")),
    &graph,
    &NexWriteOptions::default(),
    providers,
  )?;

  json_write_file(outdir.join(format!("{stem}.graph.json")), &graph, JsonPretty(true))?;

  graphviz_write_file(outdir.join(format!("{stem}.dot")), &graph)?;

  Ok(())
}
