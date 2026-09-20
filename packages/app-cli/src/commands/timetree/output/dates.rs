use eyre::Report;
use std::path::Path;
use treetime_graph::graph::Graph;

pub fn write_node_dates(_graph: &Graph, _out_base: &Path) -> Result<(), Report> {
  todo!("Write node dates to TSV file")
}
