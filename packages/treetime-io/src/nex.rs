use crate::nwk::{NwkNodeComments, NwkWriteOptions, newick_graph, newick_write_options};
use eyre::{Report, WrapErr};
use std::collections::BTreeMap;
use std::io::Write;
use std::path::Path;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::node::GraphNodeKey;
use treetime_graph::tree_view::TreeView;
use treetime_utils::io::file::write_file_with;
use util_newick::{NexusTreeRef, NexusWriteOptions, nexus_to_writer};

pub fn nex_write_file(
  filepath: impl AsRef<Path>,
  tree: &TreeView<'_>,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  weights: &BTreeMap<GraphEdgeKey, Option<f64>>,
  options: &NwkWriteOptions,
  comments: &NwkNodeComments,
) -> Result<(), Report> {
  write_file_with(filepath, |writer| {
    nex_write(writer, tree, names, weights, options, comments)
  })
}

pub fn nex_write(
  mut writer: impl Write,
  tree: &TreeView<'_>,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  weights: &BTreeMap<GraphEdgeKey, Option<f64>>,
  options: &NwkWriteOptions,
  comments: &NwkNodeComments,
) -> Result<(), Report> {
  let graph = newick_graph(tree, names, weights, comments, options.style).wrap_err("When writing Nexus")?;
  let nexus_options = NexusWriteOptions {
    newick: newick_write_options(options),
    translate: false,
  };
  nexus_to_writer(&mut writer, &[NexusTreeRef::new("tree1", &graph)], &nexus_options)?;
  Ok(())
}
