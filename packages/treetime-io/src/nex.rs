use crate::nwk::{NwkNodeComments, NwkWriteOptions, write_nwk_tree};
use eyre::{Report, WrapErr};
use std::collections::BTreeMap;
use std::io::Write;
use std::path::Path;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_graph::tree_view::TreeView;
use treetime_utils::io::file::write_file_with;
use util_newick::write_label;

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
  write_nex(&mut writer, tree, names, weights, options, comments).wrap_err("When writing Nexus")
}

fn write_nex(
  writer: &mut impl Write,
  tree: &TreeView<'_>,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  weights: &BTreeMap<GraphEdgeKey, Option<f64>>,
  options: &NwkWriteOptions,
  comments: &NwkNodeComments,
) -> Result<(), Report> {
  let leaf_names = leaf_names(tree.graph(), names);
  writeln!(writer, "#NEXUS")?;
  writeln!(writer, "Begin Taxa;")?;
  writeln!(writer, "  Dimensions NTax={};", leaf_names.len())?;
  write!(writer, "  TaxLabels ")?;
  for (i, name) in leaf_names.iter().enumerate() {
    if i > 0 {
      write!(writer, " ")?;
    }
    write_label(writer, name)?;
  }
  writeln!(writer, ";")?;
  writeln!(writer, "End;")?;
  writeln!(writer, "Begin Trees;")?;
  write!(writer, "  Tree tree1=")?;
  write_nwk_tree(writer, tree, names, weights, options, comments)?;
  writeln!(writer, ";")?;
  writeln!(writer, "End;")?;
  Ok(())
}

fn leaf_names<'a>(graph: &Graph, names: &'a BTreeMap<GraphNodeKey, Option<String>>) -> Vec<&'a str> {
  graph
    .get_leaves()
    .filter_map(|leaf| names[&leaf.key()].as_deref())
    .collect()
}
