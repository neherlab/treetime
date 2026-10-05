use crate::nwk::{NwkNodeComments, NwkWriteOptions, nwk_write_str};
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
  let leaf_names = leaf_names(tree.graph(), names);
  let n_leaves = leaf_names.len();
  let leaf_names = tax_labels(&leaf_names)?;
  let nwk = nwk_write_str(tree, names, weights, options, comments)?;
  let nwk = nwk.strip_suffix(';').unwrap_or(&nwk);

  write!(
    writer,
    r#"#NEXUS
Begin Taxa;
  Dimensions NTax={n_leaves};
  TaxLabels {leaf_names};
End;
Begin Trees;
  Tree tree1={nwk};
End;
"#
  )
  .wrap_err("When writing Nexus")
}

fn leaf_names<'a>(graph: &Graph, names: &'a BTreeMap<GraphNodeKey, Option<String>>) -> Vec<&'a str> {
  graph
    .get_leaves()
    .filter_map(|leaf| names[&leaf.key()].as_deref())
    .collect()
}

fn tax_labels(leaf_names: &[&str]) -> Result<String, Report> {
  let mut labels = String::new();
  for name in leaf_names {
    if !labels.is_empty() {
      labels.push(' ');
    }
    write_label(&mut labels, name)?;
  }
  Ok(labels)
}
