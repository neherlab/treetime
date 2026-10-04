use crate::nwk::{CommentProviders, NwkWriteOptions, nwk_write_str};
use eyre::{Report, WrapErr};
use std::collections::BTreeMap;
use std::io::Write;
use std::path::Path;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_utils::io::file::write_file_with;
use util_newick::write_label;

pub fn nex_write_file(
  filepath: impl AsRef<Path>,
  graph: &Graph,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  weights: &BTreeMap<GraphEdgeKey, Option<f64>>,
  options: &NwkWriteOptions,
  providers: &CommentProviders,
) -> Result<(), Report> {
  write_file_with(filepath, |writer| {
    nex_write(writer, graph, names, weights, options, providers)
  })
}

pub fn nex_write(
  mut writer: impl Write,
  graph: &Graph,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  weights: &BTreeMap<GraphEdgeKey, Option<f64>>,
  options: &NwkWriteOptions,
  providers: &CommentProviders,
) -> Result<(), Report> {
  let n_leaves = graph.num_leaves();
  let leaf_names = tax_labels(graph, names)?;
  let nwk = nwk_write_str(graph, names, weights, options, providers)?;
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

fn tax_labels(graph: &Graph, names: &BTreeMap<GraphNodeKey, Option<String>>) -> Result<String, Report> {
  let mut labels = String::new();
  for name in graph.get_leaves().filter_map(|leaf| names[&leaf.key()].as_deref()) {
    if !labels.is_empty() {
      labels.push(' ');
    }
    write_label(&mut labels, name)?;
  }
  Ok(labels)
}
