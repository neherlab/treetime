use crate::nwk::{CommentProviders, NwkWriteOptions, nwk_write_str_with};
use eyre::{Report, WrapErr};
use itertools::Itertools;
use smart_default::SmartDefault;
use std::collections::BTreeMap;
use std::io::Write;
use std::path::Path;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_utils::io::file::create_file_or_stdout;
use util_newick::NwkStyle;

#[derive(Clone, SmartDefault)]
pub struct NexWriteOptions {
  #[default(NwkStyle::Plain)]
  pub style: NwkStyle,

  pub weight_significant_digits: Option<u8>,

  pub weight_decimal_digits: Option<i8>,
}

pub fn nex_write_file(
  filepath: impl AsRef<Path>,
  graph: &Graph,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  weights: &BTreeMap<GraphEdgeKey, Option<f64>>,
  options: &NexWriteOptions,
) -> Result<(), Report> {
  let filepath = filepath.as_ref();
  let context = || format!("When writing Nexus file '{}'", filepath.display());
  let mut f = create_file_or_stdout(filepath)?;
  nex_write(&mut f, graph, names, weights, options).wrap_err_with(context)?;
  writeln!(f).wrap_err_with(context)
}

pub fn nex_write_file_with(
  filepath: impl AsRef<Path>,
  graph: &Graph,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  weights: &BTreeMap<GraphEdgeKey, Option<f64>>,
  options: &NexWriteOptions,
  providers: &CommentProviders,
) -> Result<(), Report> {
  let filepath = filepath.as_ref();
  let context = || format!("When writing Nexus file '{}'", filepath.display());
  let mut f = create_file_or_stdout(filepath)?;
  nex_write_with(&mut f, graph, names, weights, options, providers).wrap_err_with(context)?;
  writeln!(f).wrap_err_with(context)
}

pub fn nex_write_str(
  graph: &Graph,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  weights: &BTreeMap<GraphEdgeKey, Option<f64>>,
  options: &NexWriteOptions,
) -> Result<String, Report> {
  let providers = CommentProviders::new();
  nex_write_str_with(graph, names, weights, options, &providers)
}

pub fn nex_write_str_with(
  graph: &Graph,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  weights: &BTreeMap<GraphEdgeKey, Option<f64>>,
  options: &NexWriteOptions,
  providers: &CommentProviders,
) -> Result<String, Report> {
  let mut buf = Vec::new();
  nex_write_with(&mut buf, graph, names, weights, options, providers)?;
  Ok(String::from_utf8(buf)?)
}

pub fn nex_write(
  w: &mut impl Write,
  graph: &Graph,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  weights: &BTreeMap<GraphEdgeKey, Option<f64>>,
  options: &NexWriteOptions,
) -> Result<(), Report> {
  let providers = CommentProviders::new();
  nex_write_with(w, graph, names, weights, options, &providers)
}

fn nex_write_with(
  w: &mut impl Write,
  graph: &Graph,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  weights: &BTreeMap<GraphEdgeKey, Option<f64>>,
  options: &NexWriteOptions,
  providers: &CommentProviders,
) -> Result<(), Report> {
  let n_leaves = graph.num_leaves();
  let leaf_names = graph.get_leaves().filter_map(|n| names[&n.key()].clone()).join(" ");
  let nwk = nwk_write_str_with(
    graph,
    names,
    weights,
    &NwkWriteOptions {
      style: options.style,
      weight_significant_digits: options.weight_significant_digits,
      weight_decimal_digits: options.weight_decimal_digits,
    },
    providers,
  )?;
  let nwk = nwk.strip_suffix(';').unwrap_or(&nwk);

  writeln!(
    w,
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
