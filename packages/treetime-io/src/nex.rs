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
use treetime_utils::io::file::write_file_with;
use util_newick::NwkStyle;

pub fn nex_write_file_with(
  filepath: impl AsRef<Path>,
  graph: &Graph,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  weights: &BTreeMap<GraphEdgeKey, Option<f64>>,
  options: &NexWriteOptions,
  providers: &CommentProviders,
) -> Result<(), Report> {
  let filepath = filepath.as_ref();
  write_file_with(filepath, |f| {
    nex_write_with(&mut *f, graph, names, weights, options, providers)?;
    writeln!(f).map_err(Report::new)
  })
  .wrap_err_with(|| format!("When writing Nexus file '{}'", filepath.display()))
}

#[cfg_attr(
  dylint_lib = "treetime_lints",
  allow(
    pub_unused_in_workspace,
    reason = "used only by tests of other workspace crates, which a cfg(test) item cannot reach"
  )
)]
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

#[derive(Clone, SmartDefault)]
pub struct NexWriteOptions {
  #[default(NwkStyle::Plain)]
  pub style: NwkStyle,

  pub weight_significant_digits: Option<u8>,

  pub weight_decimal_digits: Option<i8>,
}
