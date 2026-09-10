use crate::nwk::{CommentProviders, NwkWriteOptions, nwk_write_with};
use eyre::Report;
use itertools::Itertools;
use smart_default::SmartDefault;
use std::collections::BTreeMap;
use std::io::Write;
use std::path::Path;
use treetime_graph::edge::{GraphEdge, GraphEdgeKey};
use treetime_graph::graph::Graph;
use treetime_graph::node::{GraphNode, GraphNodeKey};
use treetime_utils::io::file::create_file_or_stdout;
use util_newick::NwkStyle;

#[derive(Clone, SmartDefault)]
pub struct NexWriteOptions {
  /// Annotation style: Plain suppresses annotations, Beast/Nhx emit structured comments.
  #[default(NwkStyle::Plain)]
  pub style: NwkStyle,

  /// Format node weights keeping this many significant digits
  pub weight_significant_digits: Option<u8>,

  /// Format node weights keeping this many decimal digits
  pub weight_decimal_digits: Option<i8>,
}

pub fn nex_write_file<N, E, D>(
  filepath: impl AsRef<Path>,
  graph: &Graph<N, E, D>,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  weights: &BTreeMap<GraphEdgeKey, Option<f64>>,
  options: &NexWriteOptions,
) -> Result<(), Report>
where
  N: GraphNode,
  E: GraphEdge,
  D: Sync + Send,
{
  let mut f = create_file_or_stdout(filepath)?;
  nex_write(&mut f, graph, names, weights, options)?;
  writeln!(f)?;
  Ok(())
}

pub fn nex_write_file_with<N, E, D>(
  filepath: impl AsRef<Path>,
  graph: &Graph<N, E, D>,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  weights: &BTreeMap<GraphEdgeKey, Option<f64>>,
  options: &NexWriteOptions,
  providers: &CommentProviders,
) -> Result<(), Report>
where
  N: GraphNode,
  E: GraphEdge,
  D: Sync + Send,
{
  let mut f = create_file_or_stdout(filepath)?;
  nex_write_with(&mut f, graph, names, weights, options, providers)?;
  writeln!(f)?;
  Ok(())
}

pub fn nex_write_str<N, E, D>(
  graph: &Graph<N, E, D>,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  weights: &BTreeMap<GraphEdgeKey, Option<f64>>,
  options: &NexWriteOptions,
) -> Result<String, Report>
where
  N: GraphNode,
  E: GraphEdge,
  D: Sync + Send,
{
  let providers = CommentProviders::new();
  nex_write_str_with(graph, names, weights, options, &providers)
}

/// Return the Nexus representation of a graph, augmented by external node comment providers.
pub fn nex_write_str_with<N, E, D>(
  graph: &Graph<N, E, D>,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  weights: &BTreeMap<GraphEdgeKey, Option<f64>>,
  options: &NexWriteOptions,
  providers: &CommentProviders,
) -> Result<String, Report>
where
  N: GraphNode,
  E: GraphEdge,
  D: Sync + Send,
{
  let mut buf = Vec::new();
  nex_write_with(&mut buf, graph, names, weights, options, providers)?;
  Ok(String::from_utf8(buf)?)
}

pub fn nex_write<N, E, D>(
  w: &mut impl Write,
  graph: &Graph<N, E, D>,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  weights: &BTreeMap<GraphEdgeKey, Option<f64>>,
  options: &NexWriteOptions,
) -> Result<(), Report>
where
  N: GraphNode,
  E: GraphEdge,
  D: Sync + Send,
{
  let providers = CommentProviders::new();
  nex_write_with(w, graph, names, weights, options, &providers)
}

/// Write a graph in Nexus format, passing node names, edge weights, and comment providers through to
/// the embedded Newick tree.
pub fn nex_write_with<N, E, D>(
  w: &mut impl Write,
  graph: &Graph<N, E, D>,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  weights: &BTreeMap<GraphEdgeKey, Option<f64>>,
  options: &NexWriteOptions,
  providers: &CommentProviders,
) -> Result<(), Report>
where
  N: GraphNode,
  E: GraphEdge,
  D: Sync + Send,
{
  let n_leaves = graph.num_leaves();
  let leaf_names = graph
    .get_leaves()
    .iter()
    .filter_map(|n| names[&n.read_arc().key()].clone())
    .join(" ");
  let mut nwk = Vec::new();
  nwk_write_with(
    &mut nwk,
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
  let nwk = String::from_utf8(nwk)?;
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
  )?;

  Ok(())
}
