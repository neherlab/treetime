use eyre::{Report, WrapErr};
use log::warn;
use maplit::btreemap;
use smart_default::SmartDefault;
use std::collections::BTreeMap;
use std::io::{Read, Write};
use std::path::Path;
use std::sync::Arc;
use treetime_graph::assign_node_names::assign_node_names;
use treetime_graph::edge::GraphEdge;
use treetime_graph::graph::{Graph, SafeEdge, SafeNode};
use treetime_graph::node::GraphNodeKey;
use treetime_graph::node::{GraphNode, Named};
use treetime_utils::fmt::float::float_to_digits;
use treetime_utils::io::file::create_file_or_stdout;
use treetime_utils::io::file::open_file_or_stdin;
use treetime_utils::make_error;
use treetime_utils::make_report;
pub use util_newick::NwkStyle;
use util_newick::{newick_from_reader, newick_from_string, NewickGraph, NewickValue, write_beast_attrs, write_label, write_nhx_attrs};

pub fn nwk_read_file<N, E, D>(filepath: impl AsRef<Path>) -> Result<Graph<N, E, D>, Report>
where
  N: GraphNode + NodeFromNwk + Named,
  E: GraphEdge + EdgeFromNwk,
  D: Sync + Send + Default,
{
  let filepath = filepath.as_ref();
  nwk_read(open_file_or_stdin(&Some(filepath))?).wrap_err_with(|| format!("When reading file '{}'", filepath.display()))
}

pub fn nwk_read_str<N, E, D>(nwk_string: impl AsRef<str>) -> Result<Graph<N, E, D>, Report>
where
  N: GraphNode + NodeFromNwk + Named,
  E: GraphEdge + EdgeFromNwk,
  D: Sync + Send + Default,
{
  let nwk_graph = newick_from_string(nwk_string.as_ref()).wrap_err("When parsing Newick string")?;
  graph_from_newick(&nwk_graph)
}

pub fn nwk_read<N, E, D>(reader: impl Read) -> Result<Graph<N, E, D>, Report>
where
  N: GraphNode + NodeFromNwk + Named,
  E: GraphEdge + EdgeFromNwk,
  D: Sync + Send + Default,
{
  let nwk_graph = newick_from_reader(reader)?;
  graph_from_newick(&nwk_graph)
}

fn graph_from_newick<N, E, D>(nwk_graph: &NewickGraph) -> Result<Graph<N, E, D>, Report>
where
  N: GraphNode + NodeFromNwk + Named,
  E: GraphEdge + EdgeFromNwk,
  D: Sync + Send + Default,
{
  for (idx, node) in nwk_graph.nodes.iter().enumerate() {
    if node.hybrid.is_some() {
      return make_error!(
        "eNewick hybrid/reticulate node #{idx} '{}' is not supported. Treetime requires tree structure.",
        node.name.as_deref().unwrap_or("")
      );
    }
  }

  let mut graph = Graph::<N, E, D>::new();

  let mut node_keys: Vec<GraphNodeKey> = Vec::with_capacity(nwk_graph.nodes.len());
  for (nwk_idx, nwk_node) in nwk_graph.nodes.iter().enumerate() {
    let name: Option<&str> = nwk_node.name.as_deref().filter(|n| !n.is_empty());

    // Annotations parsed by util-newick but not yet wired through (kb/tickets/io-nwk-comment-dialect-parsing.md)
    let comments = btreemap! {};
    let node = N::from_nwk(name, &comments)
      .wrap_err_with(|| format!("When reading node #{nwk_idx} '{}'", name.unwrap_or_default()))?;
    node_keys.push(graph.add_node(node));
  }

  for (nwk_idx, nwk_edge) in nwk_graph.edges.iter().enumerate() {
    let source = node_keys
      .get(nwk_edge.parent)
      .ok_or_else(|| make_report!("When inserting edge {nwk_idx}: Node with index {} not found.", nwk_edge.parent))?;

    let target = node_keys
      .get(nwk_edge.child)
      .ok_or_else(|| make_report!("When inserting edge {nwk_idx}: Node with index {} not found.", nwk_edge.child))?;

    let edge = E::from_nwk(nwk_edge.data.branch_length)?;
    graph.add_edge(*source, *target, edge)?;
  }

  graph.build()?;

  assign_node_names(&graph)?;

  Ok(graph)
}

#[derive(Clone, SmartDefault)]
pub struct NwkWriteOptions {
  /// Annotation style: Plain suppresses annotations, Beast/Nhx emit structured comments.
  #[default(NwkStyle::Plain)]
  pub style: NwkStyle,

  /// Format node weights keeping this many significant digits
  pub weight_significant_digits: Option<u8>,

  /// Format node weights keeping this many decimal digits
  pub weight_decimal_digits: Option<i8>,
}

pub fn nwk_write_file<N, E, D>(
  filepath: impl AsRef<Path>,
  graph: &Graph<N, E, D>,
  options: &NwkWriteOptions,
) -> Result<(), Report>
where
  N: GraphNode + NodeToNwk,
  E: GraphEdge + EdgeToNwk,
  D: Sync + Send + Default,
{
  let mut f = create_file_or_stdout(filepath)?;
  nwk_write(&mut f, graph, options)?;
  writeln!(f)?;
  Ok(())
}

pub fn nwk_write_file_with<N, E, D>(
  filepath: impl AsRef<Path>,
  graph: &Graph<N, E, D>,
  options: &NwkWriteOptions,
  providers: &CommentProviders,
) -> Result<(), Report>
where
  N: GraphNode + NodeToNwk,
  E: GraphEdge + EdgeToNwk,
  D: Sync + Send + Default,
{
  let mut f = create_file_or_stdout(filepath)?;
  nwk_write_with(&mut f, graph, options, providers)?;
  writeln!(f)?;
  Ok(())
}

pub fn nwk_write_str<N, E, D>(graph: &Graph<N, E, D>, options: &NwkWriteOptions) -> Result<String, Report>
where
  N: GraphNode + NodeToNwk,
  E: GraphEdge + EdgeToNwk,
  D: Sync + Send + Default,
{
  let providers = CommentProviders::new();
  nwk_write_str_with(graph, options, &providers)
}

/// Return the Newick representation of a graph, augmented by external node comment providers.
pub fn nwk_write_str_with<N, E, D>(
  graph: &Graph<N, E, D>,
  options: &NwkWriteOptions,
  providers: &CommentProviders,
) -> Result<String, Report>
where
  N: GraphNode + NodeToNwk,
  E: GraphEdge + EdgeToNwk,
  D: Sync + Send + Default,
{
  let mut buf = Vec::new();
  nwk_write_with(&mut buf, graph, options, providers)?;
  Ok(String::from_utf8(buf)?)
}

pub fn nwk_write<N, E, D>(
  writer: &mut impl Write,
  graph: &Graph<N, E, D>,
  options: &NwkWriteOptions,
) -> Result<(), Report>
where
  N: GraphNode + NodeToNwk,
  E: GraphEdge + EdgeToNwk,
  D: Sync + Send + Default,
{
  let providers = CommentProviders::new();
  nwk_write_with(writer, graph, options, &providers)
}

/// Write a graph in Newick format, merging payload comments with external provider comments.
///
/// Provider comments override payload comments with the same key.
pub fn nwk_write_with<N, E, D>(
  writer: &mut impl Write,
  graph: &Graph<N, E, D>,
  options: &NwkWriteOptions,
  providers: &CommentProviders,
) -> Result<(), Report>
where
  N: GraphNode + NodeToNwk,
  E: GraphEdge + EdgeToNwk,
  D: Sync + Send + Default,
{
  let roots = graph.get_roots();
  if roots.is_empty() {
    return make_error!("When converting graph to Newick format: No roots found.");
  }

  if roots.len() > 1 {
    return make_error!("Multiple roots are not supported. Found {} roots", roots.len());
  }
  let root = &roots[0];

  let mut stack: Vec<(SafeNode<N>, Option<SafeEdge<E>>, usize)> = vec![(Arc::clone(root), None, 0)];
  while let Some((node, edge, child_visit)) = stack.pop() {
    let children: Vec<_> = graph.children_of(&node.read()).into_iter().collect();

    if child_visit < children.len() {
      stack.push((node, edge, child_visit + 1));

      if child_visit == 0 {
        write!(writer, "(")?;
      } else {
        write!(writer, ",")?;
      }

      let (child, child_edge) = &children[child_visit];
      stack.push((Arc::clone(child), Some(Arc::clone(child_edge)), 0));
    } else {
      if child_visit > 0 {
        write!(writer, ")")?;
      }

      let (node_key, name, mut comments) = {
        let node = node.read_arc();
        let node_key = node.key();
        let node_payload = node.payload().read_arc();
        let comments = node_payload.nwk_comments();
        let name = node_payload.nwk_name().map(|n| n.as_ref().to_owned());
        (node_key, name, comments)
      };
      comments.extend(providers.merged_comments(node_key)?);

      let weight = edge.and_then(|edge| edge.read_arc().payload().read_arc().nwk_weight());

      if let Some(name) = &name {
        write_label(writer, name)?;
      }

      if options.style != NwkStyle::Plain && !comments.is_empty() {
        let attrs = str_comments_to_newick_values(&comments);
        match options.style {
          NwkStyle::Beast => write_beast_attrs(writer, &attrs)?,
          NwkStyle::Nhx => write_nhx_attrs(writer, &attrs)?,
          NwkStyle::Plain => {},
        }
      }

      if let Some(weight) = weight {
        write!(writer, ":{}", format_weight(weight, options))?;
      }
    }
  }

  write!(writer, ";")?;

  Ok(())
}

fn str_comments_to_newick_values(comments: &BTreeMap<String, String>) -> BTreeMap<String, NewickValue> {
  comments
    .iter()
    .filter(|(_, val)| !val.is_empty())
    .map(|(key, val)| {
      let nwk_val = if val.eq_ignore_ascii_case("true") {
        NewickValue::Boolean(true)
      } else if val.eq_ignore_ascii_case("false") {
        NewickValue::Boolean(false)
      } else if let Ok(n) = val.parse::<f64>() {
        NewickValue::Number(n)
      } else {
        NewickValue::String(val.clone())
      };
      (key.clone(), nwk_val)
    })
    .collect()
}

pub fn format_weight(weight: f64, options: &NwkWriteOptions) -> String {
  if !weight.is_finite() {
    warn!("When converting graph to Newick: Weight is invalid: '{weight}'");
  }
  float_to_digits(
    weight,
    options.weight_significant_digits.or(Some(3)),
    options.weight_decimal_digits,
  )
}

/// Defines how to construct node when reading from Newick and Nexus files
pub trait NodeFromNwk: Sized {
  fn from_nwk(name: Option<impl AsRef<str>>, _: &BTreeMap<String, String>) -> Result<Self, Report>;
}

/// Defines how to display node information when writing to Newick and Nexus files
pub trait NodeToNwk {
  fn nwk_name(&self) -> Option<impl AsRef<str>>;

  fn nwk_comments(&self) -> BTreeMap<String, String> {
    BTreeMap::<String, String>::new()
  }
}

/// Return extra node comments for a graph node during Newick or Nexus serialization.
pub trait NodeCommentProvider {
  fn node_comments(&self, key: GraphNodeKey) -> Result<BTreeMap<String, String>, Report>;
}

/// Compose multiple node comment providers.
///
/// Providers are queried in insertion order. Later providers override earlier providers on key conflicts.
#[must_use]
#[derive(Default)]
pub struct CommentProviders<'a> {
  providers: Vec<&'a dyn NodeCommentProvider>,
}

impl<'a> CommentProviders<'a> {
  pub fn new() -> Self {
    Self::default()
  }

  pub fn with(mut self, provider: &'a dyn NodeCommentProvider) -> Self {
    self.providers.push(provider);
    self
  }

  pub fn merged_comments(&self, key: GraphNodeKey) -> Result<BTreeMap<String, String>, Report> {
    let mut comments = BTreeMap::new();
    for provider in &self.providers {
      comments.extend(provider.node_comments(key)?);
    }
    Ok(comments)
  }
}

/// Defines how to construct edge when reading from Newick and Nexus files
pub trait EdgeFromNwk: Sized {
  fn from_nwk(weight: Option<f64>) -> Result<Self, Report>;
}

/// Defines how to display edge information when writing to Newick and Nexus files
pub trait EdgeToNwk {
  fn nwk_weight(&self) -> Option<f64>;
}
