use eyre::{Report, WrapErr};
use log::warn;
use smart_default::SmartDefault;
use std::collections::BTreeMap;
use std::io::{Read, Write};
use std::path::Path;
use std::sync::Arc;
use treetime_graph::assign_node_names::assign_node_names;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::{Graph, SafeEdge, SafeNode};
use treetime_graph::node::GraphNodeKey;
use treetime_utils::fmt::float::float_to_digits;
use treetime_utils::io::file::create_file_or_stdout;
use treetime_utils::io::file::open_file_or_stdin;
use treetime_utils::make_error;
use treetime_utils::make_report;
pub use util_newick::NwkStyle;
use util_newick::{
  NewickGraph, NewickValue, newick_from_reader, newick_from_string, write_beast_attrs, write_label, write_nhx_attrs,
};

/// A parsed Newick tree: the graph together with the per-node input-tree branch support and names.
///
/// `confidences` is keyed by the graph's own node keys and holds each node's Newick branch support
/// (bootstrap or posterior), with `None` where a node carried no confidence annotation. It lets a
/// consumer read each node's input branch support as a value threaded from the parse rather than off
/// the node payload.
///
/// `names` is keyed by the graph's own node keys and holds each node's name: the parsed name for
/// named nodes and the synthetic `NODE_xxxxx` name that `assign_node_names` assigns to internals,
/// with `None` where a node has no name. It lets a consumer read each node's name as a value
/// threaded from the parse rather than off the node payload.
#[derive(Debug)]
pub struct NwkParse<D = ()>
where
  D: Sync + Send,
{
  pub graph: Graph<D>,
  pub confidences: BTreeMap<GraphNodeKey, Option<f64>>,
  pub names: BTreeMap<GraphNodeKey, Option<String>>,
  /// Each edge's raw input-tree branch length, keyed by the graph's own edge keys, with `None`
  /// where an edge carried no `:length`. Threaded from the parse so a consumer reads each branch
  /// length as a value rather than off the edge payload.
  pub branch_lengths: BTreeMap<GraphEdgeKey, Option<f64>>,
}

pub fn nwk_read_file<D>(filepath: impl AsRef<Path>) -> Result<NwkParse<D>, Report>
where
  D: Sync + Send + Default,
{
  let filepath = filepath.as_ref();
  nwk_read(open_file_or_stdin(&Some(filepath))?).wrap_err_with(|| format!("When reading file '{}'", filepath.display()))
}

pub fn nwk_read_str<D>(nwk_string: impl AsRef<str>) -> Result<NwkParse<D>, Report>
where
  D: Sync + Send + Default,
{
  let nwk_graph = newick_from_string(nwk_string.as_ref()).wrap_err("When parsing Newick string")?;
  graph_from_newick(&nwk_graph)
}

pub fn nwk_read<D>(reader: impl Read) -> Result<NwkParse<D>, Report>
where
  D: Sync + Send + Default,
{
  let nwk_graph = newick_from_reader(reader)?;
  graph_from_newick(&nwk_graph)
}

fn graph_from_newick<D>(nwk_graph: &NewickGraph) -> Result<NwkParse<D>, Report>
where
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

  let mut graph = Graph::<D>::new();

  let mut node_keys: Vec<GraphNodeKey> = Vec::with_capacity(nwk_graph.nodes.len());
  let mut confidences: BTreeMap<GraphNodeKey, Option<f64>> = BTreeMap::new();
  let mut names: BTreeMap<GraphNodeKey, Option<String>> = BTreeMap::new();
  let mut branch_lengths: BTreeMap<GraphEdgeKey, Option<f64>> = BTreeMap::new();
  for nwk_node in &nwk_graph.nodes {
    let name: Option<&str> = nwk_node.name.as_deref().filter(|n| !n.is_empty());

    let key = graph.add_node();
    confidences.insert(key, nwk_node.confidence);
    names.insert(key, name.map(ToOwned::to_owned));
    node_keys.push(key);
  }

  for (nwk_idx, nwk_edge) in nwk_graph.edges.iter().enumerate() {
    let source = node_keys.get(nwk_edge.parent).ok_or_else(|| {
      make_report!(
        "When inserting edge {nwk_idx}: Node with index {} not found.",
        nwk_edge.parent
      )
    })?;

    let target = node_keys.get(nwk_edge.child).ok_or_else(|| {
      make_report!(
        "When inserting edge {nwk_idx}: Node with index {} not found.",
        nwk_edge.child
      )
    })?;

    let edge_key = graph.add_edge(*source, *target)?;
    branch_lengths.insert(edge_key, nwk_edge.data.branch_length);
  }

  graph.build()?;

  let names = assign_node_names(names, &graph)?;

  Ok(NwkParse {
    graph,
    confidences,
    names,
    branch_lengths,
  })
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

pub fn nwk_write_file<D>(
  filepath: impl AsRef<Path>,
  graph: &Graph<D>,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  weights: &BTreeMap<GraphEdgeKey, Option<f64>>,
  options: &NwkWriteOptions,
) -> Result<(), Report>
where
  D: Sync + Send,
{
  let mut f = create_file_or_stdout(filepath)?;
  nwk_write(&mut f, graph, names, weights, options)?;
  writeln!(f)?;
  Ok(())
}

pub fn nwk_write_file_with<D>(
  filepath: impl AsRef<Path>,
  graph: &Graph<D>,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  weights: &BTreeMap<GraphEdgeKey, Option<f64>>,
  options: &NwkWriteOptions,
  providers: &CommentProviders,
) -> Result<(), Report>
where
  D: Sync + Send,
{
  let mut f = create_file_or_stdout(filepath)?;
  nwk_write_with(&mut f, graph, names, weights, options, providers)?;
  writeln!(f)?;
  Ok(())
}

pub fn nwk_write_str<D>(
  graph: &Graph<D>,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  weights: &BTreeMap<GraphEdgeKey, Option<f64>>,
  options: &NwkWriteOptions,
) -> Result<String, Report>
where
  D: Sync + Send,
{
  let providers = CommentProviders::new();
  nwk_write_str_with(graph, names, weights, options, &providers)
}

/// Return the Newick representation of a graph, augmented by external node comment providers.
pub fn nwk_write_str_with<D>(
  graph: &Graph<D>,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  weights: &BTreeMap<GraphEdgeKey, Option<f64>>,
  options: &NwkWriteOptions,
  providers: &CommentProviders,
) -> Result<String, Report>
where
  D: Sync + Send,
{
  let mut buf = Vec::new();
  nwk_write_with(&mut buf, graph, names, weights, options, providers)?;
  Ok(String::from_utf8(buf)?)
}

pub fn nwk_write<D>(
  writer: &mut impl Write,
  graph: &Graph<D>,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  weights: &BTreeMap<GraphEdgeKey, Option<f64>>,
  options: &NwkWriteOptions,
) -> Result<(), Report>
where
  D: Sync + Send,
{
  let providers = CommentProviders::new();
  nwk_write_with(writer, graph, names, weights, options, &providers)
}

/// Write a graph in Newick format, taking node names and edge weights from explicit value maps and
/// node comments from external comment providers.
///
/// `names` supplies each node's display label and `weights` each edge's branch weight, both keyed by
/// the graph's own keys and kept as `Option` so a missing label writes no name and a missing weight
/// writes no `:weight`. Comments come solely from the providers.
pub fn nwk_write_with<D>(
  writer: &mut impl Write,
  graph: &Graph<D>,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  weights: &BTreeMap<GraphEdgeKey, Option<f64>>,
  options: &NwkWriteOptions,
  providers: &CommentProviders,
) -> Result<(), Report>
where
  D: Sync + Send,
{
  let roots = graph.get_roots();
  if roots.is_empty() {
    return make_error!("When converting graph to Newick format: No roots found.");
  }

  if roots.len() > 1 {
    return make_error!("Multiple roots are not supported. Found {} roots", roots.len());
  }
  let root = &roots[0];

  let mut stack: Vec<(SafeNode, Option<SafeEdge>, usize)> = vec![(Arc::clone(root), None, 0)];
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

      let node_key = node.read_arc().key();
      let name = names[&node_key].clone();
      let comments = providers.merged_comments(node_key)?;

      let weight = edge
        .map(|edge| edge.read_arc().key())
        .and_then(|edge_key| weights[&edge_key]);

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
