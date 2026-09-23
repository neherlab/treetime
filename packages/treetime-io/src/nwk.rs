use crate::fasta::FastaRecord;
use eyre::{Report, WrapErr};
use log::warn;
use serde::{Deserialize, Serialize};
use smart_default::SmartDefault;
use std::collections::{BTreeMap, BTreeSet};
use std::io::{Read, Write};
use std::path::Path;
use treetime_graph::assign_node_names::assign_node_names;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_primitives::Seq;
use treetime_utils::fmt::float::float_to_digits;
use treetime_utils::io::file::create_file_or_stdout;
use treetime_utils::io::file::open_file_or_stdin;
use treetime_utils::make_error;
use treetime_utils::make_internal_report;
use treetime_utils::make_report;
pub use util_newick::NwkStyle;
use util_newick::{
  NewickGraph, NewickValue, newick_from_reader, newick_from_string, write_beast_attrs, write_label, write_nhx_attrs,
};

#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct NwkNodeMeta {
  name: Option<String>,
  confidence: Option<f64>,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct NwkFastaNodeInput {
  name: Option<String>,
  confidence: Option<f64>,
  aln: Option<Seq>,
  desc: Option<String>,
}

#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct NwkFastaEdgeInput {
  branch_length: Option<f64>,
}

#[derive(Debug)]
pub struct NwkFastaInput {
  pub graph: Graph,
  pub nodes: BTreeMap<GraphNodeKey, NwkFastaNodeInput>,
  pub edges: BTreeMap<GraphEdgeKey, NwkFastaEdgeInput>,
}

impl NwkFastaInput {
  pub fn from_parse_and_aln(parse: NwkParse, aln: Vec<FastaRecord>) -> Self {
    let NwkParse {
      graph,
      nodes: metas,
      branch_lengths,
    } = parse;

    let nodes = build_node_inputs(&graph, metas, aln);
    let edges = branch_lengths
      .into_iter()
      .map(|(key, branch_length)| (key, NwkFastaEdgeInput { branch_length }))
      .collect();

    Self { graph, nodes, edges }
  }

  pub fn names(&self) -> BTreeMap<GraphNodeKey, Option<String>> {
    self.nodes.iter().map(|(key, node)| (*key, node.name.clone())).collect()
  }

  pub fn branch_lengths(&self) -> BTreeMap<GraphEdgeKey, Option<f64>> {
    self
      .edges
      .iter()
      .map(|(key, edge)| (*key, edge.branch_length))
      .collect()
  }
}

pub fn nwk_fasta_node_inputs(
  graph: &Graph,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  aln: Vec<FastaRecord>,
) -> BTreeMap<GraphNodeKey, NwkFastaNodeInput> {
  let metas = names
    .iter()
    .map(|(key, name)| {
      (
        *key,
        NwkNodeMeta {
          name: name.clone(),
          confidence: None,
        },
      )
    })
    .collect();
  build_node_inputs(graph, metas, aln)
}

fn build_node_inputs(
  graph: &Graph,
  metas: BTreeMap<GraphNodeKey, NwkNodeMeta>,
  aln: Vec<FastaRecord>,
) -> BTreeMap<GraphNodeKey, NwkFastaNodeInput> {
  let mut records_by_name: BTreeMap<String, FastaRecord> = BTreeMap::new();
  for record in aln {
    records_by_name.entry(record.seq_name.clone()).or_insert(record);
  }

  let leaf_keys: BTreeSet<GraphNodeKey> = graph.get_leaves().map(|leaf| leaf.key()).collect();

  metas
    .into_iter()
    .map(|(key, meta)| {
      let record = if leaf_keys.contains(&key) {
        meta.name.as_deref().and_then(|name| records_by_name.remove(name))
      } else {
        None
      };
      let (aln, desc) = match record {
        Some(record) => (Some(record.seq), record.desc),
        None => (None, None),
      };
      (
        key,
        NwkFastaNodeInput {
          name: meta.name,
          confidence: meta.confidence,
          aln,
          desc,
        },
      )
    })
    .collect()
}

#[derive(Debug)]
pub struct NwkParse {
  pub graph: Graph,
  pub nodes: BTreeMap<GraphNodeKey, NwkNodeMeta>,
  pub branch_lengths: BTreeMap<GraphEdgeKey, Option<f64>>,
}

impl NwkParse {
  pub fn names(&self) -> BTreeMap<GraphNodeKey, Option<String>> {
    self.nodes.iter().map(|(key, meta)| (*key, meta.name.clone())).collect()
  }

  pub fn confidences(&self) -> BTreeMap<GraphNodeKey, Option<f64>> {
    self.nodes.iter().map(|(key, meta)| (*key, meta.confidence)).collect()
  }
}

pub fn nwk_read_file(filepath: impl AsRef<Path>) -> Result<NwkParse, Report> {
  let filepath = filepath.as_ref();
  nwk_read(open_file_or_stdin(&Some(filepath))?).wrap_err_with(|| format!("When reading file '{}'", filepath.display()))
}

pub fn nwk_read_str(nwk_string: impl AsRef<str>) -> Result<NwkParse, Report> {
  let nwk_graph = newick_from_string(nwk_string.as_ref()).wrap_err("When parsing Newick string")?;
  graph_from_newick(&nwk_graph)
}

fn nwk_read(reader: impl Read) -> Result<NwkParse, Report> {
  let nwk_graph = newick_from_reader(reader)?;
  graph_from_newick(&nwk_graph)
}

fn graph_from_newick(nwk_graph: &NewickGraph) -> Result<NwkParse, Report> {
  for (idx, node) in nwk_graph.nodes.iter().enumerate() {
    if node.hybrid.is_some() {
      return make_error!(
        "eNewick hybrid/reticulate node #{idx} '{}' is not supported. Treetime requires tree structure.",
        node.name.as_deref().unwrap_or("")
      );
    }
  }

  let mut graph = Graph::new();

  let mut node_keys: Vec<GraphNodeKey> = Vec::with_capacity(nwk_graph.nodes.len());
  let mut nodes: BTreeMap<GraphNodeKey, NwkNodeMeta> = BTreeMap::new();
  let mut branch_lengths: BTreeMap<GraphEdgeKey, Option<f64>> = BTreeMap::new();
  for nwk_node in &nwk_graph.nodes {
    let name: Option<&str> = nwk_node.name.as_deref().filter(|n| !n.is_empty());

    let key = graph.add_node();
    nodes.insert(
      key,
      NwkNodeMeta {
        name: name.map(ToOwned::to_owned),
        confidence: nwk_node.confidence,
      },
    );
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

  let names = nodes.iter().map(|(key, meta)| (*key, meta.name.clone())).collect();
  let names = assign_node_names(names, &graph)?;
  for (key, name) in names {
    if let Some(meta) = nodes.get_mut(&key) {
      meta.name = name;
    }
  }

  Ok(NwkParse {
    graph,
    nodes,
    branch_lengths,
  })
}

#[derive(Clone, SmartDefault)]
pub struct NwkWriteOptions {
  #[default(NwkStyle::Plain)]
  pub style: NwkStyle,

  pub weight_significant_digits: Option<u8>,

  pub weight_decimal_digits: Option<i8>,
}

pub fn nwk_write_file(
  filepath: impl AsRef<Path>,
  graph: &Graph,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  weights: &BTreeMap<GraphEdgeKey, Option<f64>>,
  options: &NwkWriteOptions,
) -> Result<(), Report> {
  let mut f = create_file_or_stdout(filepath)?;
  nwk_write(&mut f, graph, names, weights, options)?;
  writeln!(f)?;
  Ok(())
}

pub fn nwk_write_file_with(
  filepath: impl AsRef<Path>,
  graph: &Graph,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  weights: &BTreeMap<GraphEdgeKey, Option<f64>>,
  options: &NwkWriteOptions,
  providers: &CommentProviders,
) -> Result<(), Report> {
  let mut f = create_file_or_stdout(filepath)?;
  nwk_write_with(&mut f, graph, names, weights, options, providers)?;
  writeln!(f)?;
  Ok(())
}

pub fn nwk_write_str(
  graph: &Graph,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  weights: &BTreeMap<GraphEdgeKey, Option<f64>>,
  options: &NwkWriteOptions,
) -> Result<String, Report> {
  let providers = CommentProviders::new();
  nwk_write_str_with(graph, names, weights, options, &providers)
}

pub(crate) fn nwk_write_str_with(
  graph: &Graph,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  weights: &BTreeMap<GraphEdgeKey, Option<f64>>,
  options: &NwkWriteOptions,
  providers: &CommentProviders,
) -> Result<String, Report> {
  let mut buf = Vec::new();
  nwk_write_with(&mut buf, graph, names, weights, options, providers)?;
  Ok(String::from_utf8(buf)?)
}

pub fn nwk_write(
  writer: &mut impl Write,
  graph: &Graph,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  weights: &BTreeMap<GraphEdgeKey, Option<f64>>,
  options: &NwkWriteOptions,
) -> Result<(), Report> {
  let providers = CommentProviders::new();
  nwk_write_with(writer, graph, names, weights, options, &providers)
}

pub(crate) fn nwk_write_with(
  writer: &mut impl Write,
  graph: &Graph,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  weights: &BTreeMap<GraphEdgeKey, Option<f64>>,
  options: &NwkWriteOptions,
  providers: &CommentProviders,
) -> Result<(), Report> {
  if graph.num_roots() == 0 {
    return make_error!("When converting graph to Newick format: No roots found.");
  }

  if graph.num_roots() > 1 {
    return make_error!("Multiple roots are not supported. Found {} roots", graph.num_roots());
  }
  let root_key = graph.root_key()?;

  let mut stack: Vec<(GraphNodeKey, Option<GraphEdgeKey>, usize)> = vec![(root_key, None, 0)];
  while let Some((node_key, edge_key, child_visit)) = stack.pop() {
    let node = graph
      .get_node(node_key)
      .ok_or_else(|| make_internal_report!("Node {node_key} not found in graph"))?;
    let children: Vec<_> = graph.children_keys_of(node).collect();

    if child_visit < children.len() {
      stack.push((node_key, edge_key, child_visit + 1));

      if child_visit == 0 {
        write!(writer, "(")?;
      } else {
        write!(writer, ",")?;
      }

      let (child_key, child_edge_key) = children[child_visit];
      stack.push((child_key, Some(child_edge_key), 0));
    } else {
      if child_visit > 0 {
        write!(writer, ")")?;
      }

      let name = names[&node_key].clone();
      let comments = providers.merged_comments(node_key)?;

      let weight = edge_key.and_then(|edge_key| weights[&edge_key]);

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

#[cfg_attr(
  dylint_lib = "treetime_lints",
  expect(
    error_dropped_by_pattern,
    reason = "a comment value that is not a number is kept as a string"
  )
)]
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

pub(crate) fn format_weight(weight: f64, options: &NwkWriteOptions) -> String {
  if !weight.is_finite() {
    warn!("When converting graph to Newick: Weight is invalid: '{weight}'");
  }
  float_to_digits(
    weight,
    options.weight_significant_digits.or(Some(3)),
    options.weight_decimal_digits,
  )
}

pub trait NodeCommentProvider {
  fn node_comments(&self, key: GraphNodeKey) -> Result<BTreeMap<String, String>, Report>;
}

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

  fn merged_comments(&self, key: GraphNodeKey) -> Result<BTreeMap<String, String>, Report> {
    let mut comments = BTreeMap::new();
    for provider in &self.providers {
      comments.extend(provider.node_comments(key)?);
    }
    Ok(comments)
  }
}
