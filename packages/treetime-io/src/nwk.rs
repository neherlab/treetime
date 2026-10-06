use eyre::{Report, WrapErr};
use serde::{Deserialize, Serialize};
use smart_default::SmartDefault;
use std::collections::BTreeMap;
use std::io::{Read, Write};
use std::path::Path;
use treetime_graph::assign_node_names::{AssignedNodeNames, assign_node_names};
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_graph::tree_view::TreeView;
use treetime_utils::fmt::float::float_to_digits;
use treetime_utils::io::file::{read_file_with, write_file_with};
use treetime_utils::make_error;
use treetime_utils::make_report;
use util_newick::{
  NewickGraph, NewickReadOptions, newick_from_reader, write_beast_attrs, write_label, write_nhx_attrs,
};
pub use util_newick::{NewickValue, NwkStyle};

pub fn nwk_read_file(filepath: impl AsRef<Path>) -> Result<NwkParse, Report> {
  read_file_with(filepath, nwk_read)
}

pub fn nwk_read(reader: impl Read) -> Result<NwkParse, Report> {
  newick_from_reader(reader, &NewickReadOptions::default())
    .and_then(|nwk_graph| graph_from_newick(&nwk_graph))
    .wrap_err("When reading Newick")
}

pub(crate) fn graph_from_newick(nwk_graph: &NewickGraph) -> Result<NwkParse, Report> {
  let mut graph = Graph::new();

  let mut node_keys: Vec<GraphNodeKey> = Vec::with_capacity(nwk_graph.nodes.len());
  let mut nodes: BTreeMap<GraphNodeKey, NwkNodeMeta> = BTreeMap::new();
  let mut branch_lengths: BTreeMap<GraphEdgeKey, Option<f64>> = BTreeMap::new();
  for nwk_node in &nwk_graph.nodes {
    let name: Option<&str> = nwk_node.name().filter(|n| !n.is_empty());

    let key = graph.add_node();
    nodes.insert(
      key,
      NwkNodeMeta {
        name: name.map(ToOwned::to_owned),
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
  let AssignedNodeNames { names, duplicate_names } = assign_node_names(names, &graph)?;
  for (key, name) in names {
    if let Some(meta) = nodes.get_mut(&key) {
      meta.name = name;
    }
  }

  Ok(NwkParse {
    graph,
    nodes,
    branch_lengths,
    duplicate_names,
  })
}

#[derive(Debug)]
pub struct NwkParse {
  pub graph: Graph,
  pub nodes: BTreeMap<GraphNodeKey, NwkNodeMeta>,
  pub branch_lengths: BTreeMap<GraphEdgeKey, Option<f64>>,
  pub duplicate_names: Vec<String>,
}

impl NwkParse {
  pub fn names(&self) -> BTreeMap<GraphNodeKey, Option<String>> {
    self.nodes.iter().map(|(key, meta)| (*key, meta.name.clone())).collect()
  }
}

#[derive(Clone, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct NwkNodeMeta {
  name: Option<String>,
}

pub fn nwk_write_file(
  filepath: impl AsRef<Path>,
  tree: &TreeView<'_>,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  weights: &BTreeMap<GraphEdgeKey, Option<f64>>,
  options: &NwkWriteOptions,
  comments: &NwkNodeComments,
) -> Result<(), Report> {
  write_file_with(filepath, |writer| {
    nwk_write(&mut *writer, tree, names, weights, options, comments)?;
    writeln!(writer)?;
    Ok(())
  })
}

pub fn nwk_write_str(
  tree: &TreeView<'_>,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  weights: &BTreeMap<GraphEdgeKey, Option<f64>>,
  options: &NwkWriteOptions,
  comments: &NwkNodeComments,
) -> Result<String, Report> {
  let mut buffer = Vec::new();
  nwk_write(&mut buffer, tree, names, weights, options, comments)?;
  Ok(String::from_utf8(buffer)?)
}

pub fn nwk_write(
  mut writer: impl Write,
  tree: &TreeView<'_>,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  weights: &BTreeMap<GraphEdgeKey, Option<f64>>,
  options: &NwkWriteOptions,
  comments: &NwkNodeComments,
) -> Result<(), Report> {
  write_nwk_tree(&mut writer, tree, names, weights, options, comments)
    .and_then(|()| Ok(writer.write_all(b";")?))
    .wrap_err("When writing Newick")
}

pub type NwkNodeComments = BTreeMap<GraphNodeKey, Vec<(String, NewickValue)>>;

pub(crate) fn write_nwk_tree(
  writer: &mut impl Write,
  tree: &TreeView<'_>,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  weights: &BTreeMap<GraphEdgeKey, Option<f64>>,
  options: &NwkWriteOptions,
  comments: &NwkNodeComments,
) -> Result<(), Report> {
  let mut stack: Vec<(GraphNodeKey, Option<GraphEdgeKey>, usize)> = vec![(tree.root(), None, 0)];
  while let Some((node_key, edge_key, child_visit)) = stack.pop() {
    let children = tree.children(node_key);

    if child_visit < children.len() {
      stack.push((node_key, edge_key, child_visit + 1));
      writer.write_all(if child_visit == 0 { b"(" } else { b"," })?;
      let (child_key, child_edge_key) = children[child_visit];
      stack.push((child_key, Some(child_edge_key), 0));
      continue;
    }

    if child_visit > 0 {
      writer.write_all(b")")?;
    }

    if let Some(name) = &names[&node_key] {
      write_label(writer, name)?;
    }

    if options.style != NwkStyle::Plain
      && let Some(node_comments) = comments.get(&node_key)
      && !node_comments.is_empty()
    {
      let attrs = node_comments.iter().map(|(key, value)| (key.as_str(), value));
      match options.style {
        NwkStyle::Beast => write_beast_attrs(writer, attrs)?,
        NwkStyle::Nhx => write_nhx_attrs(writer, attrs)?,
        NwkStyle::Plain => {},
      }
    }

    if let Some(weight) = edge_key.and_then(|edge_key| weights[&edge_key]) {
      if !weight.is_finite() {
        let name = names[&node_key].as_deref().unwrap_or("");
        return make_error!("The branch above node '{name}' has the length {weight}, which Newick cannot represent");
      }
      write!(writer, ":{}", format_weight(weight, options))?;
    }
  }
  Ok(())
}

pub(crate) fn format_weight(weight: f64, options: &NwkWriteOptions) -> String {
  float_to_digits(
    weight,
    options.weight_significant_digits.or(Some(3)),
    options.weight_decimal_digits,
  )
}

#[derive(Clone, SmartDefault)]
pub struct NwkWriteOptions {
  #[default(NwkStyle::Plain)]
  pub style: NwkStyle,

  pub weight_significant_digits: Option<u8>,

  pub weight_decimal_digits: Option<i8>,
}
