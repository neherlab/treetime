use eyre::{Report, WrapErr};
use log::warn;
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
use treetime_utils::{make_error, make_report};
use util_newick::{
  InternalLabel, LabelSide, NewickComment, NewickEdgeData, NewickGraph, NewickNodeData, NewickReadOptions,
  NewickWarning, NewickWriteOptions, NodeComment, NumberFormat, ReadMode, newick_from_reader, newick_to_writer,
};
pub use util_newick::{NewickDialect, NewickValue};

pub const TREE_DIALECT_DEFAULT: NewickDialect = NewickDialect::BEAST;

pub fn nwk_read_file(filepath: impl AsRef<Path>) -> Result<NwkParse, Report> {
  read_file_with(filepath, nwk_read)
}

pub fn nwk_read(reader: impl Read) -> Result<NwkParse, Report> {
  let tree = newick_from_reader(reader, &tree_read_options(TREE_DIALECT_DEFAULT)).wrap_err("When reading Newick")?;
  log_read_warnings(&tree.warnings);
  graph_from_newick(&tree.graph).wrap_err("When reading Newick")
}

pub(crate) fn tree_read_options(dialect: NewickDialect) -> NewickReadOptions {
  NewickReadOptions {
    dialect,
    mode: ReadMode::Tolerant,
    internal_label: InternalLabel::Auto,
    underscores_as_spaces: false,
  }
}

pub(crate) fn log_read_warnings(warnings: &[NewickWarning]) {
  for warning in warnings {
    warn!("When reading the tree: {warning}");
  }
}

pub(crate) fn graph_from_newick(nwk_graph: &NewickGraph) -> Result<NwkParse, Report> {
  reject_network_nodes(nwk_graph)?;
  let mut graph = Graph::new();

  let mut node_keys: Vec<GraphNodeKey> = Vec::with_capacity(nwk_graph.node_count());
  let mut nodes: BTreeMap<GraphNodeKey, NwkNodeMeta> = BTreeMap::new();
  let mut branch_lengths: BTreeMap<GraphEdgeKey, Option<f64>> = BTreeMap::new();
  for (_, nwk_node) in nwk_graph.nodes() {
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

  for (nwk_idx, nwk_edge) in nwk_graph.edges() {
    let source = node_keys.get(nwk_edge.parent()).ok_or_else(|| {
      make_report!(
        "When inserting edge {nwk_idx}: Node with index {} not found.",
        nwk_edge.parent()
      )
    })?;

    let target = node_keys.get(nwk_edge.child()).ok_or_else(|| {
      make_report!(
        "When inserting edge {nwk_idx}: Node with index {} not found.",
        nwk_edge.child()
      )
    })?;

    let edge_key = graph.add_edge(*source, *target)?;
    branch_lengths.insert(edge_key, nwk_edge.data().branch_length());
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

fn reject_network_nodes(nwk_graph: &NewickGraph) -> Result<(), Report> {
  let network_node = nwk_graph
    .nodes()
    .map(|(node, data)| (data, nwk_graph.parent_edges(node).len()))
    .find(|&(_, parents)| parents > 1);
  match network_node {
    Some((data, parents)) => {
      let tag = data.hybrid().map(|hybrid| hybrid.tag(false)).unwrap_or_default();
      let label = format!("{}{tag}", data.name().unwrap_or_default());
      make_error!("The tree contains the network node '{label}' with {parents} parents; TreeTime reads trees only")
    },
    None => Ok(()),
  }
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

#[derive(Clone, Debug, Default, PartialEq, Eq)]
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
  let graph = newick_graph(tree, names, weights, comments, options.style)?;
  newick_to_writer(&mut writer, &graph, &newick_write_options(options))
}

pub type NwkNodeComments = BTreeMap<GraphNodeKey, Vec<(String, NewickValue)>>;

pub(crate) fn newick_graph(
  tree: &TreeView<'_>,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  weights: &BTreeMap<GraphEdgeKey, Option<f64>>,
  comments: &NwkNodeComments,
  style: NwkStyle,
) -> Result<NewickGraph, Report> {
  let node_data = |key: GraphNodeKey| {
    let mut data = NewickNodeData::new();
    if let Some(name) = &names[&key] {
      data = data.with_name(name.clone());
    }
    match (style, comments.get(&key)) {
      (NwkStyle::Beast, Some(pairs)) if !pairs.is_empty() => data.with_comment(NodeComment::new(
        LabelSide::AfterLabel,
        NewickComment::Beast(pairs.clone()),
      )),
      (NwkStyle::Nhx, Some(pairs)) if !pairs.is_empty() => data.with_comment(NodeComment::new(
        LabelSide::AfterLabel,
        NewickComment::Nhx(pairs.clone()),
      )),
      (NwkStyle::Plain | NwkStyle::Beast | NwkStyle::Nhx, _) => data,
    }
  };
  let mut graph = NewickGraph::new(node_data(tree.root()));
  let mut indices: BTreeMap<GraphNodeKey, usize> = BTreeMap::new();
  indices.insert(tree.root(), graph.root());
  for &key in tree.preorder() {
    let parent = indices[&key];
    for &(child, edge_key) in tree.children(key) {
      let mut edge = NewickEdgeData::new();
      if let Some(weight) = weights[&edge_key] {
        edge = edge.with_length(weight);
      }
      let child_idx = graph.add_child(parent, edge, node_data(child))?;
      indices.insert(child, child_idx);
    }
  }
  Ok(graph)
}

pub(crate) fn newick_write_options(options: &NwkWriteOptions) -> NewickWriteOptions {
  NewickWriteOptions {
    numbers: NumberFormat {
      significant_digits: options.weight_significant_digits.or(Some(3)),
      decimal_digits: options.weight_decimal_digits,
      point_zero: false,
    },
    ..NewickWriteOptions::new(options.style.dialect())
  }
}

pub(crate) fn format_weight(weight: f64, options: &NwkWriteOptions) -> String {
  float_to_digits(
    weight,
    options.weight_significant_digits.or(Some(3)),
    options.weight_decimal_digits,
  )
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, SmartDefault)]
pub enum NwkStyle {
  #[default]
  Plain,
  Beast,
  Nhx,
}

impl NwkStyle {
  pub const fn dialect(self) -> NewickDialect {
    match self {
      Self::Plain => NewickDialect::CLASSIC,
      Self::Beast => NewickDialect::BEAST,
      Self::Nhx => NewickDialect::NHX,
    }
  }
}

#[derive(Clone, SmartDefault)]
pub struct NwkWriteOptions {
  #[default(NwkStyle::Plain)]
  pub style: NwkStyle,

  pub weight_significant_digits: Option<u8>,

  pub weight_decimal_digits: Option<i8>,
}
