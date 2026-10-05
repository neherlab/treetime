use eyre::{Report, WrapErr};
use log::warn;
use serde::{Deserialize, Serialize};
use smart_default::SmartDefault;
use std::collections::BTreeMap;
use std::fmt;
use std::io::{Read, Write};
use std::path::Path;
use treetime_graph::assign_node_names::assign_node_names;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_graph::tree_view::TreeView;
use treetime_utils::fmt::float::float_to_digits;
use treetime_utils::io::file::{read_file_with, write_file_with};
use treetime_utils::make_error;
use treetime_utils::make_report;
pub use util_newick::NwkStyle;
use util_newick::{NewickGraph, NewickValue, newick_from_reader, write_beast_attrs, write_label, write_nhx_attrs};

pub const NEWICK_EXTENSIONS: [&str; 4] = ["nwk", "newick", "tree", "tre"];

pub fn nwk_read_file(filepath: impl AsRef<Path>) -> Result<NwkParse, Report> {
  read_file_with(filepath, nwk_read)
}

pub fn nwk_read(reader: impl Read) -> Result<NwkParse, Report> {
  newick_from_reader(reader)
    .and_then(|nwk_graph| graph_from_newick(&nwk_graph))
    .wrap_err("When reading Newick")
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

#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct NwkNodeMeta {
  name: Option<String>,
  confidence: Option<f64>,
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
  let mut text = String::new();
  write_nwk_text(&mut text, tree, names, weights, options, comments).wrap_err("When writing Newick")?;
  Ok(text)
}

pub fn nwk_write(
  mut writer: impl Write,
  tree: &TreeView<'_>,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  weights: &BTreeMap<GraphEdgeKey, Option<f64>>,
  options: &NwkWriteOptions,
  comments: &NwkNodeComments,
) -> Result<(), Report> {
  let text = nwk_write_str(tree, names, weights, options, comments)?;
  writer.write_all(text.as_bytes()).wrap_err("When writing Newick")
}

pub type NwkNodeComments = BTreeMap<GraphNodeKey, BTreeMap<String, String>>;

fn write_nwk_text(
  writer: &mut impl fmt::Write,
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

      if let Some(name) = &names[&node_key] {
        write_label(writer, name)?;
      }

      if options.style != NwkStyle::Plain
        && let Some(node_comments) = comments.get(&node_key)
        && !node_comments.is_empty()
      {
        let attrs = str_comments_to_newick_values(node_comments);
        match options.style {
          NwkStyle::Beast => write_beast_attrs(writer, &attrs)?,
          NwkStyle::Nhx => write_nhx_attrs(writer, &attrs)?,
          NwkStyle::Plain => {},
        }
      }

      if let Some(weight) = edge_key.and_then(|edge_key| weights[&edge_key]) {
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

#[derive(Clone, SmartDefault)]
pub struct NwkWriteOptions {
  #[default(NwkStyle::Plain)]
  pub style: NwkStyle,

  pub weight_significant_digits: Option<u8>,

  pub weight_decimal_digits: Option<i8>,
}
