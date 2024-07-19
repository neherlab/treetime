use crate::graph::assign_node_names::assign_node_names;
use crate::graph::edge::GraphEdge;
use crate::graph::graph::{Graph, SafeEdge, SafeNode};
use crate::graph::node::GraphNodeKey;
use crate::graph::node::{GraphNode, Named};
use crate::io::file::create_file;
use crate::io::file::open_file_or_stdin;
use crate::make_error;
use crate::o;
use crate::utils::float_fmt::float_to_digits;
use bio::io::newick;
use eyre::{eyre, Report, WrapErr};
use indexmap::IndexMap;
use itertools::Itertools;
use log::warn;
use maplit::btreemap;
use petgraph::visit::IntoNodeReferences;
use smart_default::SmartDefault;
use std::collections::BTreeMap;
use std::io::Cursor;
use std::io::{Read, Write};
use std::path::Path;
use std::sync::Arc;

pub fn nwk_read_file<N, E>(filepath: impl AsRef<Path>) -> Result<Graph<N, E>, Report>
where
  N: GraphNode + NodeFromNwk + Named,
  E: GraphEdge + EdgeFromNwk,
{
  let filepath = filepath.as_ref();
  nwk_read(open_file_or_stdin(&Some(filepath))?).wrap_err_with(|| format!("When reading file '{filepath:#?}'"))
}

pub fn nwk_read_str<N, E>(nwk_string: impl AsRef<str>) -> Result<Graph<N, E>, Report>
where
  N: GraphNode + NodeFromNwk + Named,
  E: GraphEdge + EdgeFromNwk,
{
  let nwk_string = nwk_string.as_ref();
  nwk_read(Cursor::new(nwk_string)).wrap_err_with(|| format!("When reading Newick string:\n    '{nwk_string}'"))
}

pub fn nwk_read<N, E>(reader: impl Read) -> Result<Graph<N, E>, Report>
where
  N: GraphNode + NodeFromNwk + Named,
  E: GraphEdge + EdgeFromNwk,
{
  let mut nwk_tree = newick::read(reader)?;

  nwk_tree.g.node_weights_mut().for_each(|weight| {
    if weight == "N/A" {
      *weight = "".to_owned();
    }
  });

  let mut graph = Graph::<N, E>::new();

  // Insert nodes
  let mut index_map = IndexMap::<usize, GraphNodeKey>::new(); // Map of internal `nwk` node indices to `Graph` node indices
  for (nwk_idx, nwk_node) in nwk_tree.g.node_references() {
    // Discard node names which are parseable to a number. These are not names, but weights.
    // And we don't need them here. Weights are collected onto the edges later.
    let mut nwk_node: String = nwk_node.to_owned();
    if nwk_node.parse::<f64>().is_ok() {
      nwk_node = o!("");
    };

    let comments = btreemap! {}; // TODO: parse nwk comments
    let node = N::from_nwk(&nwk_node, &comments).wrap_err_with(|| format!("When reading node {nwk_node}"))?;
    let node_key = graph.add_node(node);
    index_map.insert(nwk_idx.index(), node_key);
  }

  // Insert edges
  for (nwk_idx, nwk_edge) in nwk_tree.g.raw_edges().iter().enumerate() {
    let weight: f64 = nwk_edge.weight as f64;
    let source: usize = nwk_edge.source().index();
    let target: usize = nwk_edge.target().index();

    let source = index_map
      .get(&source)
      .ok_or_else(|| eyre!("When inserting edge {nwk_idx}: Node with index {source} not found."))?;

    let target = index_map
      .get(&target)
      .ok_or_else(|| eyre!("When inserting edge {nwk_idx}: Node with index {target} not found."))?;

    let edge = E::from_nwk(Some(weight))?;
    graph.add_edge(*source, *target, edge)?;
  }

  graph.build()?;

  assign_node_names(&graph);

  Ok(graph)
}

#[derive(Clone, SmartDefault)]
pub struct NwkWriteOptions {
  /// Format node weights keeping this many significant digits
  pub weight_significant_digits: Option<u8>,

  /// Format node weights keeping this many decimal digits
  pub weight_decimal_digits: Option<i8>,
}

pub fn nwk_write_file<N, E>(
  filepath: impl AsRef<Path>,
  graph: &Graph<N, E>,
  options: &NwkWriteOptions,
) -> Result<(), Report>
where
  N: GraphNode + NodeToNwk,
  E: GraphEdge + EdgeToNwk,
{
  let mut f = create_file(filepath)?;
  nwk_write(&mut f, graph, options)?;
  writeln!(f)?;
  Ok(())
}

pub fn nwk_write_str<N, E>(graph: &Graph<N, E>, options: &NwkWriteOptions) -> Result<String, Report>
where
  N: GraphNode + NodeToNwk,
  E: GraphEdge + EdgeToNwk,
{
  let mut buf = Vec::new();
  nwk_write(&mut buf, graph, options)?;
  Ok(String::from_utf8(buf)?)
}

pub fn nwk_write<N, E>(writer: &mut impl Write, graph: &Graph<N, E>, options: &NwkWriteOptions) -> Result<(), Report>
where
  N: GraphNode + NodeToNwk,
  E: GraphEdge + EdgeToNwk,
{
  let roots = graph.get_roots();
  if roots.is_empty() {
    return make_error!("When converting graph to Newick format: No roots found.");
  }

  let root = {
    if roots.len() > 1 {
      unimplemented!("Multiple roots are not supported yet");
    }
    &roots[0]
  };

  let mut stack: Vec<(SafeNode<N>, Option<SafeEdge<E>>, usize)> = vec![(Arc::clone(root), None, 0)];
  while let Some((node, edge, child_visit)) = stack.pop() {
    let children = graph.children_of(&node.read()).into_iter().collect_vec();

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

      let (name, comments) = {
        let node_payload = node.read_arc().payload().read_arc();
        (
          node_payload.nwk_name().map(|n| n.as_ref().to_owned()),
          node_payload.nwk_comments(),
        )
      };

      let weight = edge.and_then(|edge| edge.read_arc().payload().read_arc().nwk_weight());

      if let Some(name) = name {
        write!(writer, "{name}")?;
      }

      if let Some(weight) = weight {
        write!(writer, ":{}", format_weight(weight, options))?;
      }

      if !comments.is_empty() {
        let comments = comments.iter().map(|(key, val)| format!("[&{key}=\"{val}\"]")).join("");
        write!(writer, "{comments}")?;
      }
    }
  }

  write!(writer, ";")?;

  Ok(())
}

pub fn format_weight(weight: f64, options: &NwkWriteOptions) -> String {
  if !weight.is_finite() {
    warn!("When converting graph to Newick: Weight is invalid: '{weight}'");
  }
  let digits = options.weight_significant_digits.unwrap_or(3);
  float_to_digits(
    weight,
    options.weight_significant_digits.or(Some(3)),
    options.weight_decimal_digits,
  )
}

/// Defines how to construct node when reading from Newick and Nexus files
pub trait NodeFromNwk: Sized {
  fn from_nwk(name: impl AsRef<str>, comments: &BTreeMap<String, String>) -> Result<Self, Report>;
}

/// Defines how to display node information when writing to Newick and Nexus files
pub trait NodeToNwk {
  fn nwk_name(&self) -> Option<impl AsRef<str>>;

  fn nwk_comments(&self) -> BTreeMap<String, String> {
    BTreeMap::<String, String>::new()
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

#[cfg(test)]
mod tests {
  use super::*;
  use crate::graph::graph::tests::{TestEdge, TestNode};
  use eyre::Report;
  use pretty_assertions::assert_eq;

  #[test]
  fn test_nwk_roundtrip() -> Result<(), Report> {
    let input = "((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root;";
    let graph = nwk_read_str::<TestNode, TestEdge>(input)?;
    let output = nwk_write_str(&graph, &NwkWriteOptions::default())?;
    assert_eq!(input, output);
    Ok(())
  }
}
