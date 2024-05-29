use crate::graph::edge::GraphEdge;
use crate::graph::graph::{Graph, SafeEdge, SafeNode};
use crate::graph::node::GraphNode;
use crate::io::file::create_file;
use crate::io::fs::read_file_to_string;
use crate::make_error;
use crate::utils::float_fmt::float_to_digits;
use bio::io::newick;
use bio_types::phylogeny::Tree;
use eyre::{Report, WrapErr};
use itertools::Itertools;
use log::warn;
use smart_default::SmartDefault;
use std::io::{Read, Write};
use std::path::Path;
use std::sync::Arc;

pub fn read_nwk_file(nwk_file_path: impl AsRef<Path>) -> Result<Tree, Report> {
  let nwk_file_path = nwk_file_path.as_ref();
  let nwk_str = read_file_to_string(nwk_file_path)?;
  read_nwk(nwk_str.as_bytes()).wrap_err_with(|| format!("When parsing Newick file {nwk_file_path:#?}"))
}

pub fn read_nwk(reader: impl Read) -> Result<Tree, Report> {
  let mut nwk_tree = newick::read(reader)?;

  nwk_tree.g.node_weights_mut().for_each(|weight| {
    if weight == "N/A" {
      *weight = "".to_owned();
    }
  });

  Ok(nwk_tree)
}

#[derive(Clone, SmartDefault)]
pub struct WriteNwkOptions {
  /// Format node weights keeping this many significant digits
  pub weight_significant_digits: Option<u8>,

  /// Format node weights keeping this many decimal digits
  pub weight_decimal_digits: Option<i8>,
}

pub fn write_nwk_file<N, E>(
  filepath: &impl AsRef<Path>,
  graph: &Graph<N, E>,
  options: &WriteNwkOptions,
) -> Result<(), Report>
where
  N: GraphNode,
  E: GraphEdge,
{
  let mut f = create_file(filepath)?;
  write_nwk_writer(&mut f, graph, options)?;
  writeln!(f)?;
  Ok(())
}

pub fn write_nwk_str<N, E>(graph: &Graph<N, E>, options: &WriteNwkOptions) -> Result<String, Report>
where
  N: GraphNode,
  E: GraphEdge,
{
  let mut buf = Vec::new();
  write_nwk_writer(&mut buf, graph, options)?;
  Ok(String::from_utf8(buf)?)
}

pub fn write_nwk_writer<N, E>(
  writer: &mut impl Write,
  graph: &Graph<N, E>,
  options: &WriteNwkOptions,
) -> Result<(), Report>
where
  N: GraphNode,
  E: GraphEdge,
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
        (node_payload.name().to_owned(), node_payload.nwk_comments())
      };

      let weight = edge.map(|edge| edge.read_arc().payload().read().weight());

      write!(writer, "{name}")?;

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

pub fn format_weight(weight: f64, options: &WriteNwkOptions) -> String {
  if !weight.is_finite() {
    warn!("When converting graph to Newick: Weight is invalid: '{weight}'");
  }

  // if let Some(precision) = options.weight_precision {
  //   return format!("{weight:.precision$}");
  // }

  let digits = options.weight_significant_digits.unwrap_or(3);
  float_to_digits(
    weight,
    options.weight_significant_digits.or(Some(3)),
    options.weight_decimal_digits,
  )
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::graph::create_graph_from_nwk::create_graph_from_nwk_str;
  use crate::graph::edge::Weighted;
  use crate::graph::node::{Named, NodeType, WithNwkComments};
  use eyre::Report;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use serde::{Deserialize, Serialize};
  use std::fmt::{Display, Formatter};

  #[derive(Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
  pub struct Node {
    pub name: String,
    pub node_type: NodeType,
  }

  impl Node {
    pub fn new(name: impl AsRef<str>, node_type: NodeType) -> Self {
      Self {
        name: name.as_ref().to_owned(),
        node_type,
      }
    }
  }

  impl GraphNode for Node {
    fn root(name: &str) -> Self {
      Self::new(name, NodeType::Root(name.to_owned()))
    }

    fn internal(name: &str) -> Self {
      Self::new(name, NodeType::Internal(name.to_owned()))
    }

    fn leaf(name: &str) -> Self {
      Self::new(name, NodeType::Leaf(name.to_owned()))
    }

    fn set_node_type(&mut self, node_type: NodeType) {
      self.node_type = node_type;
    }
  }

  impl WithNwkComments for Node {}

  impl Named for Node {
    fn name(&self) -> &str {
      &self.name
    }

    fn set_name(&mut self, name: &str) {
      self.name = name.to_owned();
    }
  }

  impl Display for Node {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
      match &self.node_type {
        NodeType::Root(name) => write!(f, "{name:}"),
        NodeType::Internal(name) => write!(f, "{name:}"),
        NodeType::Leaf(name) => write!(f, "{name}"),
      }
    }
  }

  #[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
  pub struct Edge {
    pub weight: f64,
  }

  impl GraphEdge for Edge {
    fn new(weight: f64) -> Self {
      Self { weight }
    }
  }

  impl Weighted for Edge {
    fn weight(&self) -> f64 {
      self.weight
    }
  }

  impl Display for Edge {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
      write!(f, "{:}", &self.weight)
    }
  }

  #[rstest]
  fn test_write_nwk() -> Result<(), Report> {
    let input = "((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root;";
    let graph = create_graph_from_nwk_str::<Node, Edge>(input)?;
    let output = write_nwk_str(&graph, &WriteNwkOptions::default())?;
    assert_eq!(input, output);
    Ok(())
  }
}
