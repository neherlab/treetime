use crate::graph::edge::GraphEdge;
use crate::graph::graph::Graph;
use crate::graph::node::{GraphNode, Node};
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

  let root = root.read();
  node_to_nwk_string(writer, graph, &root, None, options)?;
  write!(writer, ";")?;

  Ok(())
}

fn node_to_nwk_string<N, E>(
  writer: &mut impl Write,
  graph: &Graph<N, E>,
  node: &Node<N>,
  weight: Option<f64>,
  options: &WriteNwkOptions,
) -> Result<(), Report>
where
  N: GraphNode,
  E: GraphEdge,
{
  let outbound_edge_keys = node.outbound();

  if !outbound_edge_keys.is_empty() {
    write!(writer, "(")?;

    let mut first = true;
    for outbound_edge_key in outbound_edge_keys {
      let edge = graph.get_edge(*outbound_edge_key).unwrap();

      if first {
        first = false;
      } else {
        write!(writer, ",")?;
      }

      {
        let child_key = edge.read().target();
        let child = graph.get_node(child_key).unwrap();
        let child = child.read();
        let weight = edge.read().payload().read().weight();
        node_to_nwk_string(writer, graph, &child, Some(weight), options)?;
      }
    }
    write!(writer, ")")?;
  }

  let (name, comments) = {
    let node_payload = node.payload();
    let node_payload = node_payload.read();
    (node_payload.name().to_owned(), node_payload.nwk_comments())
  };
  write!(writer, "{name}")?;

  if let Some(weight) = weight {
    write!(writer, ":{}", format_weight(weight, options))?;
  }

  if !comments.is_empty() {
    let comments = comments.iter().map(|(key, val)| format!("[&{key}=\"{val}\"]")).join("");
    write!(writer, "{comments}")?;
  }

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
