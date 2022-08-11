use crate::graph::edge::GraphEdge;
use crate::graph::graph::Graph;
use crate::graph::node::{GraphNode, Node};
use crate::io::fs::read_file_to_string;
use crate::make_error;
use crate::utils::float_fmt::float_to_significant_digits;
use bio::io::newick;
use bio_types::phylogeny::Tree;
use eyre::{Report, WrapErr};
use std::io::{Read, Write};
use std::path::Path;

pub fn read_nwk_file(nwk_file_path: impl AsRef<Path>) -> Result<Tree, Report> {
  let nwk_file_path = nwk_file_path.as_ref();
  let nwk_str = read_file_to_string(&nwk_file_path)?;
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

pub fn write_nwk<N, E>(writer: &mut impl Write, graph: &Graph<N, E>) -> Result<(), Report>
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
  node_to_nwk_string(writer, &root)?;
  writeln!(writer, ";")?;

  Ok(())
}

pub fn to_nwk_string<N, E>(graph: &Graph<N, E>) -> Result<String, Report>
where
  N: GraphNode,
  E: GraphEdge,
{
  let mut buf = Vec::new();
  write_nwk(&mut buf, graph)?;
  Ok(String::from_utf8(buf)?)
}

fn node_to_nwk_string<N, E>(writer: &mut impl Write, node: &Node<N, E>) -> Result<(), Report>
where
  N: GraphNode,
  E: GraphEdge,
{
  let edges = node.outbound();

  if !edges.is_empty() {
    write!(writer, "(")?;

    let mut first = true;
    for edge in edges.iter() {
      let edge = edge.as_ref();

      if first {
        first = false;
      } else {
        write!(writer, ",")?;
      }

      node_to_nwk_string(writer, &edge.target().read())?;

      let weight = edge.payload().read().weight();
      if weight.is_finite() {
        let weight = float_to_significant_digits(weight, 3);
        write!(writer, ":{weight}")?;
      }
    }
    write!(writer, ")")?;
  }

  let node_payload = node.payload();
  let node_payload = node_payload.read();
  let name = node_payload.name();
  write!(writer, "{name}")?;

  Ok(())
}
