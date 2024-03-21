use crate::graph::assign_node_names::assign_node_names;
use crate::graph::edge::GraphEdge;
use crate::graph::graph::Graph;
use crate::graph::node::{GraphNode, GraphNodeKey};
use crate::io::file::open_file_or_stdin;
use crate::io::nwk::read_nwk;
use eyre::{eyre, Report, WrapErr};
use indexmap::IndexMap;
use petgraph::visit::IntoNodeReferences;
use petgraph::Direction;
use std::io::Cursor;
use std::path::Path;

pub fn create_graph_from_nwk_file<N, E>(filepath: impl AsRef<Path>) -> Result<Graph<N, E>, Report>
where
  N: GraphNode,
  E: GraphEdge,
{
  let filepath = filepath.as_ref();
  create_graph_from_nwk_reader(open_file_or_stdin(&Some(filepath))?)
    .wrap_err_with(|| format!("When reading file '{filepath:#?}'"))
}

pub fn create_graph_from_nwk_str<N, E>(nwk_string: impl AsRef<str>) -> Result<Graph<N, E>, Report>
where
  N: GraphNode,
  E: GraphEdge,
{
  let nwk_string = nwk_string.as_ref();
  create_graph_from_nwk_reader(Cursor::new(nwk_string))
    .wrap_err_with(|| format!("When reading Newick string:\n    '{nwk_string}'"))
}

pub fn create_graph_from_nwk_reader<N, E>(s: impl std::io::Read) -> Result<Graph<N, E>, Report>
where
  N: GraphNode,
  E: GraphEdge,
{
  let nwk_tree = read_nwk(s).wrap_err("When parsing Newick")?;

  let mut graph = Graph::<N, E>::new();

  // Insert nodes
  let mut index_map = IndexMap::<usize, GraphNodeKey>::new(); // Map of internal `nwk` node indices to `Graph` node indices
  for (nwk_idx, nwk_node) in nwk_tree.g.node_references() {
    let n_edges_incoming = nwk_tree.g.edges_directed(nwk_idx, Direction::Incoming).count();
    let n_edges_outgoing = nwk_tree.g.edges_directed(nwk_idx, Direction::Outgoing).count();

    let inserted_node_idx = match (n_edges_incoming, n_edges_outgoing) {
      (0, _) => graph.add_node(N::root(nwk_node)),
      (_, 0) => graph.add_node(N::leaf(nwk_node)),
      (_, _) => graph.add_node(N::internal(nwk_node)),
    };

    index_map.insert(nwk_idx.index(), inserted_node_idx);
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

    graph.add_edge(*source, *target, E::new(weight))?;
  }

  graph.build()?;

  assign_node_names(&graph);

  Ok(graph)
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::graph::edge::Weighted;
  use crate::graph::node::{Named, NodeType, WithNwkComments};
  use crate::io::nwk::{write_nwk_str, WriteNwkOptions};
  use eyre::Report;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::fmt::{Display, Formatter};

  #[derive(Clone, Debug, PartialEq, Eq)]
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

  #[derive(Clone, Debug, PartialEq)]
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
  fn test_nwk_read_write() -> Result<(), Report> {
    let input = "((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root;";

    // TODO: Try to implement weights attached to root node. In a general graph, root node cannot contain
    // incoming edges. There is nowhere they could come from. So currently root node weight is not stored.
    //  Note the 0.01 at the end in this example:
    // let input = "((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;";

    let graph = create_graph_from_nwk_str::<Node, Edge>(input)?;

    let output = write_nwk_str(&graph, &WriteNwkOptions::default())?;

    assert_eq!(input, output);
    Ok(())
  }
}
