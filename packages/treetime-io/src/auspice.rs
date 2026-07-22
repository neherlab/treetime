use crate::auspice_types::{AuspiceTree, AuspiceTreeNode};
use eyre::{Report, WrapErr};
use std::collections::VecDeque;
use std::io::Cursor;
use std::io::{Read, Write};
use std::path::Path;
use treetime_graph::edge::GraphEdge;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNode;
use treetime_utils::io::file::create_file_or_stdout;
use treetime_utils::io::file::open_file_or_stdin;
use treetime_utils::io::json::{JsonPretty, json_read, json_write};

pub fn auspice_read_file<C, N, E, D>(filepath: impl AsRef<Path>) -> Result<Graph<N, E, D>, Report>
where
  N: GraphNode,
  E: GraphEdge,
  D: Sync + Send,
  C: AuspiceRead<N, E, D>,
{
  let filepath = filepath.as_ref();
  auspice_read::<C, _, _, _>(open_file_or_stdin(&Some(filepath))?)
    .wrap_err_with(|| format!("When reading Auspice v2 JSON file '{}'", filepath.display()))
}

pub fn auspice_read_str<C, N, E, D>(auspice_string: impl AsRef<str>) -> Result<Graph<N, E, D>, Report>
where
  N: GraphNode,
  E: GraphEdge,
  D: Sync + Send,
  C: AuspiceRead<N, E, D>,
{
  let auspice_string = auspice_string.as_ref();
  auspice_read::<C, _, _, _>(Cursor::new(auspice_string)).wrap_err("When reading Auspice v2 JSON string")
}

pub fn auspice_read<C, N, E, D>(reader: impl Read) -> Result<Graph<N, E, D>, Report>
where
  N: GraphNode,
  E: GraphEdge,
  D: Sync + Send,
  C: AuspiceRead<N, E, D>,
{
  let tree: AuspiceTree = json_read(reader).wrap_err("When reading Auspice v2 JSON")?;
  auspice_to_graph::<C, _, _, _>(&tree).wrap_err("When converting Auspice v2 JSON to graph")
}

pub fn auspice_write_file(filepath: impl AsRef<Path>, tree: &AuspiceTree) -> Result<(), Report> {
  let filepath = filepath.as_ref();
  let mut f = create_file_or_stdout(filepath)?;
  auspice_write(&mut f, tree)
    .wrap_err_with(|| format!("When writing Auspice v2 JSON file '{}'", filepath.display()))?;
  writeln!(f)?;
  Ok(())
}

pub fn auspice_write_str(tree: &AuspiceTree) -> Result<String, Report> {
  let mut buf = Vec::new();
  auspice_write(&mut buf, tree).wrap_err("When writing Auspice v2 JSON string")?;
  Ok(String::from_utf8(buf)?)
}

pub fn auspice_write(writer: &mut impl Write, tree: &AuspiceTree) -> Result<(), Report> {
  json_write(writer, tree, JsonPretty(true)).wrap_err("When writing Auspice v2 JSON")
}

pub struct AuspiceTreeContext<'a> {
  pub node: &'a AuspiceTreeNode,
  pub parent: Option<&'a AuspiceTreeNode>,
  pub tree: &'a AuspiceTree,
}

impl AuspiceTreeContext<'_> {
  pub fn branch_length(&self) -> Option<f64> {
    let parent_div = self.parent.and_then(|parent| parent.node_attrs.div);
    let div = self.node.node_attrs.div;
    match (parent_div, div) {
      (Some(parent_div), Some(div)) => Some(div - parent_div),
      (None, Some(div)) => Some(div),
      (Some(_parent_div), None) => None,
      (None, None) => None,
    }
  }
}

pub trait AuspiceRead<N, E, D>: Sized
where
  N: GraphNode,
  E: GraphEdge,
  D: Sync + Send,
{
  fn new(tree: &AuspiceTree) -> Result<Self, Report>;

  fn auspice_data_to_graph_data(&mut self, tree: &AuspiceTree) -> Result<D, Report>;

  fn auspice_node_to_graph_components(&mut self, context: &AuspiceTreeContext) -> Result<(N, E), Report>;
}

/// Convert Auspice v2 JSON to graph
pub fn auspice_to_graph<C, N, E, D>(tree: &AuspiceTree) -> Result<Graph<N, E, D>, Report>
where
  N: GraphNode,
  E: GraphEdge,
  D: Sync + Send,
  C: AuspiceRead<N, E, D>,
{
  let mut converter = C::new(tree)?;
  let mut graph = Graph::<N, E, D>::with_data(converter.auspice_data_to_graph_data(tree)?);
  let mut queue = VecDeque::from([(None, &tree.tree, None)]);
  while let Some((parent_key, node, parent)) = queue.pop_front() {
    let (graph_node, graph_edge) =
      converter.auspice_node_to_graph_components(&AuspiceTreeContext { node, parent, tree })?;
    let node_key = graph.add_node(graph_node);
    if let Some(parent_key) = parent_key {
      graph.add_edge(parent_key, node_key, graph_edge)?;
    }
    for child in &node.children {
      queue.push_back((Some(node_key), child, Some(node)));
    }
  }
  graph.build()?;
  Ok(graph)
}
