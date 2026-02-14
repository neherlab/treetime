use crate::auspice_types::{AuspiceTree, AuspiceTreeData, AuspiceTreeNode};
use crate::json::{JsonPretty, json_read, json_write};
use eyre::{Report, WrapErr};
use maplit::{btreemap, btreeset};
use parking_lot::RwLock;
use std::collections::VecDeque;
use std::io::Cursor;
use std::io::{Read, Write};
use std::path::Path;
use std::sync::Arc;
use treetime_graph::edge::{Edge, GraphEdge};
use treetime_graph::graph::Graph;
use treetime_graph::node::{GraphNode, GraphNodeKey, Node};
use treetime_utils::io::file::create_file_or_stdout;
use treetime_utils::io::file::open_file_or_stdin;

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

pub fn auspice_write_file<C, N, E, D>(filepath: impl AsRef<Path>, graph: &Graph<N, E, D>) -> Result<(), Report>
where
  N: GraphNode,
  E: GraphEdge,
  D: Sync + Send,
  C: AuspiceWrite<N, E, D>,
{
  let filepath = filepath.as_ref();
  let mut f = create_file_or_stdout(filepath)?;
  auspice_write::<C, _, _, _>(&mut f, graph)
    .wrap_err_with(|| format!("When writing Auspice v2 JSON file '{}'", filepath.display()))?;
  writeln!(f)?;
  Ok(())
}

pub fn auspice_write_str<C, N, E, D>(graph: &Graph<N, E, D>) -> Result<String, Report>
where
  N: GraphNode,
  E: GraphEdge,
  D: Sync + Send,
  C: AuspiceWrite<N, E, D>,
{
  let mut buf = Vec::new();
  auspice_write::<C, _, _, _>(&mut buf, graph).wrap_err("When writing Auspice v2 JSON string")?;
  Ok(String::from_utf8(buf)?)
}

pub fn auspice_write<C, N, E, D>(writer: &mut impl Write, graph: &Graph<N, E, D>) -> Result<(), Report>
where
  N: GraphNode,
  E: GraphEdge,
  D: Sync + Send,
  C: AuspiceWrite<N, E, D>,
{
  let tree = auspice_from_graph::<C, _, _, _>(graph)?;
  json_write(writer, &tree, JsonPretty(true)).wrap_err("When writing Auspice v2 JSON")
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

pub struct AuspiceGraphContext<'a, N, E, D>
where
  N: GraphNode,
  E: GraphEdge,
  D: Sync + Send,
{
  pub node_key: GraphNodeKey,
  pub node: &'a N,
  pub parent_key: Option<GraphNodeKey>,
  pub parent: Option<&'a N>,
  pub edge: Option<&'a E>,
  pub graph: &'a Graph<N, E, D>,
}

pub trait AuspiceWrite<N, E, D>: Sized
where
  N: GraphNode,
  E: GraphEdge,
  D: Sync + Send,
{
  fn new(graph: &Graph<N, E, D>) -> Result<Self, Report>;

  fn auspice_data_from_graph_data(&self, graph: &Graph<N, E, D>) -> Result<AuspiceTreeData, Report>;

  fn auspice_node_from_graph_components(
    &mut self,
    context: &AuspiceGraphContext<N, E, D>,
  ) -> Result<AuspiceTreeNode, Report>;
}

/// Convert graph to Auspice v2 JSON
pub fn auspice_from_graph<C, N, E, D>(graph: &Graph<N, E, D>) -> Result<AuspiceTree, Report>
where
  N: GraphNode,
  E: GraphEdge,
  D: Sync + Send,
  C: AuspiceWrite<N, E, D>,
{
  let mut converter = C::new(graph)?;
  let root = graph.get_exactly_one_root().wrap_err("When writing Auspice v2 JSON")?;

  // Pre-order iteration to construct the nodes from graph nodes and edges
  let mut node_map = {
    let mut node_map = btreemap! {};
    let mut queue = VecDeque::from([(Arc::clone(&root), None, None)]);
    while let Some((current_node, current_parent, current_edge)) = queue.pop_front() {
      let node_key = current_node.read_arc().key();
      let node = &*current_node.read_arc().payload().read_arc();
      let parent_key = current_parent
        .as_ref()
        .map(|parent: &Arc<RwLock<Node<N>>>| parent.read_arc().key());
      let parent = current_parent
        .as_ref()
        .map(|parent: &Arc<RwLock<Node<N>>>| parent.read_arc().payload().read_arc());
      let parent = parent.as_deref();
      let edge = current_edge
        .as_ref()
        .map(|edge: &Arc<RwLock<Edge<E>>>| edge.read_arc().payload().read_arc());
      let edge = edge.as_deref();
      let context = AuspiceGraphContext {
        node_key,
        node,
        parent_key,
        parent,
        edge,
        graph,
      };
      let current_tree_node = converter.auspice_node_from_graph_components(&context)?;
      for (child, edge) in graph.children_of(&current_node.read_arc()) {
        queue.push_back((child, Some(Arc::clone(&current_node)), Some(edge)));
      }
      node_map.insert(current_node.read_arc().key(), current_tree_node);
    }
    node_map
  };

  // Post-order traversal to populate .children array
  let mut visited = btreeset! {};
  let mut stack = vec![Arc::clone(&root)];
  while let Some(node) = stack.pop() {
    if visited.contains(&node.read_arc().key()) {
      let mut children_to_add = vec![];
      for (child, _) in graph.children_of(&node.read_arc()) {
        let child_key = child.read_arc().key();
        if let Some(child_node) = node_map.remove(&child_key) {
          children_to_add.push(child_node);
        }
      }
      if let Some(node) = node_map.get_mut(&node.read_arc().key()) {
        node.children.extend(children_to_add);
      }
    } else {
      visited.insert(node.read_arc().key());
      stack.push(Arc::clone(&node));
      for (child, _) in graph.children_of(&node.read_arc()) {
        stack.push(child);
      }
    }
  }

  let data = converter.auspice_data_from_graph_data(graph)?;
  let tree = node_map.remove(&root.read_arc().key()).unwrap();
  Ok(AuspiceTree { data, tree })
}
