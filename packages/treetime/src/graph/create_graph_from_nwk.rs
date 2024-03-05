use crate::graph::assign_node_names::assign_node_names;
use crate::graph::edge::GraphEdge;
use crate::graph::graph::Graph;
use crate::graph::node::{GraphNode, GraphNodeKey, NodeType};
use crate::io::nwk::read_nwk_file;
use color_eyre::Section;
use eyre::{eyre, Report, WrapErr};
use indexmap::IndexMap;
use std::path::Path;

pub fn create_graph_from_nwk<N, E, P>(tree_path: P) -> Result<Graph<N, E>, Report>
where
  N: GraphNode,
  E: GraphEdge,
  P: AsRef<Path>,
{
  let nwk_tree = read_nwk_file(tree_path)
    .wrap_err("When parsing input tree")
    .with_section(|| "Note: only Newick format is currently supported")?;

  let mut graph = Graph::<N, E>::new();

  // Insert nodes
  let mut index_map = IndexMap::<usize, GraphNodeKey>::new(); // Map of internal `nwk` node indices to `Graph` node indices
  for (nwk_idx, nwk_node) in nwk_tree.g.raw_nodes().iter().enumerate() {
    // Attempt to parse weight as float. If not a float, then it's a named leaf node, otherwise - internal node.
    let inserted_node_idx = if let Ok(weight) = nwk_node.weight.parse::<f64>() {
      let node = N::internal("", 0.0);
      graph.add_node(node)
    } else {
      let name = nwk_node.weight.as_str();
      graph.add_node(N::leaf(name))
    };

    index_map.insert(nwk_idx, inserted_node_idx);
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

  // Mark roots
  for root in graph.get_roots() {
    let root = root.write().payload();
    let mut root = root.write();
    root.set_node_type(NodeType::Root(0.0));
  }

  assign_node_names(&graph);

  Ok(graph)
}
