use crate::model::data::NewickEdgeData;
use crate::model::graph::NewickGraph;
use crate::model::traverse::Postorder;
use eyre::{Report, eyre};
use std::collections::BTreeSet;

pub(crate) fn validate_graph(graph: &NewickGraph) -> Result<(), Report> {
  let node_count = graph.node_count();
  if graph.root() >= node_count {
    return Err(eyre!(
      "The root is node {}, but the graph has {node_count} nodes",
      graph.root()
    ));
  }
  validate_edges(graph)?;
  for node in 0..node_count {
    validate_node(graph, node)?;
  }
  if !graph.parent_edges(graph.root()).is_empty() {
    return Err(eyre!("The root, {}, has a parent", describe_node(graph, graph.root())));
  }
  if Postorder::has_cycle(graph) {
    return Err(eyre!("The graph contains a cycle"));
  }
  let mut reached = vec![false; node_count];
  for node in graph.preorder() {
    reached[node] = true;
  }
  if let Some(node) = reached.iter().position(|&is_reached| !is_reached) {
    return Err(eyre!("{} is not reachable from the root", describe_node(graph, node)));
  }
  Ok(())
}

pub(crate) fn describe_node(graph: &NewickGraph, node: usize) -> String {
  match (node < graph.node_count()).then(|| graph.node(node).name()).flatten() {
    Some(name) => format!("node {node} ('{name}')"),
    None => format!("node {node}"),
  }
}

fn validate_edges(graph: &NewickGraph) -> Result<(), Report> {
  let node_count = graph.node_count();
  for (idx, edge) in graph.edges() {
    if edge.parent() >= node_count || edge.child() >= node_count {
      return Err(eyre!(
        "Edge {idx} connects node {} to node {}, but the graph has {node_count} nodes",
        edge.parent(),
        edge.child()
      ));
    }
    let listed_by_parent = graph.child_edges(edge.parent()).iter().filter(|&&e| e == idx).count();
    let listed_by_child = graph.parent_edges(edge.child()).iter().filter(|&&e| e == idx).count();
    if listed_by_parent != 1 || listed_by_child != 1 {
      return Err(eyre!(
        "Edge {idx} must be listed once among the child edges of {} and once among the parent edges of {}",
        describe_node(graph, edge.parent()),
        describe_node(graph, edge.child())
      ));
    }
  }
  Ok(())
}

fn validate_node(graph: &NewickGraph, node: usize) -> Result<(), Report> {
  let edge_count = graph.edge_count();
  let lists = [(graph.child_edges(node), true), (graph.parent_edges(node), false)];
  for (edges, is_child_list) in lists {
    for &edge in edges {
      let entry = (edge < edge_count).then(|| graph.edge(edge)).ok_or_else(|| {
        eyre!(
          "{} lists edge {edge}, but the graph has {edge_count} edges",
          describe_node(graph, node)
        )
      })?;
      let end = if is_child_list { entry.parent() } else { entry.child() };
      if end != node {
        return Err(eyre!(
          "{} lists edge {edge}, which does not end at it",
          describe_node(graph, node)
        ));
      }
    }
  }
  let mut children = BTreeSet::new();
  for child in graph.children(node) {
    if !children.insert(child) {
      return Err(eyre!(
        "{} has more than one edge to {}",
        describe_node(graph, node),
        describe_node(graph, child)
      ));
    }
  }
  let parent_count = graph.parent_edges(node).len();
  let data = graph.node(node);
  if data.hybrid().is_none() && parent_count > 1 {
    return Err(eyre!(
      "{} has {parent_count} parents, but only a hybrid node can have more than one parent",
      describe_node(graph, node)
    ));
  }
  if data.hybrid().is_some() && !data.comments().is_empty() {
    return Err(eyre!(
      "{} is a hybrid node with comments; the comments of a hybrid node belong to the edges into its occurrences",
      describe_node(graph, node)
    ));
  }
  let in_edges: Vec<&NewickEdgeData> = graph
    .parent_edges(node)
    .iter()
    .map(|&edge| graph.edge(edge).data())
    .chain((node == graph.root()).then(|| graph.root_edge()))
    .collect();
  if data.hybrid().is_none() && in_edges.iter().any(|edge| !edge.occurrence_comments().is_empty()) {
    return Err(eyre!(
      "The edge into {} has occurrence comments, but only an edge into a hybrid node can have them",
      describe_node(graph, node)
    ));
  }
  Ok(())
}
