use crate::types::NewickGraph;
use eyre::{Report, eyre};
use std::collections::BTreeSet;

pub(crate) fn validate_graph(graph: &NewickGraph) -> Result<(), Report> {
  let node_count = graph.nodes.len();
  if graph.root >= node_count {
    return Err(eyre!(
      "The root is node {}, but the graph has {node_count} nodes",
      graph.root
    ));
  }

  let mut parent_counts = vec![0_usize; node_count];
  for (idx, edge) in graph.edges.iter().enumerate() {
    if edge.parent >= node_count || edge.child >= node_count {
      return Err(eyre!(
        "Edge {idx} connects node {} to node {}, but the graph has {node_count} nodes",
        edge.parent,
        edge.child
      ));
    }
    parent_counts[edge.child] += 1;
  }

  let mut listed = vec![false; graph.edges.len()];
  for (idx, node) in graph.nodes.iter().enumerate() {
    let mut children = BTreeSet::new();
    for &edge_idx in &node.children {
      let edge = graph.edges.get(edge_idx).ok_or_else(|| {
        eyre!(
          "{} lists edge {edge_idx} as a child edge, but the graph has {} edges",
          describe_node(graph, idx),
          graph.edges.len()
        )
      })?;
      if edge.parent != idx {
        return Err(eyre!(
          "{} lists edge {edge_idx} as a child edge, but the edge starts at {}",
          describe_node(graph, idx),
          describe_node(graph, edge.parent)
        ));
      }
      if !children.insert(edge.child) {
        return Err(eyre!(
          "{} has more than one edge to {}",
          describe_node(graph, idx),
          describe_node(graph, edge.child)
        ));
      }
      listed[edge_idx] = true;
    }
    if node.hybrid.is_none() && parent_counts[idx] > 1 {
      return Err(eyre!(
        "{} has {} parents, but only a hybrid node can have more than one parent",
        describe_node(graph, idx),
        parent_counts[idx]
      ));
    }
  }

  if let Some(edge_idx) = listed.iter().position(|is_listed| !is_listed) {
    return Err(eyre!(
      "Edge {edge_idx} is not listed among the child edges of {}",
      describe_node(graph, graph.edges[edge_idx].parent)
    ));
  }

  if parent_counts[graph.root] > 0 {
    return Err(eyre!("The root, {}, has a parent", describe_node(graph, graph.root)));
  }

  if postorder(graph).is_none() {
    return Err(eyre!("The graph contains a cycle"));
  }

  Ok(())
}

pub(crate) fn postorder(graph: &NewickGraph) -> Option<Vec<usize>> {
  let mut marks = vec![Mark::New; graph.nodes.len()];
  let mut order = Vec::with_capacity(graph.nodes.len());
  let mut stack = vec![graph.root];
  while let Some(&node) = stack.last() {
    match marks[node] {
      Mark::New => {
        marks[node] = Mark::Open;
        for &edge_idx in &graph.nodes[node].children {
          let child = graph.edges[edge_idx].child;
          match marks[child] {
            Mark::Open => return None,
            Mark::New => stack.push(child),
            Mark::Done => {},
          }
        }
      },
      Mark::Open => {
        stack.pop();
        marks[node] = Mark::Done;
        order.push(node);
      },
      Mark::Done => {
        stack.pop();
      },
    }
  }
  Some(order)
}

pub(crate) fn describe_node(graph: &NewickGraph, idx: usize) -> String {
  match graph.nodes.get(idx).and_then(|node| node.name()) {
    Some(name) => format!("node {idx} ('{name}')"),
    None => format!("node {idx}"),
  }
}

#[derive(Clone, Copy)]
enum Mark {
  New,
  Open,
  Done,
}
