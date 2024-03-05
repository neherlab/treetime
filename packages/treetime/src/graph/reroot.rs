use crate::graph::edge::{invert_edge, GraphEdge};
use crate::graph::find_paths::find_paths;
use crate::graph::graph::Graph;
use crate::graph::node::{GraphNode, Node};
use eyre::Report;
use parking_lot::RwLock;
use std::sync::Arc;

/// Change root of the graph.
pub fn reroot<N: GraphNode, E: GraphEdge>(
  graph: &mut Graph<N, E>,
  old_root: &Arc<RwLock<Node<N>>>,
  new_root: &Arc<RwLock<Node<N>>>,
) -> Result<(), Report> {
  assert!(graph.num_roots() <= 1, "Multiple roots are not supported yet");

  let new_root_key = new_root.read().key();
  let old_root_key = old_root.read().key();

  if new_root_key == old_root_key {
    // The new root is the same as old. Nothing to do.
    return Ok(());
  }

  // Find paths from the old root to the new desired root
  let paths = find_paths(graph, old_root, new_root)?;

  // Invert every edge on the path from old to new root.
  // This will make the desired new root into an actual root. The old root might no longer be a root.
  for edge in &paths {
    invert_edge(graph, edge);
  }

  // Recalculate some bookkeeping
  graph.build()
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::graph::examples::get_example_tree;
  use crate::graph::node::GraphNodeKey;
  use crate::io::file::create_file;
  use eyre::Report;
  use pretty_assertions::assert_eq;
  use rstest::rstest;

  #[rstest]
  fn reroots() -> Result<(), Report> {
    let mut graph = get_example_tree()?;
    graph.print_graph(create_file("../../tmp__/graph_input.dot")?)?;

    let roots = graph.get_roots();
    let old_root = roots.first().unwrap();
    let old_root_key = old_root.read().key();

    let new_root = graph.get_node(GraphNodeKey(5)).unwrap();
    let new_root_key = new_root.read().key();

    reroot(&mut graph, old_root, &new_root).unwrap();
    graph.print_graph(create_file("../../tmp__/graph_output.dot")?)?;

    let roots = graph.get_roots();
    let actual_root = roots.first().unwrap();
    let actual_root_key = actual_root.read().key();

    assert_eq!(*new_root.read(), *actual_root.read());

    Ok(())
  }
}
