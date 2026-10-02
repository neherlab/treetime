use eyre::Report;
use serde::Serialize;
use std::collections::BTreeMap;
use treetime_graph::graph::Graph;
use treetime_graph::graph_traverse::GraphNodeForward;
use treetime_graph::node::GraphNodeKey;
use treetime_primitives::Seq;

pub(crate) fn reconstruct_preorder(
  graph: &Graph,
  include_leaves: bool,
  mut node_sequence: impl FnMut(&GraphNodeForward) -> Result<Seq, Report>,
) -> Result<Reconstruction, Report> {
  let mut sequences = BTreeMap::new();
  let mut emitted_nodes = Vec::new();
  graph.iter_depth_first_preorder_forward(|node| {
    let seq = node_sequence(&node)?;
    if include_leaves || !node.is_leaf {
      emitted_nodes.push(node.key);
    }
    sequences.insert(node.key, seq);
    Ok(())
  })?;
  Ok(Reconstruction {
    sequences,
    emitted_nodes,
  })
}

#[derive(Clone, Debug, Serialize)]
pub(crate) struct Reconstruction {
  pub sequences: BTreeMap<GraphNodeKey, Seq>,
  pub emitted_nodes: Vec<GraphNodeKey>,
}
