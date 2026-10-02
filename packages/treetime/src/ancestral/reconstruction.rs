use crate::ancestral::sample::SampleMode;
use eyre::Report;
use serde::Serialize;
use std::collections::BTreeMap;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_primitives::Seq;

pub(crate) fn emitted_nodes(graph: &Graph, include_leaves: bool) -> Result<Vec<GraphNodeKey>, Report> {
  let mut emitted = Vec::new();
  graph.iter_depth_first_preorder_forward(|node| {
    if include_leaves || !node.is_leaf {
      emitted.push(node.key);
    }
    Ok(())
  })?;
  Ok(emitted)
}

pub(crate) fn sample_internal_sequences(
  graph: &Graph,
  sample_mode: SampleMode,
  mut sample: impl FnMut(GraphNodeKey) -> Seq,
) -> Result<BTreeMap<GraphNodeKey, Seq>, Report> {
  let mut sampled = BTreeMap::new();
  graph.iter_depth_first_preorder_forward(|node| {
    if !node.is_leaf && sample_mode.samples_node(node.is_root) {
      sampled.insert(node.key, sample(node.key));
    }
    Ok(())
  })?;
  Ok(sampled)
}

#[derive(Clone, Debug, Default, Serialize)]
pub(crate) struct ReconstructedSequences {
  pub sampled: BTreeMap<GraphNodeKey, Seq>,
  pub emitted_nodes: Vec<GraphNodeKey>,
}
