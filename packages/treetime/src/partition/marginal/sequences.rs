use crate::partition::marginal::sample::SampleMode;
use eyre::Report;
use std::collections::BTreeMap;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_primitives::Seq;

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

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(crate) struct TipStates {
  pub(crate) include_leaves: bool,
  pub(crate) impute: bool,
}
