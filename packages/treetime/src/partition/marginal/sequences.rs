use crate::partition::marginal::sample::SampleMode;
use eyre::Report;
use treetime_graph::graph::Graph;
use treetime_graph::key_map::NodeMap;
use treetime_graph::node::GraphNodeKey;
use treetime_graph::schedule::GraphSchedule;
use treetime_primitives::Seq;

pub(crate) fn sample_internal_sequences(
  graph: &Graph,
  sample_mode: SampleMode,
  mut sample: impl FnMut(GraphNodeKey) -> Seq,
) -> Result<NodeMap<Seq>, Report> {
  let schedule = GraphSchedule::new(graph)?;
  let sampled = schedule
    .preorder()
    .iter()
    .filter(|&&key| !schedule.children(key).is_empty() && sample_mode.samples_node(schedule.parents(key).is_empty()))
    .map(|&key| (key, sample(key)))
    .collect();
  Ok(sampled)
}
