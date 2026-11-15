use crate::partition::marginal::sample::SampleMode;
use eyre::Report;
use treetime_graph::graph::Graph;
use treetime_graph::key_map::NodeMap;
use treetime_graph::node::GraphNodeKey;
use treetime_graph::schedule::GraphSchedule;
use treetime_graph::scope::Inheritance;
use treetime_primitives::Seq;

pub(crate) fn sample_internal_sequences(
  graph: &Graph,
  inheritance: &Inheritance,
  sample_mode: SampleMode,
  mut sample: impl FnMut(GraphNodeKey) -> Seq,
) -> Result<NodeMap<Seq>, Report> {
  let schedule = GraphSchedule::new(graph)?;
  let sampled = schedule
    .preorder()
    .iter()
    .filter(|&&key| {
      inheritance.children(&schedule, key).next().is_some()
        && sample_mode.samples_node(inheritance.parent(&schedule, key).is_none())
    })
    .map(|&key| (key, sample(key)))
    .collect();
  Ok(sampled)
}
