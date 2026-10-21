use crate::alphabet::state_alphabet::StateAlphabet;
use crate::ancestral::attach::{LeafData, warn_leaves_without_data};
use crate::make_internal_report;
use crate::partition::marginal::discrete::partition::StateObservation;
use crate::progress::LogSink;
use eyre::Report;
use itertools::Itertools;
use std::collections::BTreeMap;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_graph::node_names::NodeNames;

pub(crate) fn leaf_observations(
  graph: &Graph,
  traits: &BTreeMap<GraphNodeKey, String>,
  names: &NodeNames,
  states: &StateAlphabet,
  missing_data: &str,
  source: &str,
  log: &dyn LogSink,
) -> Result<BTreeMap<GraphNodeKey, StateObservation>, Report> {
  let leaves = graph.get_leaves().map(|leaf| leaf.key()).collect_vec();
  let without_value = leaves
    .iter()
    .filter(|key| !traits.contains_key(key))
    .map(|key| names[*key].to_owned())
    .collect_vec();
  warn_leaves_without_data(&without_value, leaves.len(), LeafData::TraitValue, source, log);

  leaves
    .into_iter()
    .map(|key| match traits.get(&key) {
      Some(value) if value != missing_data => {
        let index = states.index_of(value).ok_or_else(|| {
          let name = &names[key];
          make_internal_report!("Trait value '{value}' of leaf '{name}' is not one of the trait states")
        })?;
        Ok((key, StateObservation::Known(index)))
      },
      _ => Ok((key, StateObservation::Unknown)),
    })
    .try_collect()
}
