use crate::alphabet::state_alphabet::StateAlphabet;
use crate::ancestral::attach::{LeafData, warn_leaves_without_data};
use crate::make_internal_report;
use crate::partition::marginal::discrete::partition::StateObservation;
use crate::progress::LogSink;
use eyre::Report;
use itertools::Itertools;
use std::collections::BTreeMap;
use treetime_graph::assign_node_names::node_name_or_key;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;

pub(crate) fn leaf_observations(
  graph: &Graph,
  traits: &BTreeMap<GraphNodeKey, String>,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  states: &StateAlphabet,
  missing_data: &str,
  source: &str,
  log: &dyn LogSink,
) -> Result<BTreeMap<GraphNodeKey, StateObservation>, Report> {
  let leaves = graph.get_leaves().map(|leaf| leaf.key()).collect_vec();
  let without_value = leaves
    .iter()
    .filter(|key| !traits.contains_key(key))
    .map(|key| node_name_or_key(*key, names[key].as_deref()))
    .collect_vec();
  warn_leaves_without_data(&without_value, leaves.len(), LeafData::TraitValue, source, log);

  leaves
    .into_iter()
    .map(|key| match traits.get(&key) {
      Some(value) if value != missing_data => {
        let index = states.index_of(value).ok_or_else(|| {
          let name = node_name_or_key(key, names[&key].as_deref());
          make_internal_report!("Trait value '{value}' of leaf '{name}' is not one of the trait states")
        })?;
        Ok((key, StateObservation::Known(index)))
      },
      _ => Ok((key, StateObservation::Unknown)),
    })
    .try_collect()
}
