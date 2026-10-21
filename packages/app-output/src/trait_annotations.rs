use crate::annotated_graph::TreeTraits;
use eyre::Report;
use itertools::izip;
use ndarray::Array1;
use std::collections::BTreeMap;
use treetime::alphabet::state_alphabet::StateAlphabet;
use treetime_graph::node::GraphNodeKey;

pub(crate) fn build_confidence_map(states: &StateAlphabet, profile: &Array1<f64>) -> BTreeMap<String, f64> {
  izip!(states.iter(), profile)
    .filter(|&(_, &probability)| probability > 0.001)
    .map(|(state, &probability)| (state.to_owned(), probability))
    .collect()
}

pub(crate) fn trait_change(
  traits: &TreeTraits<'_>,
  parent: GraphNodeKey,
  child: GraphNodeKey,
) -> Result<Option<String>, Report> {
  if traits.posteriors[&parent].state() == traits.posteriors[&child].state() {
    return Ok(None);
  }
  let parent = traits.state_name(parent)?;
  let child = traits.state_name(child)?;
  Ok(Some(format!("{parent} \u{2192} {child}")))
}
