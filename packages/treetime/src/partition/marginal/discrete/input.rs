use crate::make_error;
use eyre::Report;
use indexmap::IndexSet;
use itertools::Itertools;
use ndarray::Array2;
use std::collections::BTreeMap;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;

pub(crate) fn one_hot_profile(index: usize, n_states: usize) -> Array2<f64> {
  let mut profile = Array2::zeros((1, n_states));
  profile[[0, index]] = 1.0;
  profile
}

pub(crate) fn missing_trait_profile(n_states: usize) -> Array2<f64> {
  Array2::ones((1, n_states))
}

pub(crate) fn validate_trait_leaves(
  graph: &Graph,
  traits: &BTreeMap<GraphNodeKey, String>,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
) -> Result<(), Report> {
  let missing_in_metadata: IndexSet<String> = graph
    .get_leaves()
    .map(|leaf| leaf.key())
    .filter(|key| !traits.contains_key(key))
    .map(|key| names[&key].clone().unwrap_or_default())
    .collect();
  if !missing_in_metadata.is_empty() {
    return make_error!(
      "Mugration: tree leaves missing from metadata: {}",
      missing_in_metadata.iter().join(", ")
    );
  }
  Ok(())
}
