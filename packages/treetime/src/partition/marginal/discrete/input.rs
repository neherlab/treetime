use crate::make_error;
use eyre::Report;
use indexmap::IndexSet;
use itertools::Itertools;
use log::warn;
use ndarray::Array2;
use std::collections::BTreeMap;
use treetime_graph::edge::EdgeOptimizeOps;
use treetime_graph::graph::Graph;
use treetime_graph::node::{GraphNode, Named};

pub(crate) fn one_hot_profile(index: usize, n_states: usize) -> Array2<f64> {
  let mut profile = Array2::zeros((1, n_states));
  profile[[0, index]] = 1.0;
  profile
}

pub(crate) fn uniform_profile(n_states: usize) -> Array2<f64> {
  Array2::from_elem((1, n_states), 1.0 / n_states as f64)
}

pub(crate) fn validate_trait_names<N, E>(
  graph: &Graph<N, E, ()>,
  traits: &BTreeMap<String, String>,
) -> Result<(), Report>
where
  N: GraphNode + Named,
  E: EdgeOptimizeOps,
{
  let leaf_names: IndexSet<String> = graph
    .get_leaves()
    .iter()
    .map(|leaf| {
      let leaf = leaf.read_arc();
      let payload = leaf.payload().read_arc();
      payload.name().map(|name| name.as_ref().to_owned()).unwrap_or_default()
    })
    .collect();
  let trait_names: IndexSet<String> = traits.keys().cloned().collect();

  let missing_in_metadata: IndexSet<String> = leaf_names.difference(&trait_names).cloned().collect();
  if !missing_in_metadata.is_empty() {
    return make_error!(
      "Mugration: tree leaves missing from metadata: {}",
      missing_in_metadata.iter().join(", ")
    );
  }

  // Metadata naming samples absent from the tree is the common case: metadata files are shared
  // across analyses and routinely list more samples than a pruned or subsampled tree. Warn for
  // visibility (matching the dates subsystem) instead of rejecting the run.
  let missing_in_tree: IndexSet<String> = trait_names.difference(&leaf_names).cloned().collect();
  if !missing_in_tree.is_empty() {
    let sample = missing_in_tree.iter().take(10).join(", ");
    let suffix = if missing_in_tree.len() > 10 { "..." } else { "" };
    warn!(
      "Mugration: {} metadata names not present in tree: {sample}{suffix}",
      missing_in_tree.len()
    );
  }

  Ok(())
}
