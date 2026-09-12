use crate::constants::SUPERTINY_NUMBER;
use crate::gtr::gtr::GTR;
use crate::gtr::infer_gtr::common::{
  MutationCounts, accumulate_mutation_counts, get_branch_mutation_matrix, is_profile_informative,
};
use crate::partition::storage::dense::{DenseEdgeBackward, DenseEdgeForward, DenseNodeState};
use eyre::Report;
use ndarray::prelude::*;
use serde::Serialize;
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_utils::array::ndarray::argmax_first;

/// Durable model inputs shared by the dense and discrete marginal representations: the substitution
/// model and the numerical guards. The partition owns these; the stage-filled node and edge result
/// maps are owned separately by the passes' returned values.
#[derive(Clone, Debug, Serialize)]
pub struct DenseInputs {
  pub gtr: GTR,
  pub min_branch_length: f64,
  /// When `true`, root positions whose posterior profile is essentially uniform
  /// are excluded from the equilibrium-frequency prior in `count_transitions`.
  /// v0 never filters (always folds in the root MAP state); the nucleotide
  /// ancestral path enables filtering to drop signal-free gap-only columns.
  pub filter_uninformative_root: bool,
}

impl DenseInputs {
  pub fn effective_branch_length(&self, raw: f64) -> f64 {
    raw.max(self.min_branch_length)
  }
}

/// Count posterior-weighted transitions from dense profile matrices, reading the backward and forward
/// edge messages and the node states by their distinct owners.
///
/// Shared by dense and discrete partitions (both store full profile matrices).
pub fn count_transitions_dense(
  inputs: &DenseInputs,
  graph: &Graph,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  node_states: &BTreeMap<GraphNodeKey, DenseNodeState>,
  backward: &BTreeMap<GraphEdgeKey, DenseEdgeBackward>,
  forward: &BTreeMap<GraphEdgeKey, DenseEdgeForward>,
) -> Result<MutationCounts, Report> {
  let n_states = inputs.gtr.pi.len();
  let mut nij = Array2::zeros((n_states, n_states));
  let mut Ti = Array1::zeros(n_states);

  for edge in graph.get_edges() {
    let edge_arc = edge.read_arc();
    let edge_key = edge_arc.key();
    let branch_length = inputs.effective_branch_length(branch_lengths[&edge_key].unwrap_or(0.0));

    let msg_to_child = &forward[&edge_key].msg_to_child;
    let msg_to_parent = &backward[&edge_key].msg_to_parent;

    let exp_qt = inputs.gtr.expQt(branch_length) + SUPERTINY_NUMBER;
    let mut_stack = get_branch_mutation_matrix(&msg_to_child.dis, &msg_to_parent.dis, &exp_qt);
    accumulate_mutation_counts(&mut_stack, branch_length, &mut nij, &mut Ti);
  }

  let root = graph.get_exactly_one_root()?;
  let root_key = root.read_arc().key();
  let root_profile = &node_states[&root_key].profile.dis;
  let mut root_state = Array1::zeros(n_states);
  for row in root_profile.rows() {
    if !inputs.filter_uninformative_root || is_profile_informative(row, n_states) {
      if let Some(root_idx) = argmax_first(&row) {
        root_state[root_idx] += 1.0;
      }
    }
  }

  nij.diag_mut().fill(0.0);

  Ok(MutationCounts { nij, Ti, root_state })
}
