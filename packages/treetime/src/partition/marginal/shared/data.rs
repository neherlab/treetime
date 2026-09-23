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

pub(crate) fn count_transitions_dense(
  inputs: &DenseInputs,
  gtr: &GTR,
  graph: &Graph,
  branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
  node_states: &BTreeMap<GraphNodeKey, DenseNodeState>,
  backward: &BTreeMap<GraphEdgeKey, DenseEdgeBackward>,
  forward: &BTreeMap<GraphEdgeKey, DenseEdgeForward>,
) -> Result<MutationCounts, Report> {
  let n_states = gtr.pi.len();
  let mut nij = Array2::zeros((n_states, n_states));
  let mut Ti = Array1::zeros(n_states);

  for edge in graph.get_edges() {
    let edge_arc = edge;
    let edge_key = edge_arc.key();
    let branch_length = inputs.effective_branch_length(branch_lengths[&edge_key]);

    let msg_to_child = &forward[&edge_key].msg_to_child;
    let msg_to_parent = &backward[&edge_key].msg_to_parent;

    let exp_qt = gtr.expQt(branch_length) + SUPERTINY_NUMBER;
    let mut_stack = get_branch_mutation_matrix(&msg_to_child.dis, &msg_to_parent.dis, &exp_qt);
    accumulate_mutation_counts(&mut_stack, branch_length, &mut nij, &mut Ti);
  }

  let root = graph.get_exactly_one_root()?;
  let root_key = root.key();
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

#[derive(Clone, Debug, Serialize)]
pub struct DenseInputs {
  pub min_branch_length: f64,
  pub filter_uninformative_root: bool,
}

impl DenseInputs {
  fn effective_branch_length(&self, raw: f64) -> f64 {
    raw.max(self.min_branch_length)
  }
}
