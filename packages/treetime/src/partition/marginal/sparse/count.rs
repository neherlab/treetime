use crate::constants::{MIN_BRANCH_LENGTH_FRACTION, SUPERTINY_NUMBER};
use crate::gtr::gtr::GTR;
use crate::gtr::infer_gtr::common::{MutationCounts, is_profile_informative};
use crate::partition::storage::sparse::{
  SparseEdgeBackward, SparseEdgeForward, SparseNodeState, SparseSeqDistribution, VarPos,
};
use eyre::Report;
use ndarray::{Array1, Array2, aview1};
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_utils::array::ndarray::argmax_first;
use treetime_utils::cast::ToF64Exact;

pub(crate) fn count_transitions_sparse(
  gtr: &GTR,
  length: usize,
  graph: &Graph,
  branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
  node_states: &BTreeMap<GraphNodeKey, SparseNodeState>,
  backward: &BTreeMap<GraphEdgeKey, SparseEdgeBackward>,
  forward: &BTreeMap<GraphEdgeKey, SparseEdgeForward>,
) -> Result<MutationCounts, Report> {
  let n_states = gtr.pi().len();
  let min_bl = MIN_BRANCH_LENGTH_FRACTION / length.to_f64_exact();
  let mut nij = Array2::zeros((n_states, n_states));
  let mut Ti = Array1::zeros(n_states);

  for edge in graph.get_edges() {
    let edge_arc = edge;
    let edge_key = edge_arc.key();
    let branch_length = branch_lengths[&edge_key].max(min_bl);
    let msg_to_child = &forward[&edge_key].msg_to_child;
    let msg_to_parent = &backward[&edge_key].msg_to_parent;

    let exp_qt = gtr.expQt(branch_length) + SUPERTINY_NUMBER;

    accumulate_sparse_transitions(
      msg_to_child,
      msg_to_parent,
      &exp_qt,
      branch_length,
      n_states,
      &mut nij,
      &mut Ti,
    );
  }

  let root = graph.get_exactly_one_root()?;
  let root_key = root.key();
  let root_profile = &node_states[&root_key].profile;
  let mut root_state = Array1::zeros(n_states);
  let root_dis = aggregate_sparse_profile(root_profile, n_states);
  if is_profile_informative(root_dis.view(), n_states) {
    if let Some(root_idx) = argmax_first(&root_dis.view()) {
      root_state[root_idx] = 1.0;
    }
  }

  nij.diag_mut().fill(0.0);

  Ok(MutationCounts { nij, Ti, root_state })
}

fn accumulate_sparse_transitions(
  msg_to_child: &SparseSeqDistribution,
  msg_to_parent: &SparseSeqDistribution,
  exp_qt: &Array2<f64>,
  branch_length: f64,
  n_states: usize,
  nij: &mut Array2<f64>,
  Ti: &mut Array1<f64>,
) {
  for (pos, pp_var) in &msg_to_parent.variable {
    let pc = msg_to_child
      .variable
      .get(*pos)
      .map_or_else(|| fixed_profile_for_var(msg_to_child, pp_var), |v| v.dis);
    accumulate_site_transition_weighted(pp_var.dis, pc, exp_qt, branch_length, n_states, nij, Ti, 1);
  }

  for parent_fixed in msg_to_parent.fixed.iter() {
    let count = parent_fixed.count;
    if count == 0 {
      continue;
    }
    let child_fixed = msg_to_child
      .fixed
      .get(parent_fixed.state)
      .map_or(&parent_fixed.dis, |child_fixed| &child_fixed.dis);
    accumulate_site_transition_weighted(
      &parent_fixed.dis,
      child_fixed,
      exp_qt,
      branch_length,
      n_states,
      nij,
      Ti,
      count,
    );
  }
}

fn fixed_profile_for_var<'a>(dist: &'a SparseSeqDistribution, var: VarPos<'a>) -> &'a [f64] {
  dist.fixed.get(var.state).map_or(var.dis, |fixed| fixed.dis.as_slice())
}

#[expect(
  clippy::too_many_arguments,
  reason = "each argument is an independent input of this step; a parameter struct would be built only for this call"
)]
fn accumulate_site_transition_weighted(
  pp: &[f64],
  pc: &[f64],
  exp_qt: &Array2<f64>,
  branch_length: f64,
  n_states: usize,
  nij: &mut Array2<f64>,
  Ti: &mut Array1<f64>,
  weight: usize,
) {
  let weight = weight.to_f64_exact();
  let mut site_sum = 0.0;
  let mut joint = Array2::zeros((n_states, n_states));

  for i in 0..n_states {
    for j in 0..n_states {
      let val = pp[i] * exp_qt[[i, j]] * pc[j];
      joint[[i, j]] = val;
      site_sum += val;
    }
  }

  if site_sum > 0.0 {
    joint /= site_sum;
  }

  for i in 0..n_states {
    for j in 0..n_states {
      nij[[i, j]] += weight * joint[[i, j]];
    }
  }

  for k in 0..n_states {
    let mut parent_sum = 0.0;
    let mut child_sum = 0.0;
    for s in 0..n_states {
      parent_sum += joint[[s, k]];
      child_sum += joint[[k, s]];
    }
    Ti[k] += weight * 0.5 * branch_length * (parent_sum + child_sum);
  }
}

fn aggregate_sparse_profile(profile: &SparseSeqDistribution, n_states: usize) -> Array1<f64> {
  let mut result = Array1::zeros(n_states);
  for var in profile.variable.values() {
    if is_profile_informative(aview1(var.dis), n_states) {
      if let Some(idx) = argmax_first(&aview1(var.dis)) {
        result[idx] += 1.0;
      }
    }
  }
  for fixed in profile.fixed.iter() {
    let count = fixed.count.to_f64_exact();
    if count > 0.0 && is_profile_informative(aview1(&fixed.dis), n_states) {
      if let Some(idx) = argmax_first(&aview1(&fixed.dis)) {
        result[idx] += count;
      }
    }
  }
  result
}
