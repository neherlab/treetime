use crate::constants::{MIN_BRANCH_LENGTH_FRACTION, SUPERTINY_NUMBER};
use crate::gtr::infer_gtr::common::{MutationCounts, is_profile_informative};
use crate::partition::marginal::sparse::partition::PartitionMarginalSparse;
use crate::partition::storage::sparse::{SparseSeqDistribution, VarPos};
use crate::partition::traits::TransitionCounting;
use eyre::Report;
use ndarray::{Array1, Array2};
use treetime_graph::edge::EdgeOptimizeOps;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNode;
use treetime_utils::array::ndarray::argmax_first;

// Split by codepath stage: these helpers sit beside the trait impl that consumes them.
#[allow(clippy::multiple_inherent_impl)]
impl PartitionMarginalSparse {
  fn count_transitions_impl<N, E>(&self, graph: &Graph<N, E, ()>) -> Result<MutationCounts, Report>
  where
    N: GraphNode,
    E: EdgeOptimizeOps,
  {
    let n_states = self.gtr.pi.len();
    let min_bl = MIN_BRANCH_LENGTH_FRACTION / self.length as f64;
    let mut nij = Array2::zeros((n_states, n_states));
    let mut Ti = Array1::zeros(n_states);

    for edge in graph.get_edges() {
      let edge_arc = edge.read_arc();
      let branch_length = edge_arc.payload().read_arc().branch_length().unwrap_or(0.0).max(min_bl);
      let edge_key = edge_arc.key();
      let edge_data = &self.edges[&edge_key];

      let exp_qt = self.gtr.expQt(branch_length) + SUPERTINY_NUMBER;

      Self::accumulate_sparse_transitions(
        &edge_data.msg_to_child,
        &edge_data.msg_to_parent,
        &exp_qt,
        branch_length,
        n_states,
        &mut nij,
        &mut Ti,
      );
    }

    let root = graph.get_exactly_one_root()?;
    let root_key = root.read_arc().key();
    let root_profile = &self.nodes[&root_key].profile;
    let mut root_state = Array1::zeros(n_states);
    let root_dis = Self::aggregate_sparse_profile(root_profile, n_states);
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
        .get(pos)
        .map_or_else(|| Self::fixed_profile_for_var(msg_to_child, pp_var), |v| &v.dis);
      Self::accumulate_site_transition_weighted(&pp_var.dis, pc, exp_qt, branch_length, n_states, nij, Ti, 1);
    }

    for (ch, parent_fixed) in &msg_to_parent.fixed {
      let count = msg_to_parent.fixed_counts.get(*ch).unwrap_or(0);
      if count == 0 {
        continue;
      }
      let child_fixed = msg_to_child.fixed.get(ch).unwrap_or(parent_fixed);
      Self::accumulate_site_transition_weighted(
        parent_fixed,
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

  fn fixed_profile_for_var<'a>(dist: &'a SparseSeqDistribution, var: &'a VarPos) -> &'a Array1<f64> {
    dist.fixed.get(&var.state).unwrap_or(&var.dis)
  }

  fn accumulate_site_transition_weighted(
    pp: &Array1<f64>,
    pc: &Array1<f64>,
    exp_qt: &Array2<f64>,
    branch_length: f64,
    n_states: usize,
    nij: &mut Array2<f64>,
    Ti: &mut Array1<f64>,
    weight: usize,
  ) {
    let weight = weight as f64;
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
      if is_profile_informative(var.dis.view(), n_states) {
        if let Some(idx) = argmax_first(&var.dis.view()) {
          result[idx] += 1.0;
        }
      }
    }
    for (ch, fixed_profile) in &profile.fixed {
      let count = profile.fixed_counts.get(*ch).unwrap_or(0) as f64;
      if count > 0.0 && is_profile_informative(fixed_profile.view(), n_states) {
        if let Some(idx) = argmax_first(&fixed_profile.view()) {
          result[idx] += count;
        }
      }
    }
    result
  }
}

impl<N, E> TransitionCounting<N, E> for PartitionMarginalSparse
where
  N: GraphNode,
  E: EdgeOptimizeOps,
{
  fn count_transitions(&self, graph: &Graph<N, E, ()>) -> Result<MutationCounts, Report> {
    self.count_transitions_impl(graph)
  }
}
