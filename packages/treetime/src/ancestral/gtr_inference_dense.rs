use crate::constants::SUPERTINY_NUMBER;
use crate::gtr::gtr::{GTR, GTRParams};
use crate::gtr::infer_gtr::common::{
  InferGtrOptions, InferGtrResult, MutationCounts, infer_gtr_impl, is_profile_informative,
};
use crate::hacks::fix_branch_length::fix_branch_length;
use crate::make_internal_report;
use crate::partition::marginal_dense::PartitionMarginalDense;
use eyre::Report;
use ndarray::{Array1, Array2, Array3};
use parking_lot::RwLock;
use std::sync::Arc;
use treetime_graph::edge::{GraphEdge, HasBranchLength};
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNode;
use treetime_utils::array::ndarray::argmax_first;

/// Infer GTR model from dense partition data.
///
/// Requires profiles to be populated via `initialize_marginal` + `update_marginal` before calling.
/// When `--model=infer --dense=true`, marginal reconstruction runs twice: once to populate profiles
/// for GTR inference, and again after inferring the GTR. Sparse inference avoids this by tracking
/// mutations incrementally during Fitch compression.
///
/// Not wired into production since PR #623 replaced this path with Fitch-based GTR inference
/// via `PartitionFitch::infer_gtr`. Retained for potential future use as an alternative to the
/// Fitch approach when full marginal posterior counts are preferred over parsimony counts.
#[allow(dead_code)]
pub fn infer_gtr_dense<N, E, D>(
  partition: &Arc<RwLock<PartitionMarginalDense>>,
  graph: &Graph<N, E, D>,
) -> Result<GTR, Report>
where
  N: GraphNode,
  E: GraphEdge + HasBranchLength,
  D: Send + Sync,
{
  let counts = get_mutation_counts_dense(graph, partition)?;
  let InferGtrResult { W, pi, mu } = infer_gtr_impl(&counts, &InferGtrOptions::default())?;
  let n_states = partition.read_arc().alphabet.n_canonical();
  let W = Some(W);
  GTR::new(GTRParams { n_states, mu, W, pi })
}

/// Compute branch joint distribution for GTR inference.
///
/// Given edge messages and transition matrix, computes the posterior probability
/// P(child=i, parent=j | site a) for each site.
///
/// Arguments:
/// - `msg_to_child`: Outgroup likelihood at child, indexed by parent state (L x n)
/// - `msg_to_parent`: Subtree likelihood at child, indexed by child state (L x n)
/// - `exp_qt`: Transition matrix P(child=i | parent=j) (n x n)
///
/// Returns: Array3 of shape (L, n, n) where result[a, i, j] = P(child=i, parent=j | site a)
pub fn get_branch_mutation_matrix(
  msg_to_child: &Array2<f64>,
  msg_to_parent: &Array2<f64>,
  exp_qt: &Array2<f64>,
) -> Array3<f64> {
  let (n_sites, n_states) = msg_to_parent.dim();
  let mut result = Array3::zeros((n_sites, n_states, n_states));

  for a in 0..n_sites {
    let mut site_sum = 0.0;

    // Compute unnormalized joint: pc[a,i] * exp_qt[i,j] * pp[a,j]
    for i in 0..n_states {
      for j in 0..n_states {
        let val = msg_to_parent[[a, i]] * exp_qt[[i, j]] * msg_to_child[[a, j]];
        result[[a, i, j]] = val;
        site_sum += val;
      }
    }

    // Normalize to get posterior probabilities
    if site_sum > 0.0 {
      for i in 0..n_states {
        for j in 0..n_states {
          result[[a, i, j]] /= site_sum;
        }
      }
    }
  }

  result
}

/// Accumulate mutation counts and time-in-state from branch joint distribution.
///
/// Updates:
/// - `nij`: adds expected substitution counts summed over sites
/// - `Ti`: adds time spent in each state using midpoint approximation
///
/// The midpoint formula: Ti += 0.5 * branch_length * (P(parent=i) + P(child=i))
/// where marginals are computed by summing mut_stack over child/parent states.
pub fn accumulate_mutation_counts(
  mut_stack: &Array3<f64>,
  branch_length: f64,
  nij: &mut Array2<f64>,
  Ti: &mut Array1<f64>,
) {
  let (n_sites, n_states, _) = mut_stack.dim();

  // nij += sum over sites: mut_stack.sum(axis=0)
  for i in 0..n_states {
    for j in 0..n_states {
      let mut sum = 0.0;
      for a in 0..n_sites {
        sum += mut_stack[[a, i, j]];
      }
      nij[[i, j]] += sum;
    }
  }

  // Compute marginals and accumulate Ti
  // parent_marginal[a, j] = sum_i mut_stack[a, i, j] (sum over child states)
  // child_marginal[a, i] = sum_j mut_stack[a, i, j] (sum over parent states)
  // Ti[k] += 0.5 * branch_length * (sum_a parent_marginal[a,k] + sum_a child_marginal[a,k])
  for k in 0..n_states {
    let mut parent_sum = 0.0;
    let mut child_sum = 0.0;
    for a in 0..n_sites {
      // parent_marginal: sum over child states i
      for i in 0..n_states {
        parent_sum += mut_stack[[a, i, k]];
      }
      // child_marginal: sum over parent states j
      for j in 0..n_states {
        child_sum += mut_stack[[a, k, j]];
      }
    }
    Ti[k] += 0.5 * branch_length * (parent_sum + child_sum);
  }
}

/// Get mutation counts from dense partition for GTR inference.
///
/// Uses fractional expected counts from branch joint distributions:
/// - `nij`: expected substitutions weighted by posterior P(child=i, parent=j | site)
/// - `Ti`: time in state using midpoint approximation
/// - `root_state`: state counts from argmax sequence (matching v0 `cseq` counts)
///
/// See `infer_gtr_dense` for production status.
#[allow(dead_code)]
pub fn get_mutation_counts_dense<N, E, D>(
  graph: &Graph<N, E, D>,
  partition: &Arc<RwLock<PartitionMarginalDense>>,
) -> Result<MutationCounts, Report>
where
  N: GraphNode,
  E: GraphEdge + HasBranchLength,
  D: Send + Sync,
{
  let partition = partition.read_arc();
  let alphabet = &partition.alphabet;
  let n_states = alphabet.n_canonical();

  // root_state: count of positions in each state at the tree root, used as a prior
  // on equilibrium frequencies (pi) in the iterative GTR solver. Positions where the
  // marginal profile is near-uniform (max <= 1/n + eps) carry no phylogenetic signal
  // and are excluded to avoid biasing pi toward an arbitrary argmax state.
  let root_state = {
    let root = graph.get_exactly_one_root()?;
    let root_key = root.read_arc().key();
    let root_profile = &partition.data.nodes[&root_key].profile.dis;
    let mut counts = Array1::zeros(n_states);
    for row in root_profile.rows() {
      if is_profile_informative(row, n_states) {
        let state =
          argmax_first(&row).ok_or_else(|| make_internal_report!("Empty profile row in root marginal distribution"))?;
        counts[state] += 1.0;
      }
    }
    counts
  };

  // TODO(dense-branch-clamp): clamping here matches v0's `_branch_length_to_gtr()` in
  // `infer_gtr` (`treeanc.py:1562-1589`) and is consistent with the marginal passes that
  // computed the profiles this function reads (`marginal_dense.rs:198`, `marginal_passes.rs:105`).
  //
  // This is a consistency fix, not a scientific requirement. For zero-length branches the
  // mathematically correct result is `expQt(0) = I`, `Ti = 0`, zero off-diagonal `nij`.
  // Clamping exists because the marginal backward pass computes `log(expQt)` which diverges
  // at zero. The counting code here operates in probability space where zero BL is not
  // numerically problematic on its own - but using `expQt(0)` when profiles were propagated
  // through `expQt(clamped_BL)` mixes two different transition models on the same branch.
  //
  // Practical impact: none for real datasets (branch lengths are always positive). The
  // clamped minimum (`1e-3 / seq_length`) is negligibly small. All golden master tests
  // produce identical output with or without this clamping.
  //
  // If `fix_branch_length` is removed from the marginal passes, remove it here too.
  let length = partition.length;
  let mut nij = Array2::zeros((n_states, n_states));
  let mut Ti = Array1::zeros(n_states);

  for edge in graph.get_edges() {
    let edge_arc = edge.read_arc();
    let branch_length = fix_branch_length(length, edge_arc.payload().read_arc().branch_length().unwrap_or(0.0));
    let edge_key = edge_arc.key();

    let edge_partition = &partition.data.edges[&edge_key];
    let msg_to_child = &edge_partition.msg_to_child.dis;
    let msg_to_parent = &edge_partition.msg_to_parent.dis;

    let exp_qt = partition.data.gtr.expQt(branch_length) + SUPERTINY_NUMBER;
    let mut_stack = get_branch_mutation_matrix(msg_to_child, msg_to_parent, &exp_qt);
    accumulate_mutation_counts(&mut_stack, branch_length, &mut nij, &mut Ti);
  }

  // Zero diagonal - diagonal represents no-change events, not mutations
  for i in 0..n_states {
    nij[[i, i]] = 0.0;
  }

  Ok(MutationCounts { nij, Ti, root_state })
}
