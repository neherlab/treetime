use crate::gtr::gtr::{GTR, GTRParams};
use crate::gtr::infer_gtr::common::{InferGtrOptions, InferGtrResult, MutationCounts, infer_gtr_impl};
use crate::representation::partition::marginal_dense::PartitionMarginalDense;
use crate::representation::payload::ancestral::GraphAncestral;
use eyre::Report;
use ndarray::{Array1, Array2, Array3};
use parking_lot::RwLock;
use std::sync::Arc;
use treetime_graph::edge::HasBranchLength;
use crate::make_internal_report;
use treetime_utils::array::ndarray::argmax_first;

/// Infer GTR model from dense partition data.
///
/// Requires profiles to be populated via `initialize_marginal` + `update_marginal` before calling.
/// When `--model=infer --dense=true`, marginal reconstruction runs twice: once to populate profiles
/// for GTR inference, and again after inferring the GTR. Sparse inference avoids this by tracking
/// mutations incrementally during Fitch compression.
pub fn infer_gtr_dense(partition: &Arc<RwLock<PartitionMarginalDense>>, graph: &GraphAncestral) -> Result<GTR, Report> {
  let counts = get_mutation_counts_dense(graph, partition)?;
  let InferGtrResult { W, pi, mu } = infer_gtr_impl(&counts, &InferGtrOptions::default())?;
  let alphabet = partition.read_arc().alphabet.clone();
  let W = Some(W);
  GTR::new(GTRParams { alphabet, mu, W, pi })
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
fn get_branch_mutation_matrix(
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
fn accumulate_mutation_counts(
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
pub fn get_mutation_counts_dense(
  graph: &GraphAncestral,
  partition: &Arc<RwLock<PartitionMarginalDense>>,
) -> Result<MutationCounts, Report> {
  let partition = partition.read_arc();
  let alphabet = &partition.alphabet;
  let n_states = alphabet.n_canonical();

  // TODO(dense-root-state): root_state diverges from v0 on alignments with many gap-only positions.
  //
  // root_state feeds into pi (equilibrium frequencies) via infer_gtr_impl. Both v0 and v1
  // use argmax of the marginal profile per site, but at uninformative positions (all leaves
  // are gap/N), the profile is uniform [1/n, 1/n, ..., 1/n]. Argmax on uniform rows is
  // arbitrary - tiny floating-point differences between BLAS implementations (NumPy vs
  // ndarray) cause v0 and v1 to pick different states.
  //
  // Impact: lassa_L_50 has 557 gap-only positions out of ~10k. The 557 counts shift
  // between states (e.g. all to A in v1 vs all to T in v0), causing ~29% error in W.
  // tb_20 and mpox_clade_ii_20 show the same pattern at smaller scale.
  //
  // Possible fixes:
  //
  // A. Skip uninformative positions: exclude rows where max(profile) <= 1/n + eps from
  //    root_state counting. These positions carry no phylogenetic signal and should not
  //    bias equilibrium frequency estimates. Requires matching change in capture script.
  //
  // B. Use expected counts: replace argmax with profile.sum_axis(Axis(0)), accumulating
  //    fractional posterior mass per state. A uniform row contributes [1/n, ..., 1/n]
  //    identically in both v0 and v1 (profile sums already match). Changes root_state
  //    semantics from "most likely state counts" to "posterior probability mass".
  //
  // C. Use Fitch consensus sequence counts (like sparse variant does via composition).
  //    Fitch reconstruction resolves ties deterministically, but dense partition does not
  //    store Fitch results after marginal initialization.
  //
  // Option A is preferred: scientifically correct (uninformative sites should not
  // influence equilibrium frequencies) and minimal code change.
  let root_state = {
    let root = graph.get_exactly_one_root()?;
    let root_key = root.read_arc().key();
    let root_profile = &partition.nodes[&root_key].profile.dis;
    let mut counts = Array1::zeros(n_states);
    for row in root_profile.rows() {
      let state =
        argmax_first(&row).ok_or_else(|| make_internal_report!("Empty profile row in root marginal distribution"))?;
      counts[state] += 1.0;
    }
    counts
  };

  let mut nij = Array2::zeros((n_states, n_states));
  let mut Ti = Array1::zeros(n_states);

  for edge in graph.get_edges() {
    let edge_arc = edge.read_arc();
    let branch_length = edge_arc.payload().read_arc().branch_length().unwrap_or(0.0);
    let edge_key = edge_arc.key();

    let edge_partition = &partition.edges[&edge_key];
    let msg_to_child = &edge_partition.msg_to_child.dis;
    let msg_to_parent = &edge_partition.msg_to_parent.dis;

    let exp_qt = partition.gtr.expQt(branch_length);
    let mut_stack = get_branch_mutation_matrix(msg_to_child, msg_to_parent, &exp_qt);
    accumulate_mutation_counts(&mut_stack, branch_length, &mut nij, &mut Ti);
  }

  // Zero diagonal - diagonal represents no-change events, not mutations
  for i in 0..n_states {
    nij[[i, i]] = 0.0;
  }

  Ok(MutationCounts { nij, Ti, root_state })
}
