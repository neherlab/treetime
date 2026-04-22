use crate::alphabet::alphabet::Alphabet;
use crate::commands::optimize::partition_ops::PartitionOptimizeOps;
use crate::commands::timetree::partition_ops::PartitionRerootOps;
use crate::gtr::gtr::GTR;
use crate::hacks::fix_branch_length::fix_branch_length;
use crate::make_report;
use crate::representation::partition::traits::BranchTopology;
use crate::representation::partition::traits::HasLogLh;
use crate::representation::partition::traits::PartitionBranchOps;
use crate::representation::partition::traits::{PartitionMarginal, PartitionMarginalOps};
use crate::representation::payload::dense::{DenseEdgePartition, DenseNodePartition, DenseSeqDis, DenseSeqInfo};
use crate::seq::mutation::Sub;
use eyre::Report;
use itertools::{Itertools, izip};
use ndarray::prelude::*;
use std::collections::BTreeMap;
use treetime_graph::edge::{EdgeOptimizeOps, GraphEdgeKey};
use treetime_graph::graph::Graph;
use treetime_graph::graph_traverse::{GraphNodeBackward, GraphNodeForward};
use treetime_graph::node::{GraphNode, GraphNodeKey, Named, NodeAncestralOps};
use treetime_io::fasta::FastaRecord;
use treetime_primitives::{Seq, seq};
use treetime_utils::array::ndarray::argmax_first;
use treetime_utils::array::softmax_with_log_norm::softmax_with_log_norm;
use treetime_utils::collections::container::get_exactly_one;
use treetime_utils::interval::range::range_contains;
use treetime_utils::interval::range_intersection::range_intersection;
use treetime_utils::interval::range_union::range_union;

#[derive(Clone, Debug)]
pub struct PartitionMarginalDense {
  pub index: usize,
  pub gtr: GTR,
  pub alphabet: Alphabet,
  pub length: usize,
  pub nodes: BTreeMap<GraphNodeKey, DenseNodePartition>,
  pub edges: BTreeMap<GraphEdgeKey, DenseEdgePartition>,
}

impl PartitionMarginalDense {
  #[allow(clippy::same_name_method)]
  pub fn get_sequence_length(&self) -> usize {
    self.length
  }
}

impl HasLogLh for PartitionMarginalDense {
  fn get_log_lh(&self, node_key: GraphNodeKey) -> f64 {
    self.nodes.get(&node_key).map_or(0.0, |node| node.profile.log_lh)
  }
}

impl PartitionMarginal for PartitionMarginalDense {}

impl PartitionRerootOps for PartitionMarginalDense {}

impl PartitionBranchOps for PartitionMarginalDense {
  fn sequence_length(&self) -> usize {
    self.length
  }

  /// Count branch mutations by comparing the MAP (argmax) states of the full
  /// marginal posteriors at the parent and child nodes.
  ///
  /// In dense marginal inference, each node's `profile.dis` is the full posterior
  /// distribution over states at each position, obtained by combining all incoming
  /// messages during the forward pass (backward: subtree evidence, forward: outgroup
  /// evidence). The MAP state at each position is the single most probable state
  /// given all available data.
  ///
  /// This function compares endpoint node posteriors, not per-edge partial messages.
  /// Per-edge messages (`msg_to_parent`, `msg_to_child`) represent one-sided evidence
  /// and their argmax can disagree with the full posterior argmax. For example, a
  /// uniform `msg_to_child` (uninformative outgroup) has argmax at index 0 by tie-
  /// breaking, but the child node posterior incorporates subtree evidence and may
  /// peak at a different state. Comparing partial-message argmax would report false
  /// mutations at such positions.
  ///
  /// `initial_guess_mixed()` calls this to count discrete substitutions for
  /// initial branch length estimation.
  fn edge_subs(&self, graph: &dyn BranchTopology, edge_key: GraphEdgeKey) -> Result<Vec<Sub>, Report> {
    let (parent_key, child_key) = graph.edge_endpoints(edge_key)?;
    let parent_gaps = &self.nodes[&parent_key].seq.gaps;
    let child_gaps = &self.nodes[&child_key].seq.gaps;

    let parent_profile = &self.nodes[&parent_key].profile.dis;
    let child_profile = &self.nodes[&child_key].profile.dis;
    let mut subs = Vec::new();

    for (pos, parent, child) in izip!(0..parent_profile.nrows(), parent_profile.rows(), child_profile.rows()) {
      // Gap positions get uniform profiles under treat_gap_as_unknown, so argmax
      // returns an arbitrary canonical state. Check original gap ranges instead.
      if range_contains(parent_gaps, pos) || range_contains(child_gaps, pos) {
        continue;
      }

      let parent_state = self.alphabet.char(argmax_first(&parent).unwrap_or(0));
      let child_state = self.alphabet.char(argmax_first(&child).unwrap_or(0));
      if parent_state == child_state {
        continue;
      }
      // Non-canonical states (ambiguity codes, unknowns) are not nucleotide
      // substitutions. This matches the sparse edge_subs() filter at
      // marginal_sparse.rs which skips non-canonical parent or child states.
      if !self.alphabet.is_canonical(parent_state) || !self.alphabet.is_canonical(child_state) {
        continue;
      }
      subs.push(Sub::new(parent_state, pos, child_state)?);
    }

    Ok(subs)
  }

  fn edge_effective_length(&self, graph: &dyn BranchTopology, edge_key: GraphEdgeKey) -> Result<usize, Report> {
    let (parent_key, child_key) = graph.edge_endpoints(edge_key)?;
    let parent_gaps = &self.nodes[&parent_key].seq.gaps;
    let child_gaps = &self.nodes[&child_key].seq.gaps;

    // Use DenseSeqInfo.gaps (propagated as intersection during Fitch backward
    // pass) because profile-based detection is unreliable: treat_gap_as_unknown
    // makes gap profiles indistinguishable from genuinely uncertain positions.
    let gap_positions: usize = range_union(&[parent_gaps.clone(), child_gaps.clone()])
      .iter()
      .map(|(start, end)| end - start)
      .sum();

    Ok(self.length.saturating_sub(gap_positions))
  }
}

impl PartitionOptimizeOps for PartitionMarginalDense {
  fn create_edge_contribution(
    &self,
    edge_key: GraphEdgeKey,
  ) -> Result<crate::commands::optimize::optimize_unified::OptimizationContribution, Report> {
    Ok(crate::commands::optimize::optimize_unified::OptimizationContribution::from_dense(edge_key, self))
  }

  fn edge_indel_count(&self, edge_key: GraphEdgeKey) -> usize {
    self.edges[&edge_key].indels.len()
  }
}

impl<N, E> crate::commands::timetree::partition_ops::PartitionTimetreeOps<N, E> for PartitionMarginalDense
where
  N: GraphNode + Named,
  E: EdgeOptimizeOps,
{
  fn reconcile_topology(&mut self, graph: &Graph<N, E, ()>) {
    let graph_node_keys: std::collections::BTreeSet<GraphNodeKey> =
      graph.get_nodes().into_iter().map(|n| n.read_arc().key()).collect();
    let graph_edge_keys: std::collections::BTreeSet<GraphEdgeKey> =
      graph.get_edges().into_iter().map(|e| e.read_arc().key()).collect();

    // Add missing nodes with default partition data (backward pass will recompute)
    for &key in &graph_node_keys {
      self.nodes.entry(key).or_insert_with(|| DenseNodePartition {
        seq: DenseSeqInfo::default(),
        profile: DenseSeqDis::default(),
      });
    }

    // Add missing edges with default partition data
    for &key in &graph_edge_keys {
      self.edges.entry(key).or_default();
    }

    // Remove stale entries for nodes/edges no longer in graph
    self.nodes.retain(|k, _| graph_node_keys.contains(k));
    self.edges.retain(|k, _| graph_edge_keys.contains(k));
  }
}

impl<N, E> PartitionMarginalOps<N, E> for PartitionMarginalDense
where
  N: NodeAncestralOps,
  E: EdgeOptimizeOps,
{
  fn attach_sequences(&mut self, graph: &Graph<N, E, ()>, aln: &[FastaRecord]) -> Result<(), Report> {
    for leaf in graph.get_leaves() {
      let leaf_key = leaf.read_arc().key();
      let mut leaf = leaf.read_arc().payload().write_arc();

      let leaf_name = {
        let name = leaf.name().ok_or_else(|| {
          make_report!("Expected all leaf nodes to have names, such that they can be matched to their corresponding sequences. But found a leaf node that has no name.")
        })?;
        name.as_ref().to_owned()
      };

      let leaf_fasta = aln
        .iter()
        .find(|fasta| fasta.seq_name == leaf_name)
        .ok_or_else(|| make_report!("Leaf sequence not found: '{leaf_name}'"))?;

      leaf.set_desc(leaf_fasta.desc.clone());

      let alphabet = &self.alphabet;
      self
        .nodes
        .insert(leaf_key, DenseNodePartition::new(&leaf_fasta.seq, alphabet)?);
    }

    for edge in graph.get_edges() {
      let edge_key = edge.read_arc().key();
      self.edges.insert(edge_key, DenseEdgePartition::default());
    }

    Ok(())
  }

  fn process_node_backward(&mut self, node: &GraphNodeBackward<N, E, ()>) -> Result<(), Report> {
    let alphabet = &self.alphabet;
    let length = self.length;

    let msg_to_parent = if node.is_leaf {
      let seq_info = &self.nodes[&node.key];
      DenseSeqDis {
        dis: alphabet.seq2prof(&seq_info.seq.sequence)?,
        log_lh: 0.0,
      }
    } else {
      let gaps = node
        .child_keys
        .iter()
        .map(|(child_key, _)| self.nodes[child_key].seq.gaps.clone()) // TODO: avoid cloning
        .collect_vec();

      let gaps = range_intersection(&gaps);

      let seq = DenseSeqInfo {
        gaps,
        ..DenseSeqInfo::default()
      };

      let node_data = DenseNodePartition {
        seq,
        profile: DenseSeqDis::default(),
      };
      self.nodes.insert(node.key, node_data);

      // Combine child messages in log-space to match v0's numerical behavior.
      // v0 uses: tmp_log_subtree_LH += ch.marginal_log_Lx (treeanc.py:874)
      // where marginal_log_Lx = np.log(profile.dot(Qt))
      // Convert to log and sum (equivalent to multiply in prob space)
      let mut child_iter = node.child_keys.iter();
      let (_, first_edge_key) = child_iter.next().expect("Internal node must have children");
      let mut log_dis = self.edges[first_edge_key].msg_from_child.dis.mapv(f64::ln);
      for (_, edge_key) in child_iter {
        log_dis += &self.edges[edge_key].msg_from_child.dis.mapv(f64::ln);
      }

      // Normalize using v0's normalize_profile(log=True) algorithm:
      // 1. Subtract row max (logsumexp trick)
      // 2. Exponentiate
      // 3. Normalize by row sum
      let (dis, delta_ll) = normalize_from_log(&log_dis);
      let log_lh = node
        .child_keys
        .iter()
        .map(|(_, edge_key)| self.edges[edge_key].msg_from_child.log_lh)
        .sum::<f64>();
      DenseSeqDis {
        dis,
        log_lh: log_lh + delta_ll,
      }
    };

    if node.is_root {
      let mut dis = &msg_to_parent.dis * &self.gtr.pi;
      let delta_ll = normalize_inplace(&mut dis);

      let seq_info = self.nodes.get_mut(&node.key).unwrap();
      seq_info.profile.dis = dis;
      seq_info.profile.log_lh = msg_to_parent.log_lh + delta_ll;
    } else {
      let edge_key = get_exactly_one(&node.parent_edge_keys).expect("Only nodes with exactly one parent are supported");
      let branch_length = node.parent_edges[0].branch_length().unwrap_or(0.0);
      let branch_length = fix_branch_length(length, branch_length);
      // Remove-mutate-reinsert to preserve fields set during Fitch reconstruction
      // (indels, transmission) that are not recomputed by the marginal backward pass.
      let mut edge_data = self.edges.remove(edge_key).unwrap();

      let log_lh = msg_to_parent.log_lh;
      let dis = self.gtr.propagate_profile(&msg_to_parent.dis, branch_length, false);
      edge_data.msg_from_child = DenseSeqDis { dis, log_lh };
      edge_data.msg_to_parent = msg_to_parent;
      self.edges.insert(*edge_key, edge_data);
    }
    Ok(())
  }

  fn process_node_forward(&mut self, graph: &Graph<N, E, ()>, node: &GraphNodeForward<N, E, ()>) -> Result<(), Report> {
    if !node.is_root {
      let mut dis: Option<Array2<f64>> = None;
      let mut log_lh = 0.0;
      for (_, edge_key) in &node.parent_keys {
        let edge = &self.edges[edge_key];
        let edge_payload = graph.get_edge(*edge_key).unwrap().read_arc().payload().read_arc();
        let branch_length = edge_payload.branch_length().unwrap_or(0.0);
        let msg_child = self.gtr.evolve(&edge.msg_to_child.dis, branch_length, false);
        log_lh += edge.msg_to_parent.log_lh + edge.msg_to_child.log_lh;
        dis = Some(match dis {
          None => &edge.msg_to_parent.dis * &msg_child,
          Some(mut dis) => {
            dis *= &edge.msg_to_parent.dis;
            dis *= &msg_child;
            dis
          },
        });
      }
      let mut dis = dis.expect("Non-root node must have at least one parent");
      let delta_ll = normalize_inplace(&mut dis);
      log_lh += delta_ll;
      self.nodes.get_mut(&node.key).unwrap().profile = DenseSeqDis { dis, log_lh };
    }

    let node_data = &self.nodes[&node.key];
    for child_edge_key in &node.child_edge_keys {
      let child_edge_data = self.edges.get_mut(child_edge_key).unwrap();

      // Guard against zero divisor: clamp to smallest positive normal f64.
      // Matches the stabilization pattern in discrete.rs:206.
      let safe_child = child_edge_data.msg_from_child.dis.mapv(|v| v.max(f64::MIN_POSITIVE));
      let mut dis = &node_data.profile.dis / &safe_child;
      let delta_ll = normalize_inplace(&mut dis);
      child_edge_data.msg_to_child = DenseSeqDis {
        dis,
        log_lh: node_data.profile.log_lh - child_edge_data.msg_from_child.log_lh + delta_ll,
      };
    }
    Ok(())
  }

  fn extract_ancestral_sequence(&self, node_key: GraphNodeKey) -> Seq {
    if let Some(seq_info) = self.nodes.get(&node_key) {
      assign_sequence(seq_info, &self.alphabet)
    } else {
      seq! {}
    }
  }

  fn reconstruct_node_sequence(&mut self, node: &GraphNodeForward<N, E, ()>, include_leaves: bool) -> Option<Seq> {
    if !include_leaves && node.is_leaf {
      return None;
    }

    let seq_info = self.nodes.get(&node.key)?;
    let seq = if node.is_leaf {
      seq_info.seq.sequence.clone()
    } else {
      assign_sequence(seq_info, &self.alphabet)
    };

    Some(seq)
  }

  fn get_sequence_length(&self) -> usize {
    self.length
  }
}

fn assign_sequence(seq_info: &DenseNodePartition, alphabet: &Alphabet) -> Seq {
  let mut seq = prof2seq(&seq_info.profile, alphabet);
  for gap in &seq_info.seq.gaps {
    seq[gap.0..gap.1].fill(alphabet.gap());
  }
  seq
}

/// Extract sequence from profile by taking argmax at each position.
///
/// Uses `argmax_first` for NumPy-compatible tie-breaking: when multiple states
/// have equal probability, returns the first (lowest index) state.
fn prof2seq(profile: &DenseSeqDis, alphabet: &Alphabet) -> Seq {
  let mut seq = seq! {};
  for row in profile.dis.rows() {
    let argmax = argmax_first(&row).unwrap_or(0);
    seq.push(alphabet.char(argmax));
  }
  seq
}

/// Normalize each row of a probability matrix to sum to 1 and return the total
/// log-likelihood (sum of ln(row_sum) across all rows).
///
/// When a row sums to zero or is non-finite, falls back to a uniform distribution
/// for that row and contributes NEG_INFINITY to the log-likelihood. This matches
/// the degenerate-row semantics of `softmax_with_log_norm`.
fn normalize_inplace(dis: &mut Array2<f64>) -> f64 {
  let norm = dis.sum_axis(Axis(1));
  let n_cols = dis.ncols() as f64;
  let mut log_lh = 0.0;
  for (ri, mut row) in dis.outer_iter_mut().enumerate() {
    let n = norm[ri];
    if n > 0.0 && n.is_finite() {
      row /= n;
      log_lh += n.ln();
    } else {
      row.fill(1.0 / n_cols);
      log_lh += f64::NEG_INFINITY;
    }
  }
  log_lh
}

/// Normalize a log-probability matrix to probability matrix.
///
/// Applies per-row logsumexp normalization and accumulates the total
/// log-likelihood across all rows.
fn normalize_from_log(log_dis: &Array2<f64>) -> (Array2<f64>, f64) {
  let mut dis = Array2::zeros(log_dis.raw_dim());
  let mut total_log_lh = 0.0;

  for (mut out_row, log_row) in izip!(dis.rows_mut(), log_dis.rows()) {
    let (normalized, log_norm) = softmax_with_log_norm(log_row);
    out_row.assign(&normalized);
    total_log_lh += log_norm;
  }

  (dis, total_log_lh)
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::pretty_assert_neg_inf;
  use approx::assert_abs_diff_eq;
  use ndarray::array;

  fn assert_valid_rows(dis: &Array2<f64>) {
    for row in dis.rows() {
      assert!(
        row.iter().all(|&v| v >= 0.0 && v.is_finite()),
        "all probabilities must be non-negative and finite: {row}"
      );
      assert_abs_diff_eq!(row.sum(), 1.0, epsilon = 1e-10);
    }
  }

  #[test]
  fn test_normalize_from_log_equal_probs() {
    let log_dis = Array2::from_elem((2, 4), 0.25_f64.ln());
    let (dis, log_lh) = normalize_from_log(&log_dis);

    assert_valid_rows(&dis);
    assert_abs_diff_eq!(dis, Array2::from_elem((2, 4), 0.25), epsilon = 1e-10);
    assert!(log_lh.is_finite());
  }

  #[test]
  fn test_normalize_from_log_descending() {
    let log_dis = array![[0.0, -1.0, -2.0, -3.0]];
    let (dis, log_lh) = normalize_from_log(&log_dis);

    let sum = 1.0 + (-1.0_f64).exp() + (-2.0_f64).exp() + (-3.0_f64).exp();
    let expected = array![[
      1.0 / sum,
      (-1.0_f64).exp() / sum,
      (-2.0_f64).exp() / sum,
      (-3.0_f64).exp() / sum
    ]];

    assert_valid_rows(&dis);
    assert_abs_diff_eq!(dis, expected, epsilon = 1e-10);
    assert_abs_diff_eq!(log_lh, sum.ln(), epsilon = 1e-10);
  }

  #[test]
  fn test_normalize_from_log_all_neg_inf_single_row() {
    let log_dis = Array2::from_elem((1, 4), f64::NEG_INFINITY);
    let (dis, log_lh) = normalize_from_log(&log_dis);

    assert_valid_rows(&dis);
    assert_abs_diff_eq!(dis, Array2::from_elem((1, 4), 0.25), epsilon = 1e-10);
    pretty_assert_neg_inf!(log_lh, "log_lh should be NEG_INFINITY, got {log_lh}");
  }

  #[test]
  fn test_normalize_from_log_all_neg_inf_multiple_rows() {
    let log_dis = Array2::from_elem((3, 4), f64::NEG_INFINITY);
    let (dis, log_lh) = normalize_from_log(&log_dis);

    assert_valid_rows(&dis);
    assert_abs_diff_eq!(dis, Array2::from_elem((3, 4), 0.25), epsilon = 1e-10);
    pretty_assert_neg_inf!(log_lh, "log_lh should be NEG_INFINITY, got {log_lh}");
  }

  #[test]
  fn test_normalize_from_log_mixed_neg_inf_and_finite_rows() {
    // Row 0: all -inf (degenerate), Row 1: normal values
    let mut log_dis = Array2::from_elem((2, 4), f64::NEG_INFINITY);
    log_dis[[1, 0]] = 0.0;
    log_dis[[1, 1]] = -1.0;
    log_dis[[1, 2]] = -2.0;
    log_dis[[1, 3]] = -3.0;

    let (dis, log_lh) = normalize_from_log(&log_dis);

    assert_valid_rows(&dis);
    // Row 0: uniform fallback
    assert_abs_diff_eq!(dis.row(0), array![0.25, 0.25, 0.25, 0.25], epsilon = 1e-10);
    // Row 1: normal normalization
    let sum = 1.0 + (-1.0_f64).exp() + (-2.0_f64).exp() + (-3.0_f64).exp();
    let expected_row1 = array![
      1.0 / sum,
      (-1.0_f64).exp() / sum,
      (-2.0_f64).exp() / sum,
      (-3.0_f64).exp() / sum
    ];
    assert_abs_diff_eq!(dis.row(1), expected_row1, epsilon = 1e-10);
    // total_log_lh is -inf because one row was degenerate
    pretty_assert_neg_inf!(log_lh, "log_lh should be NEG_INFINITY, got {log_lh}");
  }

  #[test]
  fn test_normalize_from_log_mixed_finite_neg_inf_within_row() {
    // Some states finite, others -inf within a row
    let log_dis = array![[0.0, f64::NEG_INFINITY, -1.0, f64::NEG_INFINITY]];
    let (dis, log_lh) = normalize_from_log(&log_dis);

    let sum = 1.0 + (-1.0_f64).exp();
    let expected = array![[1.0 / sum, 0.0, (-1.0_f64).exp() / sum, 0.0]];

    assert_valid_rows(&dis);
    assert_abs_diff_eq!(dis, expected, epsilon = 1e-10);
    assert_abs_diff_eq!(log_lh, sum.ln(), epsilon = 1e-10);
  }

  #[test]
  fn test_normalize_from_log_large_negative_values() {
    // Without logsumexp trick, exp(-1000) underflows to 0
    let log_dis = array![[-1000.0, -1001.0, -1002.0, -1003.0]];
    let (dis, log_lh) = normalize_from_log(&log_dis);

    let reference = array![[0.0, -1.0, -2.0, -3.0]];
    let (ref_dis, ref_log_lh) = normalize_from_log(&reference);

    assert_valid_rows(&dis);
    assert_abs_diff_eq!(dis, ref_dis, epsilon = 1e-10);
    assert_abs_diff_eq!(log_lh, ref_log_lh - 1000.0, epsilon = 1e-10);
  }

  #[test]
  fn test_normalize_from_log_three_states() {
    let log_dis = Array2::from_elem((1, 3), f64::NEG_INFINITY);
    let (dis, log_lh) = normalize_from_log(&log_dis);

    assert_valid_rows(&dis);
    assert_abs_diff_eq!(dis, Array2::from_elem((1, 3), 1.0 / 3.0), epsilon = 1e-10);
    pretty_assert_neg_inf!(log_lh, "log_lh should be NEG_INFINITY, got {log_lh}");
  }

  #[test]
  fn test_normalize_inplace_normal_rows() {
    let mut dis = array![[1.0, 2.0, 3.0, 4.0], [4.0, 3.0, 2.0, 1.0]];
    let log_lh = normalize_inplace(&mut dis);

    assert_valid_rows(&dis);
    assert_abs_diff_eq!(dis.row(0), array![0.1, 0.2, 0.3, 0.4], epsilon = 1e-10);
    assert_abs_diff_eq!(dis.row(1), array![0.4, 0.3, 0.2, 0.1], epsilon = 1e-10);
    assert_abs_diff_eq!(log_lh, 10.0_f64.ln() + 10.0_f64.ln(), epsilon = 1e-10);
  }

  #[test]
  fn test_normalize_inplace_zero_row_returns_uniform() {
    let mut dis = array![[0.0, 0.0, 0.0, 0.0]];
    let log_lh = normalize_inplace(&mut dis);

    assert_valid_rows(&dis);
    assert_abs_diff_eq!(dis.row(0), array![0.25, 0.25, 0.25, 0.25], epsilon = 1e-10);
    pretty_assert_neg_inf!(log_lh, "log_lh should be NEG_INFINITY for all-zero row, got {log_lh}");
  }

  #[test]
  fn test_normalize_inplace_mixed_zero_and_normal_rows() {
    let mut dis = array![[0.0, 0.0, 0.0, 0.0], [2.0, 4.0, 2.0, 2.0]];
    let log_lh = normalize_inplace(&mut dis);

    assert_valid_rows(&dis);
    // Zero row becomes uniform
    assert_abs_diff_eq!(dis.row(0), array![0.25, 0.25, 0.25, 0.25], epsilon = 1e-10);
    // Normal row is normalized
    assert_abs_diff_eq!(dis.row(1), array![0.2, 0.4, 0.2, 0.2], epsilon = 1e-10);
    // log_lh is NEG_INFINITY because one row was degenerate
    pretty_assert_neg_inf!(
      log_lh,
      "log_lh should be NEG_INFINITY when any row is zero, got {log_lh}"
    );
  }

  #[test]
  fn test_normalize_inplace_nan_row_returns_uniform() {
    let mut dis = array![[f64::NAN, f64::NAN, f64::NAN, f64::NAN]];
    let log_lh = normalize_inplace(&mut dis);

    assert_valid_rows(&dis);
    assert_abs_diff_eq!(dis.row(0), array![0.25, 0.25, 0.25, 0.25], epsilon = 1e-10);
    pretty_assert_neg_inf!(log_lh, "log_lh should be NEG_INFINITY for NaN row, got {log_lh}");
  }

  #[test]
  fn test_normalize_inplace_inf_row_returns_uniform() {
    let mut dis = array![[f64::INFINITY, 0.0, 0.0, 0.0]];
    let log_lh = normalize_inplace(&mut dis);

    assert_valid_rows(&dis);
    assert_abs_diff_eq!(dis.row(0), array![0.25, 0.25, 0.25, 0.25], epsilon = 1e-10);
    pretty_assert_neg_inf!(log_lh, "log_lh should be NEG_INFINITY for inf row, got {log_lh}");
  }
}
