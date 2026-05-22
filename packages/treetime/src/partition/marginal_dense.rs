use crate::alphabet::alphabet::Alphabet;
use crate::ancestral::fitch_indel::{compute_node_ranges, resolve_indels_backward, resolve_indels_forward};
use crate::constants::MIN_BRANCH_LENGTH_FRACTION;
use crate::gtr::gtr::GTR;
use crate::gtr::infer_gtr::common::MutationCounts;
use crate::make_report;
use crate::partition::dense::{DenseEdgePartition, DenseNodePartition, DenseSeqDistribution, DenseSeqInfo};
use crate::partition::marginal_core::{
  MarginalData, MarginalPartition, count_transitions_from_marginal_data, marginal_process_node_backward,
  marginal_process_node_forward,
};
use crate::partition::optimization_contribution::OptimizationContribution;
use crate::partition::traits::{
  BranchTopology, HasGtr, HasLogLh, PartitionBranchOps, PartitionMarginalOps, PartitionMarginalPasses,
  PartitionOptimizeOps, PartitionRerootOps, PartitionTimetreeOps, TransitionCounting,
};
use crate::seq::mutation::Sub;
use eyre::Report;
use itertools::{Itertools, izip};
use maplit::btreemap;
use std::collections::{BTreeMap, BTreeSet};
use treetime_graph::edge::{EdgeOptimizeOps, GraphEdgeKey};
use treetime_graph::graph::Graph;
use treetime_graph::graph_traverse::{GraphNodeBackward, GraphNodeForward};
use treetime_graph::node::{GraphNode, GraphNodeKey, Named, NodeAncestralOps};
use treetime_io::fasta::FastaRecord;
use treetime_primitives::{Seq, seq};
use treetime_utils::array::ndarray::argmax_first;
use treetime_utils::collections::container::get_exactly_one;
use treetime_utils::interval::range::range_contains;
use treetime_utils::interval::range_union::range_union;

#[derive(Clone, Debug)]
pub struct PartitionMarginalDense {
  pub data: MarginalData,
  pub index: usize,
  pub alphabet: Alphabet,
  pub length: usize,
}

impl PartitionMarginalDense {
  pub fn new(index: usize, gtr: GTR, alphabet: Alphabet, length: usize) -> Self {
    let min_branch_length = MIN_BRANCH_LENGTH_FRACTION / length as f64;
    Self {
      data: MarginalData {
        gtr,
        nodes: btreemap! {},
        edges: btreemap! {},
        min_branch_length,
      },
      index,
      alphabet,
      length,
    }
  }

  #[allow(clippy::same_name_method)]
  pub fn get_sequence_length(&self) -> usize {
    self.length
  }
}

impl HasGtr for PartitionMarginalDense {
  fn gtr(&self) -> &GTR {
    &self.data.gtr
  }
  fn gtr_mut(&mut self) -> &mut GTR {
    &mut self.data.gtr
  }
  fn sequence_length(&self) -> usize {
    self.length
  }
}

impl HasLogLh for PartitionMarginalDense {
  fn get_log_lh(&self, node_key: GraphNodeKey) -> f64 {
    self.data.nodes.get(&node_key).map_or(0.0, |node| node.profile.log_lh)
  }

  fn reset_node_log_likelihoods(&mut self) {
    for node_data in self.data.nodes.values_mut() {
      node_data.profile.log_lh = 0.0;
    }
  }
}

impl PartitionRerootOps for PartitionMarginalDense {}

impl<N, E> TransitionCounting<N, E> for PartitionMarginalDense
where
  N: GraphNode,
  E: EdgeOptimizeOps,
{
  fn count_transitions(&self, graph: &Graph<N, E, ()>) -> Result<MutationCounts, Report> {
    count_transitions_from_marginal_data(&self.data, graph)
  }
}

impl PartitionBranchOps for PartitionMarginalDense {
  fn sequence_length(&self) -> usize {
    self.length
  }

  fn edge_subs(&self, graph: &dyn BranchTopology, edge_key: GraphEdgeKey) -> Result<Vec<Sub>, Report> {
    let (parent_key, child_key) = graph.edge_endpoints(edge_key)?;
    let parent_non_char = &self.data.nodes[&parent_key].seq.non_char;
    let child_non_char = &self.data.nodes[&child_key].seq.non_char;

    let parent_profile = &self.data.nodes[&parent_key].profile.dis;
    let child_profile = &self.data.nodes[&child_key].profile.dis;
    let mut subs = Vec::new();

    for (pos, parent, child) in izip!(0..parent_profile.nrows(), parent_profile.rows(), child_profile.rows()) {
      let parent_state = self.alphabet.char(argmax_first(&parent).unwrap_or(0));
      let child_state = self.alphabet.char(argmax_first(&child).unwrap_or(0));
      if parent_state == child_state {
        continue;
      }
      if !self.alphabet.is_canonical(parent_state) || !self.alphabet.is_canonical(child_state) {
        continue;
      }
      if range_contains(parent_non_char, pos) || range_contains(child_non_char, pos) {
        continue;
      }

      subs.push(Sub::new(parent_state, pos, child_state)?);
    }

    Ok(subs)
  }

  fn edge_effective_length(&self, graph: &dyn BranchTopology, edge_key: GraphEdgeKey) -> Result<usize, Report> {
    let (parent_key, child_key) = graph.edge_endpoints(edge_key)?;
    let parent_non_char = &self.data.nodes[&parent_key].seq.non_char;
    let child_non_char = &self.data.nodes[&child_key].seq.non_char;

    let non_char_positions: usize = range_union(&[parent_non_char.clone(), child_non_char.clone()])
      .iter()
      .map(|(start, end)| end - start)
      .sum();

    Ok(self.length.saturating_sub(non_char_positions))
  }
}

impl PartitionOptimizeOps for PartitionMarginalDense {
  fn create_edge_contribution(&self, edge_key: GraphEdgeKey) -> Result<OptimizationContribution, Report> {
    Ok(OptimizationContribution::from_dense(edge_key, self))
  }

  fn edge_indel_count(&self, edge_key: GraphEdgeKey) -> usize {
    self.data.edges[&edge_key].indels.len()
  }
}

impl<N, E> PartitionTimetreeOps<N, E> for PartitionMarginalDense
where
  N: GraphNode + Named,
  E: EdgeOptimizeOps,
{
  fn reconcile_topology(&mut self, graph: &Graph<N, E, ()>) {
    let graph_node_keys: BTreeSet<GraphNodeKey> = graph.get_nodes().into_iter().map(|n| n.read_arc().key()).collect();
    let graph_edge_keys: BTreeSet<GraphEdgeKey> = graph.get_edges().into_iter().map(|e| e.read_arc().key()).collect();

    for &key in &graph_node_keys {
      self.data.nodes.entry(key).or_insert_with(|| DenseNodePartition {
        seq: DenseSeqInfo::default(),
        profile: DenseSeqDistribution::default(),
      });
    }

    for &key in &graph_edge_keys {
      self.data.edges.entry(key).or_default();
    }

    self.data.nodes.retain(|k, _| graph_node_keys.contains(k));
    self.data.edges.retain(|k, _| graph_edge_keys.contains(k));
  }
}

impl<N, E> MarginalPartition<N, E> for PartitionMarginalDense
where
  N: NodeAncestralOps,
  E: EdgeOptimizeOps,
{
  fn marginal_data(&self) -> &MarginalData {
    &self.data
  }

  fn marginal_data_mut(&mut self) -> &mut MarginalData {
    &mut self.data
  }

  fn leaf_profile(&self, node_key: GraphNodeKey) -> Result<DenseSeqDistribution, Report> {
    let seq_info = &self.data.nodes[&node_key];
    Ok(DenseSeqDistribution {
      dis: self.alphabet.seq2prof(&seq_info.seq.sequence)?,
      log_lh: 0.0,
    })
  }

  fn backward_internal_pre(&mut self, node: &GraphNodeBackward<N, E, ()>) {
    let child_non_chars: Vec<&Vec<(usize, usize)>> = node
      .child_keys
      .iter()
      .map(|(child_key, _)| &self.data.nodes[child_key].seq.non_char)
      .collect_vec();
    let child_gaps: Vec<&Vec<(usize, usize)>> = node
      .child_keys
      .iter()
      .map(|(child_key, _)| &self.data.nodes[child_key].seq.gaps)
      .collect_vec();

    let ranges = compute_node_ranges(&child_non_chars, &child_gaps);
    let non_char = ranges.non_char;
    let unknown = ranges.unknown;

    let child_unknown: Vec<&Vec<(usize, usize)>> = node
      .child_keys
      .iter()
      .map(|(child_key, _)| &self.data.nodes[child_key].seq.unknown)
      .collect_vec();
    let child_variable_indels: Vec<&BTreeMap<(usize, usize), _>> = node
      .child_keys
      .iter()
      .map(|(child_key, _)| &self.data.nodes[child_key].seq.variable_indel)
      .collect_vec();

    let indels_bw = resolve_indels_backward(&child_gaps, &child_unknown, &child_variable_indels, self.length);

    let seq = DenseSeqInfo {
      gaps: indels_bw.resolved_gaps,
      unknown,
      non_char,
      variable_indel: indels_bw.variable_indel,
      ..DenseSeqInfo::default()
    };

    self.data.nodes.insert(
      node.key,
      DenseNodePartition {
        seq,
        profile: DenseSeqDistribution::default(),
      },
    );
  }

  fn forward_post(&mut self, graph: &Graph<N, E, ()>, node: &GraphNodeForward<N, E, ()>) -> Result<(), Report> {
    if node.is_root {
      let node_data = self.data.nodes.get_mut(&node.key).unwrap();
      for (r, indel) in &node_data.seq.variable_indel {
        if indel.deleted > indel.present {
          node_data.seq.gaps.push(*r);
        }
      }
      node_data.seq.variable_indel.clear();
      let seq = assign_sequence(node_data, &self.alphabet);
      node_data.seq.sequence = seq;
      node_data.seq.non_char = range_union(&[node_data.seq.gaps.clone(), node_data.seq.unknown.clone()]);
    } else {
      if !node.is_leaf {
        let node_data = self.data.nodes.get_mut(&node.key).unwrap();
        node_data.seq.sequence = assign_sequence(node_data, &self.alphabet);
      }

      let (parent_key, edge_key) =
        get_exactly_one(&node.parent_keys).expect("Non-root dense node must have exactly one parent");
      let parent_gaps = self.data.nodes[parent_key].seq.gaps.clone();
      let parent_sequence = self.data.nodes[parent_key].seq.sequence.clone();

      let node_data = self.data.nodes.get_mut(&node.key).unwrap();
      let variable_indel = std::mem::take(&mut node_data.seq.variable_indel);

      let (indels, new_gaps) = resolve_indels_forward(
        &variable_indel,
        &node_data.seq.gaps,
        &node_data.seq.non_char,
        &parent_gaps,
        &parent_sequence,
        &node_data.seq.sequence,
      );
      node_data.seq.gaps = new_gaps;
      node_data.seq.non_char = range_union(&[node_data.seq.gaps.clone(), node_data.seq.unknown.clone()]);
      for gap in &node_data.seq.gaps {
        node_data.seq.sequence[gap.0..gap.1].fill(self.alphabet.gap());
      }
      for unk in &node_data.seq.unknown {
        node_data.seq.sequence[unk.0..unk.1].fill(self.alphabet.unknown());
      }

      let edge_data = self.data.edges.get_mut(edge_key).unwrap();
      edge_data.indels = indels;
    }
    Ok(())
  }
}

impl<N, E> PartitionMarginalPasses<N, E> for PartitionMarginalDense
where
  N: NodeAncestralOps,
  E: EdgeOptimizeOps,
{
  fn process_node_backward(&mut self, node: &GraphNodeBackward<N, E, ()>) -> Result<(), Report> {
    marginal_process_node_backward(self, node)
  }

  fn process_node_forward(&mut self, graph: &Graph<N, E, ()>, node: &GraphNodeForward<N, E, ()>) -> Result<(), Report> {
    marginal_process_node_forward(self, graph, node)
  }

  fn get_sequence_length(&self) -> usize {
    self.length
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
        .data
        .nodes
        .insert(leaf_key, DenseNodePartition::new(&leaf_fasta.seq, alphabet)?);
    }

    for edge in graph.get_edges() {
      let edge_key = edge.read_arc().key();
      self.data.edges.insert(edge_key, DenseEdgePartition::default());
    }

    Ok(())
  }

  fn extract_ancestral_sequence(&self, node_key: GraphNodeKey) -> Seq {
    if let Some(seq_info) = self.data.nodes.get(&node_key) {
      assign_sequence(seq_info, &self.alphabet)
    } else {
      seq! {}
    }
  }

  fn reconstruct_node_sequence(&mut self, node: &GraphNodeForward<N, E, ()>, include_leaves: bool) -> Option<Seq> {
    if !include_leaves && node.is_leaf {
      return None;
    }

    let seq_info = self.data.nodes.get(&node.key)?;
    let seq = if node.is_leaf {
      seq_info.seq.sequence.clone()
    } else {
      assign_sequence(seq_info, &self.alphabet)
    };

    Some(seq)
  }
}

fn assign_sequence(seq_info: &DenseNodePartition, alphabet: &Alphabet) -> Seq {
  let mut seq = prof2seq(&seq_info.profile, alphabet);
  for gap in &seq_info.seq.gaps {
    seq[gap.0..gap.1].fill(alphabet.gap());
  }
  for unk in &seq_info.seq.unknown {
    seq[unk.0..unk.1].fill(alphabet.unknown());
  }
  seq
}

fn prof2seq(profile: &DenseSeqDistribution, alphabet: &Alphabet) -> Seq {
  let mut seq = seq! {};
  for row in profile.dis.rows() {
    let argmax = argmax_first(&row).unwrap_or(0);
    seq.push(alphabet.char(argmax));
  }
  seq
}

#[cfg(test)]
mod tests {
  use crate::partition::marginal_core::{normalize_from_log, normalize_inplace};
  use crate::pretty_assert_neg_inf;
  use approx::assert_abs_diff_eq;
  use ndarray::{Array2, array};

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
    let mut log_dis = Array2::from_elem((2, 4), f64::NEG_INFINITY);
    log_dis[[1, 0]] = 0.0;
    log_dis[[1, 1]] = -1.0;
    log_dis[[1, 2]] = -2.0;
    log_dis[[1, 3]] = -3.0;

    let (dis, log_lh) = normalize_from_log(&log_dis);

    assert_valid_rows(&dis);
    assert_abs_diff_eq!(dis.row(0), array![0.25, 0.25, 0.25, 0.25], epsilon = 1e-10);
    let sum = 1.0 + (-1.0_f64).exp() + (-2.0_f64).exp() + (-3.0_f64).exp();
    let expected_row1 = array![
      1.0 / sum,
      (-1.0_f64).exp() / sum,
      (-2.0_f64).exp() / sum,
      (-3.0_f64).exp() / sum
    ];
    assert_abs_diff_eq!(dis.row(1), expected_row1, epsilon = 1e-10);
    pretty_assert_neg_inf!(log_lh, "log_lh should be NEG_INFINITY, got {log_lh}");
  }

  #[test]
  fn test_normalize_from_log_mixed_finite_neg_inf_within_row() {
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
    assert_abs_diff_eq!(dis.row(0), array![0.25, 0.25, 0.25, 0.25], epsilon = 1e-10);
    assert_abs_diff_eq!(dis.row(1), array![0.2, 0.4, 0.2, 0.2], epsilon = 1e-10);
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
