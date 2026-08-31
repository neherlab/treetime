use crate::ancestral::fitch_indel::{compute_node_ranges, resolve_indels_backward, resolve_indels_forward};
use crate::partition::marginal::dense::partition::{PartitionMarginalDense, assign_sequence};
use crate::partition::marginal::shared::data::{IndexedMarginalPartition, MarginalData, MarginalPartition};
use crate::partition::marginal::shared::pass::{marginal_process_backward_indexed, marginal_process_forward_indexed};
use crate::partition::storage::dense::{DenseEdgePartition, DenseNodePartition, DenseSeqDistribution, DenseSeqInfo};
use crate::partition::traits::PartitionMarginalPasses;
use eyre::Report;
use itertools::Itertools;
use std::collections::BTreeMap;
use treetime_graph::edge::{EdgeOptimizeOps, GraphEdgeKey};
use treetime_graph::graph::Graph;
use treetime_graph::node::{GraphNodeKey, NodeAncestralOps};
use treetime_primitives::LogLh;
use treetime_utils::interval::range_union::range_union;

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

  fn indexed_storage_mut(
    &mut self,
  ) -> (
    &mut BTreeMap<GraphNodeKey, DenseNodePartition>,
    &mut BTreeMap<GraphEdgeKey, DenseEdgePartition>,
  ) {
    (&mut self.data.nodes, &mut self.data.edges)
  }
}

impl<N, E> IndexedMarginalPartition<N, E> for PartitionMarginalDense
where
  N: NodeAncestralOps,
  E: EdgeOptimizeOps,
{
  fn indexed_missing_node(&self, _key: GraphNodeKey) -> Result<DenseNodePartition, Report> {
    Ok(DenseNodePartition {
      seq: DenseSeqInfo::default(),
      profile: DenseSeqDistribution::default(),
    })
  }

  fn indexed_leaf_profile(&self, node: &DenseNodePartition) -> Result<DenseSeqDistribution, Report> {
    Ok(DenseSeqDistribution {
      dis: self.alphabet.seq2prof(&node.seq.sequence)?,
      log_lh: LogLh::ZERO,
    })
  }

  fn indexed_backward_internal(&self, children: &[&DenseNodePartition]) -> Result<DenseSeqInfo, Report> {
    let child_non_chars = children.iter().map(|child| &child.seq.non_char).collect_vec();
    let child_gaps = children.iter().map(|child| &child.seq.gaps).collect_vec();
    let ranges = compute_node_ranges(&child_non_chars, &child_gaps);
    let child_unknown = children.iter().map(|child| &child.seq.unknown).collect_vec();
    let child_variable_indels = children.iter().map(|child| &child.seq.variable_indel).collect_vec();
    let indels = resolve_indels_backward(&child_gaps, &child_unknown, &child_variable_indels, self.length);
    Ok(DenseSeqInfo {
      gaps: indels.resolved_gaps,
      unknown: ranges.unknown,
      non_char: ranges.non_char,
      variable_indel: indels.variable_indel,
      ..DenseSeqInfo::default()
    })
  }

  fn indexed_forward_post(
    &self,
    is_root: bool,
    is_leaf: bool,
    parent: Option<&DenseNodePartition>,
    node: &mut DenseNodePartition,
    parent_edge: Option<&mut DenseEdgePartition>,
  ) -> Result<(), Report> {
    if is_root {
      node.seq.variable_indel.clear();
      node.seq.sequence = assign_sequence(node, &self.alphabet);
      node.seq.non_char = range_union(&[node.seq.gaps.clone(), node.seq.unknown.clone()]);
      return Ok(());
    }

    if !is_leaf {
      node.seq.sequence = assign_sequence(node, &self.alphabet);
    }
    let parent = parent.expect("Non-root dense node must have a parent");
    let variable_indel = std::mem::take(&mut node.seq.variable_indel);
    let (indels, new_gaps) = resolve_indels_forward(
      &variable_indel,
      &node.seq.gaps,
      &node.seq.non_char,
      &parent.seq.gaps,
      &parent.seq.sequence,
      &node.seq.sequence,
    );
    node.seq.gaps = new_gaps;
    node.seq.non_char = range_union(&[node.seq.gaps.clone(), node.seq.unknown.clone()]);
    for gap in &node.seq.gaps {
      node.seq.sequence[gap.0..gap.1].fill(self.alphabet.gap());
    }
    for unknown in &node.seq.unknown {
      node.seq.sequence[unknown.0..unknown.1].fill(self.alphabet.unknown());
    }
    parent_edge
      .expect("Non-root dense node must own its parent edge")
      .indels = indels;
    Ok(())
  }
}

impl<N, E> PartitionMarginalPasses<N, E> for PartitionMarginalDense
where
  N: NodeAncestralOps,
  E: EdgeOptimizeOps,
{
  fn process_backward_pass(&mut self, graph: &Graph<N, E, ()>) -> Result<(), Report> {
    marginal_process_backward_indexed(self, graph)
  }

  fn process_forward_pass(&mut self, graph: &Graph<N, E, ()>) -> Result<(), Report> {
    marginal_process_forward_indexed(self, graph)
  }

  fn get_sequence_length(&self) -> usize {
    self.length
  }
}
