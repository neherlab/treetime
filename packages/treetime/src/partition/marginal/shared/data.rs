use crate::constants::SUPERTINY_NUMBER;
use crate::gtr::gtr::GTR;
use crate::gtr::infer_gtr::common::{
  MutationCounts, accumulate_mutation_counts, get_branch_mutation_matrix, is_profile_informative,
};
use crate::partition::storage::dense::{DenseEdgePartition, DenseNodePartition, DenseSeqDistribution, DenseSeqInfo};
use eyre::Report;
use ndarray::prelude::*;
use serde::Serialize;
use std::collections::BTreeMap;
use treetime_graph::edge::{EdgeOptimizeOps, GraphEdge, GraphEdgeKey, HasBranchLength};
use treetime_graph::graph::Graph;
use treetime_graph::node::{GraphNode, GraphNodeKey, Named};
use treetime_utils::array::ndarray::argmax_first;

#[derive(Clone, Debug, Serialize)]
pub struct MarginalData {
  pub gtr: GTR,
  pub nodes: BTreeMap<GraphNodeKey, DenseNodePartition>,
  pub edges: BTreeMap<GraphEdgeKey, DenseEdgePartition>,
  pub min_branch_length: f64,
  /// When `true`, root positions whose posterior profile is essentially uniform
  /// are excluded from the equilibrium-frequency prior in `count_transitions`.
  /// v0 never filters (always folds in the root MAP state); the nucleotide
  /// ancestral path enables filtering to drop signal-free gap-only columns.
  pub filter_uninformative_root: bool,
}

impl MarginalData {
  pub fn effective_branch_length(&self, raw: f64) -> f64 {
    raw.max(self.min_branch_length)
  }

  /// Count posterior-weighted transitions from dense profile matrices.
  ///
  /// Shared by dense and discrete partitions (both store full profile matrices).
  pub fn count_transitions<N, E>(&self, graph: &Graph<N, E, ()>) -> Result<MutationCounts, Report>
  where
    N: GraphNode,
    E: GraphEdge + HasBranchLength,
  {
    let n_states = self.gtr.pi.len();
    let mut nij = Array2::zeros((n_states, n_states));
    let mut Ti = Array1::zeros(n_states);

    for edge in graph.get_edges() {
      let edge_arc = edge.read_arc();
      let branch_length = self.effective_branch_length(edge_arc.payload().read_arc().branch_length().unwrap_or(0.0));
      let edge_key = edge_arc.key();

      let edge_data = &self.edges[&edge_key];

      let exp_qt = self.gtr.expQt(branch_length) + SUPERTINY_NUMBER;
      let mut_stack = get_branch_mutation_matrix(&edge_data.msg_to_child.dis, &edge_data.msg_to_parent.dis, &exp_qt);
      accumulate_mutation_counts(&mut_stack, branch_length, &mut nij, &mut Ti);
    }

    let root = graph.get_exactly_one_root()?;
    let root_key = root.read_arc().key();
    let root_profile = &self.nodes[&root_key].profile.dis;
    let mut root_state = Array1::zeros(n_states);
    for row in root_profile.rows() {
      if !self.filter_uninformative_root || is_profile_informative(row, n_states) {
        if let Some(root_idx) = argmax_first(&row) {
          root_state[root_idx] += 1.0;
        }
      }
    }

    nij.diag_mut().fill(0.0);

    Ok(MutationCounts { nij, Ti, root_state })
  }
}

pub trait MarginalPartition<N, E>: Send + Sync
where
  N: GraphNode + Named,
  E: EdgeOptimizeOps,
{
  fn marginal_data(&self) -> &MarginalData;
  fn marginal_data_mut(&mut self) -> &mut MarginalData;

  fn indexed_storage_mut(
    &mut self,
  ) -> (
    &mut BTreeMap<GraphNodeKey, DenseNodePartition>,
    &mut BTreeMap<GraphEdgeKey, DenseEdgePartition>,
  );
}

pub trait IndexedMarginalPartition<N, E>: MarginalPartition<N, E>
where
  N: GraphNode + Named,
  E: EdgeOptimizeOps,
{
  fn indexed_missing_node(&self, key: GraphNodeKey) -> Result<DenseNodePartition, Report>;

  fn indexed_leaf_profile(&self, node: &DenseNodePartition) -> Result<DenseSeqDistribution, Report>;

  fn indexed_backward_internal(&self, _children: &[&DenseNodePartition]) -> Result<DenseSeqInfo, Report> {
    Ok(DenseSeqInfo::default())
  }

  fn indexed_forward_post(
    &self,
    is_root: bool,
    is_leaf: bool,
    parent: Option<&DenseNodePartition>,
    node: &mut DenseNodePartition,
    parent_edge: Option<&mut DenseEdgePartition>,
  ) -> Result<(), Report>;
}
