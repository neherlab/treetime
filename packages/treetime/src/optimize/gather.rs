//! Pre-gathered per-edge value maps the branch-length optimizer consumes.
//!
//! The optimizer reads, per edge, each partition's substitution-likelihood contribution and indel
//! count, and one run-level total sequence length. Gathering these into value maps once, off the
//! reconstructions, lets the optimizer and the indel-likelihood evaluators operate on plain maps
//! instead of reaching back into a partition edge by edge.
//!
//! Contribution order is significant: the mixed-partition evaluator folds per-site likelihoods across
//! partitions, so the per-edge contribution vector is built dense entries first then sparse, matching
//! the order the optimizer summed them before. Indel counts, substitution counts, and effective
//! lengths are integer sums and so order-independent.

use crate::ancestral::pipeline::{DenseReconstruction, SparseReconstruction};
use crate::partition::optimize::contribution::OptimizationContribution;
use crate::partition::timetree::partition::PartitionTimetree;
use eyre::Report;
use rayon::prelude::*;
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;

/// Per-edge branch-length optimization contributions from every partition, dense entries first then
/// sparse.
pub fn gather_edge_contributions(
  graph: &Graph,
  dense: &[DenseReconstruction],
  sparse: &[SparseReconstruction],
) -> Result<BTreeMap<GraphEdgeKey, Vec<OptimizationContribution>>, Report> {
  graph
    .get_edges()
    .collect::<Vec<_>>()
    .into_par_iter()
    .map(
      |edge_ref| -> Result<(GraphEdgeKey, Vec<OptimizationContribution>), Report> {
        let edge_key = edge_ref.key();
        let mut contributions = Vec::with_capacity(dense.len() + sparse.len());
        for family in dense {
          contributions.push(family.create_edge_contribution(edge_key));
        }
        for family in sparse {
          contributions.push(family.create_edge_contribution(edge_key)?);
        }
        Ok((edge_key, contributions))
      },
    )
    .collect()
}

/// Summed indel event count per edge across every partition.
pub fn gather_edge_indel_counts(
  graph: &Graph,
  dense: &[DenseReconstruction],
  sparse: &[SparseReconstruction],
) -> BTreeMap<GraphEdgeKey, usize> {
  graph
    .get_edges()
    .collect::<Vec<_>>()
    .into_par_iter()
    .map(|edge_ref| {
      let edge_key = edge_ref.key();
      let count = dense
        .iter()
        .map(|family| family.edge_indel_count(edge_key))
        .sum::<usize>()
        + sparse
          .iter()
          .map(|family| family.edge_indel_count(edge_key))
          .sum::<usize>();
      (edge_key, count)
    })
    .collect()
}

/// Summed substitution count per edge (`edge_subs().len()`) across every partition.
pub fn gather_edge_sub_counts(
  graph: &Graph,
  dense: &[DenseReconstruction],
  sparse: &[SparseReconstruction],
) -> Result<BTreeMap<GraphEdgeKey, usize>, Report> {
  graph
    .get_edges()
    .collect::<Vec<_>>()
    .into_par_iter()
    .map(|edge_ref| -> Result<(GraphEdgeKey, usize), Report> {
      let edge_key = edge_ref.key();
      let mut count = 0;
      for family in dense {
        count += family.edge_subs(graph, edge_key)?.len();
      }
      for family in sparse {
        count += family.edge_subs(edge_key)?.len();
      }
      Ok((edge_key, count))
    })
    .collect()
}

/// Summed effective (canonical) alignment length per edge across every partition.
pub fn gather_edge_effective_lengths(
  graph: &Graph,
  dense: &[DenseReconstruction],
  sparse: &[SparseReconstruction],
) -> Result<BTreeMap<GraphEdgeKey, usize>, Report> {
  graph
    .get_edges()
    .collect::<Vec<_>>()
    .into_par_iter()
    .map(|edge_ref| -> Result<(GraphEdgeKey, usize), Report> {
      let edge_key = edge_ref.key();
      let mut total = 0;
      for family in dense {
        total += family.edge_effective_length(graph, edge_key)?;
      }
      for family in sparse {
        total += family.edge_effective_length(graph, edge_key)?;
      }
      Ok((edge_key, total))
    })
    .collect()
}

/// Total sequence length across every partition.
pub fn total_sequence_length(dense: &[DenseReconstruction], sparse: &[SparseReconstruction]) -> usize {
  dense.iter().map(DenseReconstruction::sequence_length).sum::<usize>()
    + sparse.iter().map(SparseReconstruction::sequence_length).sum::<usize>()
}

/// Per-edge optimization contributions from a timetree partition set, in partition order.
pub fn gather_timetree_edge_contributions(
  graph: &Graph,
  partitions: &[PartitionTimetree],
) -> Result<BTreeMap<GraphEdgeKey, Vec<OptimizationContribution>>, Report> {
  graph
    .get_edges()
    .collect::<Vec<_>>()
    .into_par_iter()
    .map(
      |edge_ref| -> Result<(GraphEdgeKey, Vec<OptimizationContribution>), Report> {
        let edge_key = edge_ref.key();
        let contributions = partitions
          .iter()
          .map(|partition| partition.create_edge_contribution(edge_key))
          .collect::<Result<Vec<_>, _>>()?;
        Ok((edge_key, contributions))
      },
    )
    .collect()
}

/// Summed indel event count per edge across a timetree partition set.
pub fn gather_timetree_edge_indel_counts(
  graph: &Graph,
  partitions: &[PartitionTimetree],
) -> BTreeMap<GraphEdgeKey, usize> {
  graph
    .get_edges()
    .collect::<Vec<_>>()
    .into_par_iter()
    .map(|edge_ref| {
      let edge_key = edge_ref.key();
      (
        edge_key,
        partitions
          .iter()
          .map(|partition| partition.edge_indel_count(edge_key))
          .sum(),
      )
    })
    .collect()
}

/// Total sequence length across a timetree partition set.
pub fn timetree_total_sequence_length(partitions: &[PartitionTimetree]) -> usize {
  partitions.iter().map(PartitionTimetree::sequence_length).sum()
}
