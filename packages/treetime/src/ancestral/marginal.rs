use crate::ancestral::sample::SampleMode;
use crate::partition::marginal::shared::pass::{marginal_process_backward_indexed, marginal_process_forward_indexed};
use crate::partition::marginal::sparse::{backward, forward};
use crate::partition::traits::{MarginalPass, PartitionMarginalOps, PartitionMarginalPasses, graph_log_lh};
use eyre::Report;
use log::trace;
use std::collections::BTreeMap;
use treetime_graph::edge::{GraphEdge, GraphEdgeKey};
use treetime_graph::graph::Graph;
use treetime_graph::node::{GraphNode, GraphNodeKey};
use treetime_io::fasta::FastaRecord;
use treetime_primitives::{LogLh, Seq, seq};

/// The branch length that propagates sequence profiles along an edge: the clock-constrained length
/// when one has been committed, otherwise the raw ML or input length.
///
/// This is the domain choice that a timetree edge makes (`clock_branch_length` over the raw length);
/// every other command has no clock length and falls back to the raw length. Kept as a free function
/// so the choice stays named and testable where a profile map is derived.
pub fn profile_branch_length(clock: Option<f64>, raw: Option<f64>) -> Option<f64> {
  clock.or(raw)
}

/// Derive the per-edge profile branch length map (`f64`) each marginal pass propagates sequence
/// profiles along, from the raw input-tree branch length value map.
///
/// The value is `raw.unwrap_or(0.0)`: an edge with no length resolves to `0.0` for the passes. For a
/// timetree the clock-constrained length is combined in separately via
/// [`timetree_branch_lengths`](crate::timetree::inference::runner::timetree_branch_lengths); this
/// helper serves the non-timetree marginal passes (ancestral, optimize, clock, mugration).
pub fn profile_branch_lengths(branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>) -> BTreeMap<GraphEdgeKey, f64> {
  branch_lengths
    .iter()
    .map(|(key, raw)| (*key, profile_branch_length(None, *raw).unwrap_or(0.0)))
    .collect()
}

pub fn initialize_marginal<N, E, P>(
  graph: &Graph<N, E, ()>,
  branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
  partitions: &mut [P],
  aln: &[FastaRecord],
  names: &BTreeMap<GraphNodeKey, Option<String>>,
) -> Result<LogLh, Report>
where
  N: GraphNode,
  E: GraphEdge,
  P: PartitionMarginalOps<N, E>,
{
  for partition in partitions.iter_mut() {
    partition.attach_sequences(graph, aln, names)?;
  }
  marginal_update(graph, branch_lengths, partitions)
}

/// Run the marginal backward and forward passes over the given partitions, propagating profiles
/// along the supplied per-edge branch lengths, and return the substitution log likelihood.
///
/// Branch lengths are an explicit input: the caller decides whether they come from the parsed
/// input tree (via [`profile_branch_lengths`]) or from its own store. Each partition contributes its own
/// substitution model; the boundary dispatches dense and sparse representations to their separate
/// tails via [`MarginalPass`].
pub fn marginal_update<N, E, P>(
  graph: &Graph<N, E, ()>,
  branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
  partitions: &mut [P],
) -> Result<LogLh, Report>
where
  N: GraphNode,
  E: GraphEdge,
  P: PartitionMarginalPasses<N, E>,
{
  marginal_backward(graph, branch_lengths, partitions)?;
  let log_lh = graph_log_lh(graph, partitions)?;
  trace!("Marginal log likelihood (substitution): {}", log_lh.value());
  marginal_forward(graph, branch_lengths, partitions)?;
  Ok(log_lh)
}

pub fn marginal_backward<N, E, P>(
  graph: &Graph<N, E, ()>,
  branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
  partitions: &mut [P],
) -> Result<(), Report>
where
  N: GraphNode,
  E: GraphEdge,
  P: PartitionMarginalPasses<N, E>,
{
  for partition in partitions.iter_mut() {
    match partition.as_marginal_pass() {
      MarginalPass::Indexed(partition) => marginal_process_backward_indexed(partition, graph, branch_lengths)?,
      MarginalPass::Sparse(partition) => backward::process_backward_indexed(partition, graph, branch_lengths)?,
    }
  }
  Ok(())
}

fn marginal_forward<N, E, P>(
  graph: &Graph<N, E, ()>,
  branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
  partitions: &mut [P],
) -> Result<(), Report>
where
  N: GraphNode,
  E: GraphEdge,
  P: PartitionMarginalPasses<N, E>,
{
  for partition in partitions.iter_mut() {
    match partition.as_marginal_pass() {
      MarginalPass::Indexed(partition) => marginal_process_forward_indexed(partition, graph, branch_lengths)?,
      MarginalPass::Sparse(partition) => forward::process_forward_indexed(partition, graph, branch_lengths)?,
    }
  }
  Ok(())
}

/// Reconstruct ancestral sequences with marginal inference, emitting each to `visitor` and
/// returning the reconstructed sequences keyed by node id.
///
/// The returned map is the reconstruction result as a value: the command captures it instead of
/// having to read it back out of partition state. Every node whose sequence is emitted to the
/// visitor is also recorded in the map, so the two views hold the same sequences. Until the tree
/// writers read the map directly, the partition still stores each `seq.sequence` (written inside
/// `reconstruct_node_sequence`) for the node-data serializer.
pub fn ancestral_reconstruction_marginal<N, E, P>(
  graph: &Graph<N, E, ()>,
  include_leaves: bool,
  impute: bool,
  partitions: &mut [P],
  sample_mode: SampleMode,
  rng: &mut dyn rand::RngCore,
  mut visitor: impl FnMut(GraphNodeKey, &Seq) -> Result<(), Report>,
) -> Result<BTreeMap<GraphNodeKey, Seq>, Report>
where
  N: GraphNode,
  E: GraphEdge,
  P: PartitionMarginalOps<N, E>,
{
  // Preorder traversal is sequential, so a single threaded RNG yields deterministic output under a
  // fixed seed: every node draws from the profile in a fixed traversal order.
  // Reconstruct every node so each partition's stored `seq.sequence` reflects the final, flag-aware
  // reconstruction that the node-data serializer reads back. `include_leaves` gates only whether tip
  // sequences are emitted to the visitor (the reconstructed FASTA), not whether they are computed:
  // `reconstruct_node_sequence` returns `None` for a suppressed tip after writing its `seq.sequence`.
  let mut node_sequences = BTreeMap::new();
  graph.iter_depth_first_preorder_forward(|node| {
    if partitions.is_empty() {
      if !include_leaves && node.is_leaf {
        return Ok(());
      }
      let seq = seq![];
      visitor(node.key, &seq)?;
      node_sequences.insert(node.key, seq);
      return Ok(());
    }

    let reconstructed = partitions[0].reconstruct_node_sequence(&node, include_leaves, impute, sample_mode, rng);

    match reconstructed {
      Some(seq) => {
        visitor(node.key, &seq)?;
        node_sequences.insert(node.key, seq);
        Ok(())
      },
      None => Ok(()),
    }
  })?;
  Ok(node_sequences)
}
