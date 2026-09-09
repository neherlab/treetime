use crate::ancestral::sample::SampleMode;
use crate::partition::marginal::shared::pass::{marginal_process_backward_indexed, marginal_process_forward_indexed};
use crate::partition::marginal::sparse::{backward, forward};
use crate::partition::traits::{MarginalPass, PartitionMarginalOps, PartitionMarginalPasses, graph_log_lh};
use eyre::Report;
use log::trace;
use parking_lot::RwLock;
use std::collections::BTreeMap;
use std::sync::Arc;
use treetime_graph::edge::{EdgeOptimizeOps, GraphEdge, GraphEdgeKey, HasBranchLength};
use treetime_graph::graph::Graph;
use treetime_graph::node::{GraphNode, Named};
use treetime_io::fasta::FastaRecord;
use treetime_primitives::{LogLh, Seq, seq};

/// The per-edge branch length each marginal pass propagates sequence profiles along, keyed by edge.
///
/// The value is `profile_branch_length().unwrap_or(0.0)`: the clock-constrained length for a
/// timetree edge (`clock_branch_length` when committed, else the ML or input length) and the ML or
/// input length for every other edge. Collected once at the boundary and threaded into the passes,
/// so the branch length is an explicit pass input rather than a value reached back off the graph
/// edge inside each parallel worker. Callers that maintain their own branch lengths (branch-length
/// optimization, timetree commit) supply their own map instead of this one.
pub fn profile_branch_lengths<N, E, D>(graph: &Graph<N, E, D>) -> BTreeMap<GraphEdgeKey, f64>
where
  N: GraphNode,
  E: GraphEdge + HasBranchLength,
  D: Send + Sync,
{
  graph
    .get_edges()
    .iter()
    .map(|edge| {
      let edge = edge.read_arc();
      let key = edge.key();
      let branch_length = edge.payload().read_arc().profile_branch_length().unwrap_or(0.0);
      (key, branch_length)
    })
    .collect()
}

pub fn initialize_marginal<N, E, P>(
  graph: &Graph<N, E, ()>,
  partitions: &[Arc<RwLock<P>>],
  aln: &[FastaRecord],
) -> Result<LogLh, Report>
where
  N: GraphNode + Named,
  E: EdgeOptimizeOps,
  P: PartitionMarginalOps<N, E> + ?Sized,
{
  for partition in partitions {
    partition.write_arc().attach_sequences(graph, aln)?;
  }
  marginal_update(graph, &profile_branch_lengths(graph), partitions)
}

/// Run the marginal backward and forward passes over the given partitions, propagating profiles
/// along the supplied per-edge branch lengths, and return the substitution log likelihood.
///
/// Branch lengths are an explicit input: the caller decides whether they come from the graph (via
/// [`profile_branch_lengths`]) or from its own store. Each partition contributes its own
/// substitution model; the boundary dispatches dense and sparse representations to their separate
/// tails via [`MarginalPass`].
pub fn marginal_update<N, E, P>(
  graph: &Graph<N, E, ()>,
  branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
  partitions: &[Arc<RwLock<P>>],
) -> Result<LogLh, Report>
where
  N: GraphNode + Named,
  E: EdgeOptimizeOps,
  P: PartitionMarginalPasses<N, E> + ?Sized,
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
  partitions: &[Arc<RwLock<P>>],
) -> Result<(), Report>
where
  N: GraphNode + Named,
  E: EdgeOptimizeOps,
  P: PartitionMarginalPasses<N, E> + ?Sized,
{
  for partition in partitions {
    let mut partition = partition.write_arc();
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
  partitions: &[Arc<RwLock<P>>],
) -> Result<(), Report>
where
  N: GraphNode + Named,
  E: EdgeOptimizeOps,
  P: PartitionMarginalPasses<N, E> + ?Sized,
{
  for partition in partitions {
    let mut partition = partition.write_arc();
    match partition.as_marginal_pass() {
      MarginalPass::Indexed(partition) => marginal_process_forward_indexed(partition, graph, branch_lengths)?,
      MarginalPass::Sparse(partition) => forward::process_forward_indexed(partition, graph, branch_lengths)?,
    }
  }
  Ok(())
}

pub fn ancestral_reconstruction_marginal<N, E, P>(
  graph: &Graph<N, E, ()>,
  include_leaves: bool,
  impute: bool,
  partitions: &[Arc<RwLock<P>>],
  sample_mode: SampleMode,
  rng: &mut dyn rand::RngCore,
  mut visitor: impl FnMut(&N, &Seq) -> Result<(), Report>,
) -> Result<(), Report>
where
  N: GraphNode + Named,
  E: EdgeOptimizeOps,
  P: PartitionMarginalOps<N, E> + ?Sized,
{
  // Preorder traversal is sequential, so a single threaded RNG yields deterministic output under a
  // fixed seed: every node draws from the profile in a fixed traversal order.
  // Reconstruct every node so each partition's stored `seq.sequence` reflects the final, flag-aware
  // reconstruction that the node-data serializer reads back. `include_leaves` gates only whether tip
  // sequences are emitted to the visitor (the reconstructed FASTA), not whether they are computed:
  // `reconstruct_node_sequence` returns `None` for a suppressed tip after writing its `seq.sequence`.
  graph.iter_depth_first_preorder_forward(|node| {
    if partitions.is_empty() {
      if !include_leaves && node.is_leaf {
        return Ok(());
      }
      return visitor(&node.payload, &seq![]);
    }

    let reconstructed = {
      let mut partition = partitions[0].write_arc();
      partition.reconstruct_node_sequence(&node, include_leaves, impute, sample_mode, rng)
    };

    match reconstructed {
      Some(seq) => visitor(&node.payload, &seq),
      None => Ok(()),
    }
  })
}
