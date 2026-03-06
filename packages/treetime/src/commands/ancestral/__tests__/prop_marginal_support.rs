#![cfg(test)]

use crate::alphabet::alphabet::{Alphabet, AlphabetName};
use crate::commands::ancestral::__tests__::prop_generators::input::MarginalTestInput;
use crate::commands::ancestral::fitch::{compress_sequences, get_common_length};
use crate::commands::ancestral::marginal::{initialize_marginal, update_marginal};
use crate::representation::partition::marginal_dense::PartitionMarginalDense;
use crate::representation::partition::marginal_sparse::PartitionMarginalSparse;
use crate::representation::payload::ancestral::GraphAncestral;
use eyre::Report;
use maplit::btreemap;
use parking_lot::RwLock;
use std::sync::Arc;
use treetime_io::nwk::nwk_read_str;

/// Run Felsenstein's pruning algorithm (marginal likelihood) using dense representation.
///
/// Dense representation stores full probability vectors at every alignment position for
/// every node. This is the reference implementation: straightforward but memory-intensive.
///
/// The function constructs a `PartitionMarginalDense` from the test input (tree, alignment,
/// GTR model), runs `initialize_marginal` (backward pass from leaves to root), and returns
/// the total log-likelihood along with the populated partition for further inspection.
///
/// Used by proptest-based tests to verify invariants of marginal ancestral reconstruction.
pub fn run_dense_marginal(
  input: &MarginalTestInput,
) -> Result<(f64, [Arc<RwLock<PartitionMarginalDense>>; 1]), Report> {
  let graph: GraphAncestral = nwk_read_str(&input.newick)?;
  let alphabet = Alphabet::new(AlphabetName::Nuc, true)?;
  let length = get_common_length(&input.alignment)?;

  let partitions = [Arc::new(RwLock::new(PartitionMarginalDense {
    index: 0,
    gtr: input.gtr.clone(),
    alphabet,
    length,
    nodes: btreemap! {},
    edges: btreemap! {},
  }))];

  let log_lh = initialize_marginal(&graph, &partitions, &input.alignment)?;
  Ok((log_lh, partitions))
}

/// Run Felsenstein's pruning algorithm (marginal likelihood) using sparse representation.
///
/// Sparse representation stores probability vectors only at variable (mutated) positions,
/// with invariant positions handled via Fitch compression. This is the optimized path
/// for real datasets where most positions are conserved across the tree.
///
/// The function constructs a `PartitionMarginalSparse` from the test input, runs Fitch
/// compression via `compress_sequences` (to identify variable positions), then runs
/// `update_marginal` (backward pass) to compute the total log-likelihood.
///
/// Returns the log-likelihood and the populated partition. Used by proptest-based tests
/// to verify that the sparse path produces results consistent with the dense path.
pub fn run_sparse_marginal(
  input: &MarginalTestInput,
) -> Result<(f64, [Arc<RwLock<PartitionMarginalSparse>>; 1]), Report> {
  let graph: GraphAncestral = nwk_read_str(&input.newick)?;
  let alphabet = Alphabet::default();
  let length = get_common_length(&input.alignment)?;

  let partitions = [Arc::new(RwLock::new(PartitionMarginalSparse {
    index: 0,
    gtr: input.gtr.clone(),
    alphabet,
    length,
    nodes: btreemap! {},
    edges: btreemap! {},
  }))];

  compress_sequences(&graph, &partitions, &input.alignment)?;
  let log_lh = update_marginal(&graph, &partitions)?;
  Ok((log_lh, partitions))
}
