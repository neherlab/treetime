#[cfg(test)]
pub mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::commands::ancestral::__tests__::prop_generators::input::MarginalTestInput;
  use crate::seq::alignment::get_common_length;
  use crate::commands::ancestral::marginal::{initialize_marginal, update_marginal};
  use crate::representation::partition::fitch::PartitionFitch;
  use crate::representation::partition::marginal_dense::PartitionMarginalDense;
  use crate::representation::partition::marginal_sparse::PartitionMarginalSparse;
  use crate::representation::payload::ancestral::GraphAncestral;
  use eyre::Report;
  use maplit::btreemap;
  use parking_lot::RwLock;
  use std::sync::Arc;
  use treetime_io::nwk::nwk_read_str;

  /// Run marginal ancestral reconstruction using dense representation.
  ///
  /// Marginal reconstruction computes the posterior distribution P(state|data) at each
  /// internal node independently, integrating over states at all other nodes. The algorithm
  /// uses Felsenstein's pruning (sum-product belief propagation on a tree):
  ///
  ///  - Backward pass (postorder, leaves to root): computes partial likelihood vectors
  ///    (ingroup profiles) at each node from its descendants, using GTR transition
  ///    probability matrices P(t) = exp(Qt) to transform messages along branches.
  ///  - Log-likelihood: computed at the root as sum over sites of log(sum_s pi_s * L_root(s)),
  ///    where pi is the GTR equilibrium frequency vector.
  ///  - Forward pass (preorder, root to leaves): propagates outgroup profiles (information
  ///    from the rest of the tree) to produce full marginal posteriors at every node.
  ///
  /// Dense representation stores full probability vectors at every alignment position for
  /// every node. This is the reference implementation: straightforward but memory-intensive.
  ///
  /// Constructs a `PartitionMarginalDense` from the test input (tree, alignment, GTR model),
  /// runs `initialize_marginal` (attach sequences, then both passes), and returns the total
  /// log-likelihood along with the populated partition for further inspection.
  ///
  /// Used by property tests to verify invariants of marginal ancestral reconstruction.
  pub fn run_dense_marginal(
    input: &MarginalTestInput,
  ) -> Result<(f64, [Arc<RwLock<PartitionMarginalDense>>; 1]), Report> {
    let graph: GraphAncestral = nwk_read_str(&input.newick)?;
    let alphabet = Alphabet::new(AlphabetName::Nuc)?;
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

  /// Run marginal ancestral reconstruction using sparse representation.
  ///
  /// Sparse representation stores probability vectors only at variable (mutated) positions,
  /// with invariant positions handled via Fitch compression. This is the optimized path
  /// for real datasets where most positions are conserved across the tree.
  ///
  /// Two-phase process:
  ///
  ///  1. Fitch parsimony via `compress_sequences`: runs the full Fitch algorithm (attach
  ///     sequences, backward pass, forward pass, cleanup) to reconstruct ancestral states
  ///     and produce the sparse representation - each node stores only mutations relative
  ///     to its parent, not the full sequence.
  ///  2. Marginal reconstruction via `update_marginal`: runs both the backward pass
  ///     (ingroup partial likelihoods) and forward pass (outgroup profiles) on the
  ///     variable positions only, computing the log-likelihood between passes.
  ///
  /// Returns the log-likelihood and the populated partition. Used by property tests
  /// to verify that the sparse path produces results consistent with the dense path.
  pub fn run_sparse_marginal(
    input: &MarginalTestInput,
  ) -> Result<(f64, [Arc<RwLock<PartitionMarginalSparse>>; 1]), Report> {
    let graph: GraphAncestral = nwk_read_str(&input.newick)?;
    let alphabet = Alphabet::default();
    let length = get_common_length(&input.alignment)?;

    let fitch = PartitionFitch::compress(&graph, 0, alphabet, &input.alignment)?;
    let partitions = [Arc::new(RwLock::new(
      fitch.into_marginal_sparse(input.gtr.clone(), &graph)?,
    ))];
    let log_lh = update_marginal(&graph, &partitions)?;
    Ok((log_lh, partitions))
  }
}
