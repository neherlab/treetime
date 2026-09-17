#[cfg(test)]
pub mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::ancestral::__tests__::prop_generators::input::MarginalTestInput;
  use crate::ancestral::fitch::create_fitch_partition;
  use crate::ancestral::marginal::branch_lengths_or_zero;
  use crate::ancestral::pipeline::{DenseReconstruction, SparseReconstruction};
  use crate::partition::marginal::dense::partition::PartitionMarginalDense;
  use crate::seq::alignment::get_common_length;
  use crate::seq::alignment::node_seq_inputs;
  use eyre::Report;

  use treetime_graph::graph::Graph;

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
  pub fn run_dense_marginal(input: &MarginalTestInput) -> Result<(f64, DenseReconstruction), Report> {
    let nwk_parsed = nwk_read_str(&input.newick)?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let graph: Graph = graph;
    let alphabet = Alphabet::new(AlphabetName::Nuc)?;
    let length = get_common_length(&input.alignment)?;

    let partition = PartitionMarginalDense::new(0, alphabet, length);
    let node_states = partition.attach_sequences(&graph, &node_seq_inputs(&graph, &names, input.alignment.clone()))?;
    let recon = DenseReconstruction::seeded(partition, input.gtr.clone(), node_states);
    let (recon, log_lh) = recon.marginal_update(&graph, &branch_lengths_or_zero(&branch_lengths))?;
    let log_lh = log_lh.value();
    Ok((log_lh, recon))
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
  ///  2. Marginal reconstruction via `marginal_update`: runs both the backward pass
  ///     (ingroup partial likelihoods) and forward pass (outgroup profiles) on the
  ///     variable positions only, computing the log-likelihood between passes.
  ///
  /// Returns the log-likelihood and the populated partition. Used by property tests
  /// to verify that the sparse path produces results consistent with the dense path.
  pub fn run_sparse_marginal(input: &MarginalTestInput) -> Result<(f64, SparseReconstruction), Report> {
    let nwk_parsed = nwk_read_str(&input.newick)?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let graph: Graph = graph;
    let alphabet = Alphabet::default();
    let _ = get_common_length(&input.alignment)?;

    let fitch = create_fitch_partition(
      &graph,
      0,
      alphabet,
      &node_seq_inputs(&graph, &names, input.alignment.clone()),
    )?;
    let (partition, node_states) = fitch.into_marginal_sparse(&graph)?;
    let recon = SparseReconstruction::seeded(partition, input.gtr.clone(), node_states);
    let (recon, log_lh) = recon.marginal_update(&graph, &branch_lengths_or_zero(&branch_lengths))?;
    let log_lh = log_lh.value();
    Ok((log_lh, recon))
  }
}
