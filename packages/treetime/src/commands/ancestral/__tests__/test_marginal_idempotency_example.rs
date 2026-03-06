#[cfg(test)]
mod tests {
  use crate::commands::ancestral::__tests__::prop_generators::input::MarginalTestInput;
  use crate::commands::ancestral::__tests__::prop_marginal_support::{run_dense_marginal, run_sparse_marginal};
  use crate::commands::ancestral::marginal::update_marginal;
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::representation::payload::ancestral::GraphAncestral;
  use approx::assert_ulps_eq;
  use eyre::Report;
  use treetime_io::fasta::read_many_fasta_str;
  use treetime_io::nwk::nwk_read_str;

  /// Build a fixed 4-taxon test input for marginal idempotency verification.
  ///
  /// Tree topology: `((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01`
  /// Alignment: 16-character nucleotide sequences containing standard bases,
  /// ambiguity codes (N, R), and gap characters (-).
  /// Model: JC69 (Jukes-Cantor 1969) with equal equilibrium frequencies
  /// pi = [0.25, 0.25, 0.25, 0.25] and uniform substitution rates.
  fn example_input() -> Result<MarginalTestInput, Report> {
    let alignment = read_many_fasta_str(
      "
>A
ACATCGCCNNA--GAC
>B
GCATCCCTGTA-NG--
>C
CCGGCGATGTRTTG--
>D
TCGGCCGTGTRTTG--
",
      &*crate::test_utils::NUC_ALPHABET,
    )?;
    let gtr = jc69(JC69Params::default())?;
    Ok(MarginalTestInput {
      newick: "((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;".to_owned(),
      alignment,
      gtr,
      n_taxa: 4,
      seq_len: 16,
    })
  }

  /// Verify that `update_marginal` is a fixed-point operation on a dense partition.
  ///
  /// Felsenstein's pruning algorithm computes exact marginal likelihoods in a single
  /// upward pass (leaves to root) followed by a downward pass (root to leaves). If the
  /// implementation is correct, re-running `update_marginal` on already-computed partition
  /// data must produce an identical log-likelihood.
  ///
  /// Invariant: `L(update_marginal(P)) == L(update_marginal(update_marginal(P)))`,
  /// where L is the total log-likelihood and P is the partition state.
  ///
  /// Any deviation indicates state corruption, accumulating numerical drift, or
  /// incorrect in-place mutation of node/edge profiles.
  ///
  /// Uses a fixed 4-taxon tree with JC69 model and dense (full profile matrix)
  /// representation.
  #[test]
  fn test_marginal_idempotency_example_dense() -> Result<(), Report> {
    let input = example_input()?;
    let graph: GraphAncestral = nwk_read_str(&input.newick)?;
    let (_, partitions) = run_dense_marginal(&input)?;

    let log_lh_first = update_marginal(&graph, &partitions)?;
    let log_lh_second = update_marginal(&graph, &partitions)?;
    assert_ulps_eq!(log_lh_first, log_lh_second, epsilon = 1e-10);

    Ok(())
  }

  /// Verify that `update_marginal` is a fixed-point operation on a sparse partition.
  ///
  /// Same idempotency invariant as the dense variant: two consecutive calls to
  /// `update_marginal` must yield identical log-likelihoods. The sparse representation
  /// stores only variable positions and groups fixed characters, introducing different
  /// code paths for accumulation and normalization that must also be idempotent.
  ///
  /// Invariant: `L(update_marginal(P)) == L(update_marginal(update_marginal(P)))`.
  ///
  /// Uses the same fixed 4-taxon tree and JC69 model as the dense variant.
  #[test]
  fn test_marginal_idempotency_example_sparse() -> Result<(), Report> {
    let input = example_input()?;
    let graph: GraphAncestral = nwk_read_str(&input.newick)?;
    let (_, partitions) = run_sparse_marginal(&input)?;

    let log_lh_first = update_marginal(&graph, &partitions)?;
    let log_lh_second = update_marginal(&graph, &partitions)?;
    assert_ulps_eq!(log_lh_first, log_lh_second, epsilon = 1e-10);

    Ok(())
  }
}
