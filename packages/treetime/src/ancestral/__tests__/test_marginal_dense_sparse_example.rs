#[cfg(test)]
mod tests {
  use crate::ancestral::__tests__::prop_generators::input::MarginalTestInput;
  use crate::ancestral::__tests__::prop_marginal_support::tests::{run_dense_marginal, run_sparse_marginal};
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::pretty_assert_ulps_eq;
  use eyre::Report;
  use treetime_io::fasta::read_many_fasta_str;

  /// Construct a deterministic test input: 4-taxon balanced binary tree with
  /// JC69 model and a 16bp gap-free alignment.
  ///
  /// All four sequences share a 15bp invariant prefix and differ only at
  /// position 16 (A=T, B=A, C=G, D=C), creating a single variable site with
  /// all four nucleotide states represented. This exercises both the invariant-
  /// site accumulation path and the per-site Felsenstein pruning at the
  /// variable position.
  ///
  /// Branch lengths are asymmetric (0.05-0.2 subs/site), creating non-uniform
  /// transition probability matrices P(t) = exp(Q*t) across edges.
  ///
  /// JC69 (Jukes-Cantor 1969: uniform equilibrium pi=1/4, equal
  /// exchangeabilities) is the simplest GTR submodel, isolating dense/sparse
  /// agreement from GTR parameterization effects.
  fn example_gap_free_input() -> Result<MarginalTestInput, Report> {
    let alignment = read_many_fasta_str(
      "
>A
ACGTACGTACGTACGT
>B
ACGTACGTACGTACGA
>C
ACGTACGTACGTACGG
>D
ACGTACGTACGTACGC
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

  /// Deterministic dense/sparse log-likelihood consistency test on a fixed input.
  ///
  /// Companion to the property-based test
  /// `test_prop_marginal_dense_sparse_gap_free_consistency` which uses random
  /// GTR models and tree topologies. This test uses a fixed 4-taxon tree, 16bp
  /// alignment, and JC69 model to provide a stable, debuggable baseline. Fixed
  /// inputs make failures reproducible without proptest shrinking.
  ///
  /// Both implementations compute the Felsenstein pruning algorithm (Felsenstein
  /// 1981) via different code paths. Dense stores a K-state probability vector
  /// (K=4 for nucleotides) at every alignment position for every tree node.
  /// Sparse compresses invariant sites via Fitch parsimony and stores
  /// probability vectors only at variable positions, accumulating invariant-site
  /// contributions separately. The total log-likelihood
  /// ln P(D|T) = sum_i ln(sum_x pi[x] * L_root_i(x)), where i indexes sites
  /// and x indexes states, must agree between implementations.
  ///
  /// Gap-free alignment is required because sparse compression of gapped
  /// positions follows a different code path.
  ///
  /// Tolerance: `epsilon = 1e-10` absolute difference OR <= 4 ULPs (the
  /// `pretty_assert_ulps_eq!` default `max_ulps`). With JC69, both paths produce
  /// near-identical floating-point results (operation ordering differences only).
  #[test]
  fn test_marginal_dense_sparse_example_gap_free_consistency() -> Result<(), Report> {
    let input = example_gap_free_input()?;
    let (dense_log_lh, _) = run_dense_marginal(&input)?;
    let (sparse_log_lh, _) = run_sparse_marginal(&input)?;

    pretty_assert_ulps_eq!(dense_log_lh, sparse_log_lh, epsilon = 1e-10);
    Ok(())
  }
}
