#[cfg(test)]
mod tests {
  use crate::commands::ancestral::__tests__::prop_generators::input::MarginalTestInput;
  use crate::commands::ancestral::__tests__::prop_marginal_support::{run_dense_marginal, run_sparse_marginal};
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use approx::assert_ulps_eq;
  use eyre::Report;
  use treetime_io::fasta::read_many_fasta_str;

  /// Construct a deterministic test input: 4-taxon balanced tree with JC69 model
  /// and a 16bp gap-free alignment.
  ///
  /// All four sequences share a 15bp prefix and differ only at position 16
  /// (A=T, B=A, C=G, D=C), creating a single variable site. The JC69 model
  /// (uniform equilibrium, equal rates) provides the simplest non-trivial
  /// substitution model for verifying dense/sparse agreement.
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

  /// Deterministic dense/sparse log-likelihood consistency test on a known input.
  ///
  /// Companion to the property-based test
  /// `test_prop_marginal_dense_sparse_gap_free_consistency` which uses random
  /// inputs. This test uses a fixed 4-taxon tree, 16bp alignment, and JC69 model
  /// to provide a stable, debuggable baseline for dense/sparse agreement.
  ///
  /// Both implementations compute the Felsenstein likelihood via different code
  /// paths (full NxK matrices vs compressed variable positions). The total
  /// log-likelihood L = sum over sites of log(sum_s pi[s] * L_root(s)) must
  /// agree to floating-point precision.
  #[test]
  fn test_marginal_dense_sparse_example_gap_free_consistency() -> Result<(), Report> {
    let input = example_gap_free_input()?;
    let (dense_log_lh, _) = run_dense_marginal(&input)?;
    let (sparse_log_lh, _) = run_sparse_marginal(&input)?;

    assert_ulps_eq!(dense_log_lh, sparse_log_lh, epsilon = 1e-10);
    Ok(())
  }
}
