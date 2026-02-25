#[cfg(test)]
mod tests {
  use crate::commands::ancestral::__tests__::prop_generators::input::MarginalTestInput;
  use crate::commands::ancestral::__tests__::prop_marginal_support::{run_dense_marginal, run_sparse_marginal};
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use approx::assert_ulps_eq;
  use eyre::Report;
  use treetime_io::fasta::read_many_fasta_str;

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

  #[test]
  fn test_marginal_dense_sparse_example_gap_free_consistency() -> Result<(), Report> {
    let input = example_gap_free_input()?;
    let (dense_log_lh, _) = run_dense_marginal(&input)?;
    let (sparse_log_lh, _) = run_sparse_marginal(&input)?;

    assert_ulps_eq!(dense_log_lh, sparse_log_lh, epsilon = 1e-10);
    Ok(())
  }
}
