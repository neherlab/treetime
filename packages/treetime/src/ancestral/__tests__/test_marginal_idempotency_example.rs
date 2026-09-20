#[cfg(test)]
mod tests {
  use crate::ancestral::__tests__::prop_generators::input::MarginalTestInput;
  use crate::ancestral::__tests__::prop_marginal_support::tests::{run_dense_marginal, run_sparse_marginal};
  use crate::ancestral::marginal::branch_lengths_or_zero;
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::pretty_assert_ulps_eq;
  use eyre::Report;
  use treetime_graph::graph::Graph;
  use treetime_io::fasta::read_many_fasta_str;
  use treetime_io::nwk::nwk_read_str;
  use treetime_primitives::AlignmentRecord;

  fn example_input() -> Result<MarginalTestInput, Report> {
    let alignment: Vec<AlignmentRecord> = read_many_fasta_str(
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
    )?
    .into_iter()
    .map(AlignmentRecord::from)
    .collect();
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
  fn test_marginal_idempotency_example_dense() -> Result<(), Report> {
    let input = example_input()?;
    let nwk_parsed = nwk_read_str(&input.newick)?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let graph: Graph = graph;
    let (_, recon) = run_dense_marginal(&input)?;

    let (recon, log_lh_first) = recon.marginal_update(&graph, &branch_lengths_or_zero(&branch_lengths))?;
    let log_lh_first = log_lh_first.value();
    let (recon, log_lh_second) = recon.marginal_update(&graph, &branch_lengths_or_zero(&branch_lengths))?;
    let log_lh_second = log_lh_second.value();
    pretty_assert_ulps_eq!(log_lh_first, log_lh_second, epsilon = 1e-10);

    Ok(())
  }

  #[test]
  fn test_marginal_idempotency_example_sparse() -> Result<(), Report> {
    let input = example_input()?;
    let nwk_parsed = nwk_read_str(&input.newick)?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let graph: Graph = graph;
    let (_, recon) = run_sparse_marginal(&input)?;

    let (recon, log_lh_first) = recon.marginal_update(&graph, &branch_lengths_or_zero(&branch_lengths))?;
    let log_lh_first = log_lh_first.value();
    let (recon, log_lh_second) = recon.marginal_update(&graph, &branch_lengths_or_zero(&branch_lengths))?;
    let log_lh_second = log_lh_second.value();
    pretty_assert_ulps_eq!(log_lh_first, log_lh_second, epsilon = 1e-10);

    Ok(())
  }
}
