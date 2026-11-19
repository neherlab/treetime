#[cfg(test)]
mod tests {
  use crate::ancestral::__tests__::prop_marginal_support::tests::{
    run_dense_marginal_for_input, run_sparse_marginal_for_input,
  };
  use crate::branch_lengths::branch_lengths_or_zero;
  use crate::pretty_assert_ulps_eq;
  use eyre::Report;
  use treetime_io::nwk::nwk_read;

  #[test]
  fn test_marginal_idempotency_example_dense_repeat_update_keeps_log_lh() -> Result<(), Report> {
    let input = helpers::example_input()?;
    let nwk_parsed = nwk_read(input.newick.as_bytes())?;
    let branch_lengths = branch_lengths_or_zero(&nwk_parsed.branch_lengths);
    let (_, recon) = run_dense_marginal_for_input(&input)?;

    let (recon, log_lh_first) = recon.marginal_update(&nwk_parsed.graph, &branch_lengths)?;
    let (_, log_lh_second) = recon.marginal_update(&nwk_parsed.graph, &branch_lengths)?;
    pretty_assert_ulps_eq!(log_lh_first.value(), log_lh_second.value(), epsilon = 1e-10);

    Ok(())
  }

  #[test]
  fn test_marginal_idempotency_example_sparse_repeat_update_keeps_log_lh() -> Result<(), Report> {
    let input = helpers::example_input()?;
    let nwk_parsed = nwk_read(input.newick.as_bytes())?;
    let branch_lengths = branch_lengths_or_zero(&nwk_parsed.branch_lengths);
    let (_, recon) = run_sparse_marginal_for_input(&input)?;

    let (recon, log_lh_first) = recon.marginal_update(&nwk_parsed.graph, &branch_lengths)?;
    let (_, log_lh_second) = recon.marginal_update(&nwk_parsed.graph, &branch_lengths)?;
    pretty_assert_ulps_eq!(log_lh_first.value(), log_lh_second.value(), epsilon = 1e-10);

    Ok(())
  }

  mod helpers {
    use crate::ancestral::__tests__::prop_generators::input::MarginalTestInput;
    use crate::gtr::get_gtr::{JC69Params, jc69};
    use crate::test_utils::NUC_ALPHABET;
    use eyre::Report;
    use treetime_io::fasta::fasta_read;
    use treetime_primitives::AlignmentRecord;

    pub(super) fn example_input() -> Result<MarginalTestInput, Report> {
      let alignment: Vec<AlignmentRecord> = fasta_read(
        "
>A
ACATCGCCNNA--GAC
>B
GCATCCCTGTA-NG--
>C
CCGGCGATGTRTTG--
>D
TCGGCCGTGTRTTG--
"
        .as_bytes(),
        &*NUC_ALPHABET,
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
  }
}
