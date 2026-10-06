#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::branch_lengths::branch_lengths_or_zero;
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::partition::marginal::dense::partition::PartitionMarginalDense;
  use crate::partition::marginal::reconstruction::{DenseReconstruction, MarginalReconstruction};
  use crate::pretty_assert_abs_diff_eq;
  use crate::reroot::orchestrate::{RerootTopologyParams, reroot_at_node};
  use crate::test_utils::leaf_seq_inputs;
  use crate::test_utils::{NUC_ALPHABET, find_node_key_by_name};
  use eyre::Report;
  use indoc::indoc;
  use treetime_io::fasta::fasta_read;
  use treetime_io::nwk::nwk_read;
  use treetime_primitives::AlignmentRecord;

  #[test]
  fn test_reroot_dense_reconstruction_matches_a_fresh_reconstruction_on_the_rerooted_tree() -> Result<(), Report> {
    let nwk_parsed = nwk_read(b"((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root;".as_slice())?;
    let names = nwk_parsed.names();
    let mut graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
    let aln: Vec<AlignmentRecord> = fasta_read(
      indoc! {"
        >A
        ACGTACGTACGTACGT
        >B
        ACGTACCTACGTACGA
        >C
        ACTTACGTACGAACGG
        >D
        ACTTACGTACGAACGC
      "}
      .as_bytes(),
      &*NUC_ALPHABET,
    )?
    .into_iter()
    .map(AlignmentRecord::from)
    .collect();
    let gtr = jc69(JC69Params::default())?;
    let alphabet = Alphabet::new(AlphabetName::Nuc)?;
    let partition = PartitionMarginalDense::new(
      0,
      alphabet.clone(),
      &graph,
      &leaf_seq_inputs(&graph, &names, aln.clone()),
    )?;
    let reconstruction = MarginalReconstruction::Dense(DenseReconstruction::seeded(partition, gtr.clone()));
    let (reconstruction, log_lh_before) =
      reconstruction.marginal_update(&graph, &branch_lengths_or_zero(&branch_lengths))?;
    let ab_key = find_node_key_by_name(&graph, &names, "AB").expect("AB exists");

    let reroot = reroot_at_node(
      &mut graph,
      ab_key,
      RerootTopologyParams::default(),
      &mut branch_lengths,
      &names,
    )?;
    let (_, log_lh_rerooted) = reconstruction
      .apply_reroot(&reroot)?
      .marginal_update(&graph, &branch_lengths_or_zero(&branch_lengths))?;

    let fresh_partition = PartitionMarginalDense::new(0, alphabet, &graph, &leaf_seq_inputs(&graph, &names, aln))?;
    let (_, log_lh_fresh) = DenseReconstruction::seeded(fresh_partition, gtr)
      .marginal_update(&graph, &branch_lengths_or_zero(&branch_lengths))?;
    assert!(reroot.edge_merge.is_some(), "the old root is merged away");
    pretty_assert_abs_diff_eq!(log_lh_fresh.value(), log_lh_rerooted.value(), epsilon = 1e-10);
    pretty_assert_abs_diff_eq!(log_lh_before.value(), log_lh_rerooted.value(), epsilon = 1e-10);
    Ok(())
  }
}
