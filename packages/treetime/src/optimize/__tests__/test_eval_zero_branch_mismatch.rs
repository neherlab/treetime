#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::branch_lengths::branch_lengths_or_zero;
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::optimize::dispatch::run_optimize_mixed;
  use crate::optimize::gather::{gather_edge_contributions, gather_edge_indel_counts};
  use crate::optimize::params::BranchOptMethod;
  use crate::partition::fitch::passes::create_fitch_partition;
  use crate::partition::marginal::reconstruction::{MarginalReconstruction, SparseReconstruction};
  use crate::seq::alignment::node_seq_inputs;
  use treetime_primitives::AlignmentRecord;

  use eyre::Report;
  use indoc::indoc;

  use treetime_graph::graph::Graph;

  use treetime_io::fasta::read_many_fasta_str;
  use treetime_io::nwk::nwk_read_str;

  #[test]
  fn test_eval_zero_branch_mismatch_no_nan() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("((A:0.0,B:0.0)AB:0.0,(C:0.0,D:0.0)CD:0.0)root:0.0;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
    let graph: Graph = graph;

    let alphabet = Alphabet::default();
    let aln: Vec<AlignmentRecord> = read_many_fasta_str(
      indoc! {r#"
        >A
        AAAAAAAAAAAAAAAA
        >B
        CCCCCCCCCCCCCCCC
        >C
        GGGGGGGGGGGGGGGG
        >D
        TTTTTTTTTTTTTTTT
      "#},
      &alphabet,
    )?
    .into_iter()
    .map(AlignmentRecord::from)
    .collect();

    let fitch = create_fitch_partition(
      &graph,
      0,
      Alphabet::new(AlphabetName::Nuc)?,
      &node_seq_inputs(&graph, &names, aln),
    )?;
    let (sparse_partition, sparse_node_states) = fitch.into_marginal_sparse(&graph)?;
    let reconstruction = MarginalReconstruction::Sparse(SparseReconstruction::seeded(
      sparse_partition,
      jc69(JC69Params::default())?,
      sparse_node_states,
    ));
    let (reconstruction, _) = reconstruction.marginal_update(&graph, &branch_lengths_or_zero(&branch_lengths))?;

    let total_length = reconstruction.sequence_length();
    let contributions = gather_edge_contributions(&graph, &reconstruction)?;
    let indel_counts = gather_edge_indel_counts(&graph, &reconstruction);

    run_optimize_mixed(
      &graph,
      total_length,
      &contributions,
      &indel_counts,
      BranchOptMethod::Newton,
      &mut branch_lengths,
    )?;

    for edge_ref in graph.get_edges() {
      let bl = branch_lengths[&edge_ref.key()].unwrap();
      assert!(
        bl.is_finite(),
        "Branch length must be finite after optimization, got {bl}"
      );
      assert!(bl >= 0.0, "Branch length must be non-negative, got {bl}");
    }

    Ok(())
  }
}
