#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::ancestral::fitch::create_fitch_partition;
  use crate::ancestral::marginal::branch_lengths_or_zero;
  use crate::ancestral::pipeline::{DenseReconstruction, SparseReconstruction};
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::optimize::dispatch::run_optimize_mixed;
  use crate::optimize::gather::{gather_edge_contributions, gather_edge_indel_counts, total_sequence_length};
  use crate::optimize::params::BranchOptMethod;
  use crate::optimize::run_loop::{marginal_update_dense, marginal_update_sparse};
  use crate::partition::marginal::dense::partition::PartitionMarginalDense;
  use crate::seq::alignment::get_common_length;
  use treetime_io::nwk::nwk_fasta_node_inputs;

  use eyre::Report;
  use indoc::indoc;

  use treetime_graph::graph::Graph;

  use treetime_io::fasta::read_many_fasta_str;
  use treetime_io::nwk::nwk_read_str;

  // Regression: run_optimize_mixed must not produce -inf/NaN when entering
  // with branch_length=0 and mismatched certain states. Before the fix,
  // the evaluator computed ln(0) and divided by zero at t=0.
  #[test]
  fn test_eval_zero_branch_mismatch_no_nan() -> Result<(), Report> {
    // Tree with zero-length branches to force the edge case
    let nwk_parsed = nwk_read_str("((A:0.0,B:0.0)AB:0.0,(C:0.0,D:0.0)CD:0.0)root:0.0;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
    let graph: Graph = graph;

    // Alignment with mismatches: leaf A differs from leaf B at multiple positions,
    // so after marginal reconstruction some edges have sites where parent and child
    // have disjoint support at t=0.
    let alphabet = Alphabet::default();
    let aln = read_many_fasta_str(
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
    )?;

    let alphabet_dense = Alphabet::new(AlphabetName::Nuc)?;
    let alphabet_sparse = Alphabet::new(AlphabetName::Nuc)?;

    let dense_partition = PartitionMarginalDense::new(0, alphabet_dense, get_common_length(&aln)?);
    let dense_node_states =
      dense_partition.attach_sequences(&graph, &nwk_fasta_node_inputs(&graph, &names, aln.clone()))?;
    let dense_partitions = vec![DenseReconstruction::seeded(
      dense_partition,
      jc69(JC69Params::default())?,
      dense_node_states,
    )];

    let fitch = create_fitch_partition(&graph, 1, alphabet_sparse, &nwk_fasta_node_inputs(&graph, &names, aln))?;
    let (sparse_partition, sparse_node_states) = fitch.into_marginal_sparse(&graph)?;
    let sparse_partitions = vec![SparseReconstruction::seeded(
      sparse_partition,
      jc69(JC69Params::default())?,
      sparse_node_states,
    )];

    let (dense_partitions, _) =
      marginal_update_dense(&graph, &branch_lengths_or_zero(&branch_lengths), dense_partitions)?;
    let (sparse_partitions, _) =
      marginal_update_sparse(&graph, &branch_lengths_or_zero(&branch_lengths), sparse_partitions)?;

    let total_length = total_sequence_length(&dense_partitions, &sparse_partitions);
    let contributions = gather_edge_contributions(&graph, &dense_partitions, &sparse_partitions)?;
    let indel_counts = gather_edge_indel_counts(&graph, &dense_partitions, &sparse_partitions);

    // Do NOT call initial_guess_mixed -- leave branch lengths at 0.0
    // to exercise the zero-branch mismatch code path.
    run_optimize_mixed(
      &graph,
      total_length,
      &contributions,
      &indel_counts,
      BranchOptMethod::Newton,
      &mut branch_lengths,
    )?;

    // All branch lengths must be finite after optimization
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
