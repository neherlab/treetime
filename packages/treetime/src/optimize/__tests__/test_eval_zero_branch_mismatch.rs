#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::ancestral::marginal::{initialize_marginal, update_marginal};
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::optimize::params::BranchOptMethod;
  use crate::optimize::dispatch::run_optimize_mixed;
  use crate::optimize::run_loop::collect_optimize_partitions;
  use crate::ancestral::fitch::create_fitch_partition;
  use crate::partition::marginal_dense::PartitionMarginalDense;
  use crate::seq::alignment::get_common_length;

  use crate::payload::ancestral::GraphAncestral;
  use eyre::Report;
  use indoc::indoc;
  
  use parking_lot::RwLock;
  use std::sync::Arc;
  use treetime_graph::edge::HasBranchLength;
  use treetime_io::fasta::read_many_fasta_str;
  use treetime_io::nwk::nwk_read_str;

  // Regression: run_optimize_mixed must not produce -inf/NaN when entering
  // with branch_length=0 and mismatched certain states. Before the fix,
  // the evaluator computed ln(0) and divided by zero at t=0.
  #[test]
  fn test_eval_zero_branch_mismatch_no_nan() -> Result<(), Report> {
    // Tree with zero-length branches to force the edge case
    let graph: GraphAncestral = nwk_read_str("((A:0.0,B:0.0)AB:0.0,(C:0.0,D:0.0)CD:0.0)root:0.0;")?;

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

    let dense_partitions = vec![Arc::new(RwLock::new(PartitionMarginalDense::new(0, jc69(JC69Params::default())?, alphabet_dense, get_common_length(&aln)?)))];

    let fitch = create_fitch_partition(&graph, 1, alphabet_sparse, &aln)?;
    let sparse_partitions = vec![Arc::new(RwLock::new(
      fitch.into_marginal_sparse(jc69(JC69Params::default())?, &graph)?,
    ))];
    initialize_marginal(&graph, &dense_partitions, &aln)?;
    update_marginal(&graph, &sparse_partitions)?;

    let mixed_partitions = collect_optimize_partitions(&dense_partitions, &sparse_partitions);

    // Do NOT call initial_guess_mixed -- leave branch lengths at 0.0
    // to exercise the zero-branch mismatch code path.
    run_optimize_mixed(&graph, &mixed_partitions, BranchOptMethod::Newton)?;

    // All branch lengths must be finite after optimization
    for edge_ref in graph.get_edges() {
      let bl = edge_ref.read_arc().payload().read_arc().branch_length().unwrap();
      assert!(
        bl.is_finite(),
        "Branch length must be finite after optimization, got {bl}"
      );
      assert!(bl >= 0.0, "Branch length must be non-negative, got {bl}");
    }

    Ok(())
  }
}
