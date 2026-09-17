#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::Alphabet;
  use crate::seq::alignment::node_seq_inputs;
  use treetime_primitives::AlignmentRecord;

  use crate::ancestral::fitch::create_fitch_partition;
  use crate::ancestral::marginal::branch_lengths_or_zero;
  use crate::ancestral::pipeline::SparseReconstruction;
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::optimize::dispatch::initial_guess_mixed;
  use crate::optimize::gather::{
    gather_edge_effective_lengths, gather_edge_indel_counts, gather_edge_sub_counts, total_sequence_length,
  };
  use crate::optimize::params::{BranchOptMethod, TopologyOps};
  use crate::optimize::run_loop::{marginal_update_sparse, run_optimize_loop};

  use eyre::Report;

  use std::path::Path;
  use treetime_io::fasta::read_many_fasta_path;
  use treetime_io::nwk::nwk_read_file;

  /// Regression test: sparse optimize loop converges on sc2/2844 (dataset with indels).
  ///
  /// This dataset exhibits a persistent likelihood 2-cycle without the three-condition convergence check.
  /// The loop alternated between ~-143156 and ~-143157 from iteration ~15 onward,
  /// exhausting all max_iter iterations. With the three-condition convergence check,
  /// damping floor, restored indel rate estimation, and v0-aligned defaults, the loop
  /// must stop via one of the three conditions before exhausting max_iter=50.
  #[test]
  fn test_convergence_sc2_sparse_converges_on_sc2_2844() -> Result<(), Report> {
    let workspace_root = Path::new(env!("CARGO_MANIFEST_DIR"))
      .parent()
      .and_then(|p| p.parent())
      .unwrap();

    let alphabet = Alphabet::default();
    let tree_path = workspace_root.join("data/sc2/2844/tree.nwk");
    let aln_path = workspace_root.join("data/sc2/2844/aln.fasta.xz");
    let aln: Vec<AlignmentRecord> = read_many_fasta_path(&[aln_path.to_str().unwrap()], &alphabet)?
      .into_iter()
      .map(AlignmentRecord::from)
      .collect();
    let nwk_parsed = nwk_read_file(&tree_path)?;
    let names = nwk_parsed.names();
    let mut graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;

    let fitch = create_fitch_partition(&graph, 0, alphabet, &node_seq_inputs(&graph, &names, aln))?;
    let (partition, node_states) = fitch.into_marginal_sparse(&graph)?;
    let sparse_partitions = vec![SparseReconstruction::seeded(
      partition,
      jc69(JC69Params::default())?,
      node_states,
    )];
    let (sparse_partitions, _) =
      marginal_update_sparse(&graph, &branch_lengths_or_zero(&branch_lengths), sparse_partitions)?;

    let dense_partitions = vec![];
    let total_length = total_sequence_length(&dense_partitions, &sparse_partitions);
    let indel_counts = gather_edge_indel_counts(&graph, &dense_partitions, &sparse_partitions);
    let sub_counts = gather_edge_sub_counts(&graph, &dense_partitions, &sparse_partitions)?;
    let effective_lengths = gather_edge_effective_lengths(&graph, &dense_partitions, &sparse_partitions)?;
    initial_guess_mixed(
      &graph,
      total_length,
      &indel_counts,
      &sub_counts,
      &effective_lengths,
      true,
      false,
      &mut branch_lengths,
    )?;

    let max_iter = 50;
    let names_tt_2 = names;
    let result = run_optimize_loop(
      &mut graph,
      sparse_partitions,
      dense_partitions,
      max_iter,
      0.1,
      0.75,
      BranchOptMethod::BrentSqrt,
      false,
      TopologyOps::default(),
      branch_lengths,
      &names_tt_2,
    )?;
    let sparse_partitions = result.sparse_partitions;
    let dense_partitions = result.dense_partitions;

    assert!(
      result.stopped_at.is_some(),
      "Sparse optimize loop on sc2/2844 did not converge within {max_iter} iterations. LH history (last 5): {:?}",
      result.lh_history.iter().rev().take(5).collect::<Vec<_>>()
    );

    Ok(())
  }

  /// Regression test: flu/h3n2/20 converges normally (indel-free dataset).
  ///
  /// Verifies the fix does not regress convergence on a dataset where the
  /// variable position set is stable (no oscillation).
  #[test]
  fn test_convergence_sc2_flu_h3n2_20_converges() -> Result<(), Report> {
    let workspace_root = Path::new(env!("CARGO_MANIFEST_DIR"))
      .parent()
      .and_then(|p| p.parent())
      .unwrap();

    let alphabet = Alphabet::default();
    let tree_path = workspace_root.join("data/flu/h3n2/20/tree.nwk");
    let aln_path = workspace_root.join("data/flu/h3n2/20/aln.fasta.xz");
    let aln: Vec<AlignmentRecord> = read_many_fasta_path(&[aln_path.to_str().unwrap()], &alphabet)?
      .into_iter()
      .map(AlignmentRecord::from)
      .collect();
    let nwk_parsed = nwk_read_file(&tree_path)?;
    let names = nwk_parsed.names();
    let mut graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;

    let fitch = create_fitch_partition(&graph, 0, alphabet, &node_seq_inputs(&graph, &names, aln))?;
    let (partition, node_states) = fitch.into_marginal_sparse(&graph)?;
    let sparse_partitions = vec![SparseReconstruction::seeded(
      partition,
      jc69(JC69Params::default())?,
      node_states,
    )];
    let (sparse_partitions, _) =
      marginal_update_sparse(&graph, &branch_lengths_or_zero(&branch_lengths), sparse_partitions)?;

    let dense_partitions = vec![];
    let total_length = total_sequence_length(&dense_partitions, &sparse_partitions);
    let indel_counts = gather_edge_indel_counts(&graph, &dense_partitions, &sparse_partitions);
    let sub_counts = gather_edge_sub_counts(&graph, &dense_partitions, &sparse_partitions)?;
    let effective_lengths = gather_edge_effective_lengths(&graph, &dense_partitions, &sparse_partitions)?;
    initial_guess_mixed(
      &graph,
      total_length,
      &indel_counts,
      &sub_counts,
      &effective_lengths,
      true,
      false,
      &mut branch_lengths,
    )?;

    let names_tt_1 = names;
    let result = run_optimize_loop(
      &mut graph,
      sparse_partitions,
      dense_partitions,
      10,
      0.1,
      0.75,
      BranchOptMethod::BrentSqrt,
      false,
      TopologyOps::default(),
      branch_lengths,
      &names_tt_1,
    )?;
    let sparse_partitions = result.sparse_partitions;
    let dense_partitions = result.dense_partitions;

    assert!(
      result.stopped_at.is_some(),
      "flu/h3n2/20 did not converge within 10 iterations"
    );
    Ok(())
  }
}
