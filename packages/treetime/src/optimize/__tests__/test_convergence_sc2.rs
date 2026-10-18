#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::Alphabet;
  use crate::branch_lengths::branch_lengths_or_zero;
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::optimize::dispatch::guess_initial_branch_lengths;
  use crate::optimize::gather::InitialGuessEvidence;
  use crate::optimize::indel::IndelPolicy;
  use crate::optimize::params::ExistingBranchLengths;
  use crate::optimize::params::OptimizeLoopSettings;
  use crate::optimize::params::{BranchOptMethod, TopologyOps};
  use crate::optimize::partition_rates::PartitionRates;
  use crate::optimize::run_loop::run_optimize_loop;
  use crate::partition::fitch::passes::create_fitch_partition;
  use crate::partition::marginal::reconstruction::SparseReconstruction;
  use crate::partition::set::PartitionSet;
  use crate::progress::NoopProgress;
  use crate::test_utils::leaf_seq_inputs;
  use eyre::Report;
  use std::path::Path;
  use treetime_io::fasta::fasta_read_file;
  use treetime_io::nwk::nwk_read_file;
  use treetime_primitives::AlignmentRecord;

  #[test]
  fn test_convergence_sc2_sparse_converges_on_sc2_2844() -> Result<(), Report> {
    let workspace_root = Path::new(env!("CARGO_MANIFEST_DIR"))
      .parent()
      .and_then(|p| p.parent())
      .unwrap();

    let alphabet = Alphabet::default();
    let tree_path = workspace_root.join("data/sc2/2844/tree.nwk");
    let aln_path = workspace_root.join("data/sc2/2844/aln.fasta.xz");
    let aln: Vec<AlignmentRecord> = fasta_read_file(&aln_path, &alphabet)?
      .into_iter()
      .map(AlignmentRecord::from)
      .collect();
    let nwk_parsed = nwk_read_file(&tree_path)?;
    let names = nwk_parsed.names();
    let mut graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;

    let fitch = create_fitch_partition(&graph, alphabet, leaf_seq_inputs(&graph, &names, aln))?;
    let (partition, node_states) = fitch.into_marginal_sparse(&graph)?;
    let reconstruction = PartitionSet::Sparse(vec![SparseReconstruction::seeded(
      partition,
      jc69(JC69Params::default())?,
      node_states,
    )]);
    let (reconstruction, _) = reconstruction.marginal_update(&graph, &branch_lengths_or_zero(&branch_lengths))?;

    let indels = IndelPolicy::Include.edge_indels(&graph, &reconstruction, &branch_lengths);
    let evidence = InitialGuessEvidence::gather(&graph, &reconstruction)?;
    guess_initial_branch_lengths(
      &graph,
      &indels,
      &evidence,
      ExistingBranchLengths::Overwrite,
      &mut branch_lengths,
    )?;

    let max_iter = 50;
    let names_tt_2 = names;
    let result = run_optimize_loop(
      &mut graph,
      reconstruction,
      &OptimizeLoopSettings {
        max_iter,
        dp: 0.1,
        damping: 0.75,
        opt_method: BranchOptMethod::Brent,
        indel_policy: IndelPolicy::Include,
        topology_ops: TopologyOps::default(),
        partition_rates: PartitionRates::Estimate,
      },
      branch_lengths,
      &names_tt_2,
      &NoopProgress,
    )?;

    assert!(
      result.stopped_at.is_some(),
      "Sparse optimize loop on sc2/2844 did not converge within {max_iter} iterations. LH history (last 5): {:?}",
      result.lh_history.iter().rev().take(5).collect::<Vec<_>>()
    );

    Ok(())
  }

  #[test]
  fn test_convergence_sc2_flu_h3n2_20_converges() -> Result<(), Report> {
    let workspace_root = Path::new(env!("CARGO_MANIFEST_DIR"))
      .parent()
      .and_then(|p| p.parent())
      .unwrap();

    let alphabet = Alphabet::default();
    let tree_path = workspace_root.join("data/flu/h3n2/20/tree.nwk");
    let aln_path = workspace_root.join("data/flu/h3n2/20/aln.fasta.xz");
    let aln: Vec<AlignmentRecord> = fasta_read_file(&aln_path, &alphabet)?
      .into_iter()
      .map(AlignmentRecord::from)
      .collect();
    let nwk_parsed = nwk_read_file(&tree_path)?;
    let names = nwk_parsed.names();
    let mut graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;

    let fitch = create_fitch_partition(&graph, alphabet, leaf_seq_inputs(&graph, &names, aln))?;
    let (partition, node_states) = fitch.into_marginal_sparse(&graph)?;
    let reconstruction = PartitionSet::Sparse(vec![SparseReconstruction::seeded(
      partition,
      jc69(JC69Params::default())?,
      node_states,
    )]);
    let (reconstruction, _) = reconstruction.marginal_update(&graph, &branch_lengths_or_zero(&branch_lengths))?;

    let indels = IndelPolicy::Include.edge_indels(&graph, &reconstruction, &branch_lengths);
    let evidence = InitialGuessEvidence::gather(&graph, &reconstruction)?;
    guess_initial_branch_lengths(
      &graph,
      &indels,
      &evidence,
      ExistingBranchLengths::Overwrite,
      &mut branch_lengths,
    )?;

    let names_tt_1 = names;
    let result = run_optimize_loop(
      &mut graph,
      reconstruction,
      &OptimizeLoopSettings {
        max_iter: 10,
        dp: 0.1,
        damping: 0.75,
        opt_method: BranchOptMethod::Brent,
        indel_policy: IndelPolicy::Include,
        topology_ops: TopologyOps::default(),
        partition_rates: PartitionRates::Estimate,
      },
      branch_lengths,
      &names_tt_1,
      &NoopProgress,
    )?;

    assert!(
      result.stopped_at.is_some(),
      "flu/h3n2/20 did not converge within 10 iterations"
    );
    Ok(())
  }
}
