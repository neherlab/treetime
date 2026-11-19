#[cfg(test)]
mod tests {
  use crate::test_utils::TestDataset;
  use eyre::Report;
  use helpers::run_sparse_brent_loop;
  use itertools::Itertools;

  #[test]
  fn test_slow_convergence_sc2_sparse_converges_on_sc2_2844() -> Result<(), Report> {
    let max_iter = 50;
    let result = run_sparse_brent_loop(TestDataset::Sc2_2844, max_iter)?;

    assert!(
      result.stopped_at.is_some(),
      "Sparse optimize loop on sc2/2844 did not converge within {max_iter} iterations. LH history (last 5): {:?}",
      result.lh_history.iter().rev().take(5).collect_vec()
    );
    Ok(())
  }

  #[test]
  fn test_convergence_sc2_sparse_converges_on_flu_h3n2_20() -> Result<(), Report> {
    let max_iter = 10;
    let result = run_sparse_brent_loop(TestDataset::Flu20, max_iter)?;

    assert!(
      result.stopped_at.is_some(),
      "Sparse optimize loop on flu/h3n2/20 did not converge within {max_iter} iterations. LH history (last 5): {:?}",
      result.lh_history.iter().rev().take(5).collect_vec()
    );
    Ok(())
  }

  mod helpers {
    use crate::alphabet::alphabet::Alphabet;
    use crate::branch_lengths::branch_lengths_or_zero;
    use crate::gtr::get_gtr::{JC69Params, jc69};
    use crate::optimize::dispatch::guess_initial_branch_lengths;
    use crate::optimize::gather::InitialGuessEvidence;
    use crate::optimize::indel::IndelPolicy;
    use crate::optimize::params::{BranchOptMethod, ExistingBranchLengths, OptimizeLoopSettings, TopologyOps};
    use crate::optimize::partition_rates::PartitionRates;
    use crate::optimize::run_loop::{OptimizeLoopResult, run_optimize_loop};
    use crate::partition::fitch::passes::create_fitch_partition;
    use crate::partition::marginal::reconstruction::SparseReconstruction;
    use crate::partition::set::PartitionSet;
    use crate::progress::NoopProgress;
    use crate::test_utils::{TestDataset, leaf_seq_inputs};
    use eyre::Report;
    use treetime_io::nwk::nwk_read;

    pub(super) fn run_sparse_brent_loop(dataset: TestDataset, max_iter: usize) -> Result<OptimizeLoopResult, Report> {
      let aln = dataset.alignment()?;
      let nwk_parsed = nwk_read(dataset.tree().as_bytes())?;
      let names = nwk_parsed.names.clone();
      let mut graph = nwk_parsed.graph;
      let mut branch_lengths = nwk_parsed.branch_lengths;

      let fitch = create_fitch_partition(&graph, Alphabet::default(), leaf_seq_inputs(&graph, &names, aln))?;
      let (partition, node_states) = fitch.into_marginal_sparse(&graph)?;
      let reconstruction = PartitionSet::Sparse(vec![SparseReconstruction::seeded(
        partition,
        jc69(JC69Params::default())?,
        node_states,
      )]);
      let (reconstruction, _) = reconstruction.marginal_update(&graph, &branch_lengths_or_zero(&branch_lengths))?;

      let indels = IndelPolicy::Include.edge_indels(&graph, &reconstruction, &branch_lengths)?;
      let evidence = InitialGuessEvidence::gather(&graph, &reconstruction)?;
      guess_initial_branch_lengths(
        &graph,
        &indels,
        &evidence,
        ExistingBranchLengths::Overwrite,
        &mut branch_lengths,
      )?;

      run_optimize_loop(
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
        names,
        &NoopProgress,
      )
    }
  }
}
