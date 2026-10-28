#[cfg(test)]
mod tests {
  use crate::gtr::jc_distance::jukes_cantor_distance;
  use crate::optimize::indel::IndelPolicy;
  use crate::optimize::params::{BranchOptMethod, OptimizeLoopSettings, TopologyOps};
  use crate::optimize::partition_rates::PartitionRates;
  use crate::optimize::run_loop::run_optimize_loop;
  use crate::progress::NoopProgress;
  use crate::test_utils::find_edge_key;
  use eyre::Report;
  use helpers::shared_substitution_star;
  use indoc::indoc;
  use pretty_assertions::assert_eq;

  #[test]
  fn test_run_optimize_loop_optimizes_the_new_edges_of_a_topology_change_in_the_last_iteration() -> Result<(), Report> {
    let (mut graph, partitions, branch_lengths, names) = shared_substitution_star(indoc! {r#"
      >A
      TTTTACGTACGTACGTACGA
      >B
      TTTTACGTACGTACGTACGC
      >C
      ACGTACGTACGTACGTACGG
      >D
      ACGTACGTACGTACGTACGT
    "#})?;

    let result = run_optimize_loop(
      &mut graph,
      partitions,
      &OptimizeLoopSettings {
        max_iter: 1,
        dp: 0.0,
        damping: 0.0,
        opt_method: BranchOptMethod::Brent,
        indel_policy: IndelPolicy::Include,
        topology_ops: TopologyOps::default(),
        partition_rates: PartitionRates::Estimate,
      },
      branch_lengths,
      names,
      &NoopProgress,
    )?;

    let helper_edge = find_edge_key(&graph, &result.names, "root", "NODE_0000000").unwrap();
    let helper_length = result.branch_lengths[&helper_edge].unwrap();
    let merge_start_length = jukes_cantor_distance(3.0 / 20.0, 4);
    assert_eq!(3, result.lh_history.len());
    assert!(
      (helper_length - merge_start_length).abs() > 1e-6,
      "the new edge keeps its merge start length {merge_start_length}"
    );
    Ok(())
  }

  mod helpers {
    use crate::alphabet::alphabet::{Alphabet, AlphabetName};
    use crate::branch_lengths::branch_lengths_or_zero;
    use crate::gtr::get_gtr::{JC69Params, jc69};
    use crate::optimize::dispatch::guess_initial_branch_lengths;
    use crate::optimize::gather::InitialGuessEvidence;
    use crate::optimize::indel::IndelPolicy;
    use crate::optimize::params::ExistingBranchLengths;
    use crate::partition::fitch::passes::create_fitch_partition;
    use crate::partition::marginal::reconstruction::SparseReconstruction;
    use crate::partition::set::PartitionSet;
    use crate::test_utils::leaf_seq_inputs;
    use eyre::Report;
    use treetime_graph::graph::Graph;
    use treetime_graph::key_map::EdgeMap;
    use treetime_graph::node_names::NodeNames;
    use treetime_io::fasta::fasta_read;
    use treetime_io::nwk::nwk_read;
    use treetime_primitives::AlignmentRecord;

    pub(super) fn shared_substitution_star(
      fasta: &str,
    ) -> Result<(Graph, PartitionSet, EdgeMap<Option<f64>>, NodeNames), Report> {
      let nuc = Alphabet::new(AlphabetName::Nuc)?;
      let aln: Vec<AlignmentRecord> = fasta_read(fasta.as_bytes(), &nuc)?
        .into_iter()
        .map(AlignmentRecord::from)
        .collect();
      let nwk_parsed = nwk_read(b"(A:0.1,B:0.1,C:0.1,D:0.1)root:0.0;".as_slice())?;
      let names = nwk_parsed.names.clone();
      let graph = nwk_parsed.graph;
      let mut branch_lengths = nwk_parsed.branch_lengths;

      let fitch = create_fitch_partition(&graph, nuc, leaf_seq_inputs(&graph, &names, aln))?;
      let (partition, node_states) = fitch.into_marginal_sparse(&graph)?;
      let partitions = PartitionSet::Sparse(vec![SparseReconstruction::seeded(
        partition,
        jc69(JC69Params::default())?,
        node_states,
      )]);
      let (partitions, _) = partitions.marginal_update(&graph, &branch_lengths_or_zero(&branch_lengths))?;
      let indels = IndelPolicy::Include.edge_indels(&graph, &partitions, &branch_lengths)?;
      let evidence = InitialGuessEvidence::gather(&graph, &partitions)?;
      guess_initial_branch_lengths(
        &graph,
        &indels,
        &evidence,
        ExistingBranchLengths::Overwrite,
        &mut branch_lengths,
      )?;
      Ok((graph, partitions, branch_lengths, names))
    }
  }
}
