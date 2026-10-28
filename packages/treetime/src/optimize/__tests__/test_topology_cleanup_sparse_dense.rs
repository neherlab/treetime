#[cfg(test)]
mod tests {
  use eyre::Report;
  use helpers::{
    Case, HOIST_CLEAN, HOIST_GAPPED_REVERTING_CHILD, HOIST_N_RUN_IN_SIBLING, MERGE_CLEAN, MERGE_GAPPED_CHILD,
    MERGE_N_RUN_IN_SIBLING, MERGE_SHARED_DELETION, cleaned_sparse_summary, dense_summary,
  };
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime_utils::pretty_assert_map_abs_diff_eq;

  #[rustfmt::skip]
  #[rstest]
  #[case::merge_clean(                  MERGE_CLEAN,                  "((C,D,(A,B)NODE_0000000)P,E)root;"          )]
  #[case::merge_gapped_child(           MERGE_GAPPED_CHILD,           "((C,D,(A,B)NODE_0000000)P,E)root;"          )]
  #[case::merge_shared_deletion(        MERGE_SHARED_DELETION,        "((C,D,(A,B)NODE_0000000)P,E)root;"          )]
  #[case::merge_n_run_in_sibling(       MERGE_N_RUN_IN_SIBLING,       "((C,D,(A,B)NODE_0000000)P,E)root;"          )]
  #[case::hoist_clean(                  HOIST_CLEAN,                  "((S,((C2,C3,C4)V,C1)NODE_0000000)U,T)root;")]
  #[case::hoist_gapped_reverting_child( HOIST_GAPPED_REVERTING_CHILD, "((S,((C2,C3,C4)V,C1)NODE_0000000)U,T)root;")]
  #[case::hoist_n_run_in_sibling(       HOIST_N_RUN_IN_SIBLING,       "((S,((C2,C3,C4)V,C1)NODE_0000000)U,T)root;")]
  #[trace]
  fn test_topology_cleanup_sparse_matches_dense_on_the_nodes_that_merge_and_hoist_create(
    #[case] case: Case,
    #[case] expected_topology: &str,
  ) -> Result<(), Report> {
    let sparse = cleaned_sparse_summary(case)?;
    let dense = dense_summary(case, &sparse.newick)?;

    assert_eq!(
      (expected_topology, vec![]),
      (sparse.topology.as_str(), dense.sites_resolved_above_folding_threshold)
    );
    pretty_assert_map_abs_diff_eq!(dense.values, sparse.values, epsilon = 1e-10);
    Ok(())
  }

  mod helpers {
    use crate::gtr::get_gtr::{JC69Params, jc69};
    use crate::optimize::indel::IndelPolicy;
    use crate::optimize::likelihood::evaluate_partitions_log_lh_only;
    use crate::optimize::params::TopologyOps;
    use crate::optimize::run_loop::prune_and_merge_in_loop;
    use crate::optimize::topology::merge_shared_mutations::MergeLedger;
    use crate::partition::marginal::reconstruction::ReconstructionRef;
    use crate::partition::set::PartitionSet;
    use crate::test_utils::{NUC_ALPHABET, run_dense_marginal, run_sparse_marginal};
    use eyre::Report;
    use itertools::{Itertools, izip};
    use std::collections::BTreeMap;
    use std::iter::once;
    use treetime_graph::graph::Graph;
    use treetime_graph::key_map::EdgeMap;
    use treetime_graph::node_names::NodeNames;
    use treetime_graph::schedule::GraphSchedule;
    use treetime_io::fasta::fasta_read;
    use treetime_io::nwk::{NwkNodeComments, NwkParse, NwkWriteOptions, nwk_read, nwk_write_str};
    use treetime_primitives::AlignmentRecord;

    const MERGE_TREE: &str = "((A:0.3,B:0.25,C:0.35,D:0.3)P:0.2,E:0.4)root;";
    const MERGE_NAMES: [&str; 5] = ["A", "B", "C", "D", "E"];
    const HOIST_TREE: &str = "(((C1:0.2,C2:0.25,C3:0.3,C4:0.2)V:0.2,S:0.3)U:0.15,T:0.35)root;";
    const HOIST_NAMES: [&str; 6] = ["C1", "C2", "C3", "S", "T", "C4"];
    const BRANCH_LENGTHS: [f64; 3] = [0.01, 0.1, 1.0];
    const CLEANED_LENGTH: f64 = 0.3;
    const FOLDING_THRESHOLD: f64 = 1.0 - 1e-4;

    #[rustfmt::skip]
    pub(crate) const MERGE_CLEAN:                  Case = Case::merge(&["ACGGACGTACGTACGTACGA", "ACGGACGTACCTACGTACGT", "ACGTACGTACGTTCGTACGT", "ACGTACGAACGTACGTACGT", "ACGTACGTACGTACGTCCGT"]);
    #[rustfmt::skip]
    pub(crate) const MERGE_GAPPED_CHILD:           Case = Case::merge(&["ACGG------GTACGTACGA", "ACGGACGTACCTACGTACGT", "ACGTACGTACGTTCGTACGT", "ACGTACGAACGTACGTACGT", "ACGTACGTACGTACGTCCGT"]);
    #[rustfmt::skip]
    pub(crate) const MERGE_SHARED_DELETION:        Case = Case::merge(&["ACGG------GTACGTACGA", "ACGG------CTACGTACGT", "ACGTACGTACGTTCGTACGT", "ACGTACGAACGTACGTACGT", "ACGTACGTACGTACGTCCGT"]);
    #[rustfmt::skip]
    pub(crate) const MERGE_N_RUN_IN_SIBLING:       Case = Case::merge(&["ACGGACGTACGTACGTACGA", "ACGGACGTACCTACGTACGT", "ACGTNNNNNNGTTCGTACGT", "ACGTACGAACGTACGTACGT", "ACGTACGTACGTACGTCCGT"]);
    #[rustfmt::skip]
    pub(crate) const HOIST_CLEAN:                  Case = Case::hoist(&["ACGTACGTAGGTACGTACGA", "ACGTTCGTAGGTACATACGT", "ATGTTCGTAGGTACGTACGT", "ACGTACGTACGTACGTAGGT", "ACGTACGTACGTCCGTACGT", "ACGTTCCTAGGTACGTACGT"]);
    #[rustfmt::skip]
    pub(crate) const HOIST_GAPPED_REVERTING_CHILD: Case = Case::hoist(&["ACGTACGTAGGT-----CGA", "ACGTTCGTAGGTACATACGT", "ATGTTCGTAGGTACGTACGT", "ACGTACGTACGTACGTAGGT", "ACGTACGTACGTCCGTACGT", "ACGTTCCTAGGTACGTACGT"]);
    #[rustfmt::skip]
    pub(crate) const HOIST_N_RUN_IN_SIBLING:       Case = Case::hoist(&["ACGTACGTAGGTACGTACGA", "ACGTTCGTAGGTACATACGT", "ATGTTCGTAGGTACGTACGT", "ACGTACGTACGTNNNNNGGT", "ACGTACGTACGTCCGTACGT", "ACGTTCCTAGGTACGTACGT"]);

    #[derive(Clone, Copy, Debug)]
    pub(crate) struct Case {
      tree: &'static str,
      names: &'static [&'static str],
      seqs: &'static [&'static str],
    }

    impl Case {
      const fn merge(seqs: &'static [&'static str]) -> Self {
        Self {
          tree: MERGE_TREE,
          names: &MERGE_NAMES,
          seqs,
        }
      }

      const fn hoist(seqs: &'static [&'static str]) -> Self {
        Self {
          tree: HOIST_TREE,
          names: &HOIST_NAMES,
          seqs,
        }
      }

      fn alignment(self) -> Result<Vec<AlignmentRecord>, Report> {
        let fasta = izip!(self.names, self.seqs)
          .map(|(name, seq)| format!(">{name}\n{seq}\n"))
          .join("");
        Ok(
          fasta_read(fasta.as_bytes(), &*NUC_ALPHABET)?
            .into_iter()
            .map(AlignmentRecord::from)
            .collect(),
        )
      }
    }

    pub(crate) struct CleanedSparse {
      pub(crate) topology: String,
      pub(crate) newick: String,
      pub(crate) values: BTreeMap<String, f64>,
    }

    pub(crate) fn cleaned_sparse_summary(case: Case) -> Result<CleanedSparse, Report> {
      let NwkParse {
        mut graph,
        mut branch_lengths,
        mut names,
        ..
      } = nwk_read(case.tree.as_bytes())?;
      let (_, sparse) = run_sparse_marginal(
        &graph,
        &branch_lengths,
        &names,
        &case.alignment()?,
        jc69(JC69Params::default())?,
      )?;
      let cleanup = prune_and_merge_in_loop(
        &mut graph,
        PartitionSet::Sparse(vec![sparse]),
        &[],
        TopologyOps::default(),
        &mut branch_lengths,
        &mut names,
        &mut MergeLedger::new(IndelPolicy::Include),
      )?;
      let branch_lengths: EdgeMap<Option<f64>> = graph
        .get_edges()
        .map(|edge| (edge.key(), Some(CLEANED_LENGTH)))
        .collect();
      let lengths = branch_lengths
        .iter()
        .map(|(key, length)| (key, length.unwrap_or(0.0)))
        .collect();
      let (partitions, log_lh) = cleanup.partitions.marginal_update(&graph, &lengths)?;
      let reconstruction = partitions.iter().next().unwrap();
      let no_lengths: EdgeMap<Option<f64>> = graph.get_edges().map(|edge| (edge.key(), None)).collect();
      Ok(CleanedSparse {
        topology: write_newick(&graph, &names, &no_lengths)?,
        newick: write_newick(&graph, &names, &branch_lengths)?,
        values: summary(&graph, &names, &lengths, reconstruction, log_lh.value())?,
      })
    }

    pub(crate) struct DenseSummary {
      pub(crate) values: BTreeMap<String, f64>,
      pub(crate) sites_resolved_above_folding_threshold: Vec<(String, usize)>,
    }

    pub(crate) fn dense_summary(case: Case, newick: &str) -> Result<DenseSummary, Report> {
      let NwkParse {
        graph,
        branch_lengths,
        names,
        ..
      } = nwk_read(newick.as_bytes())?;
      let (log_lh, dense) = run_dense_marginal(
        &graph,
        &branch_lengths,
        &names,
        &case.alignment()?,
        jc69(JC69Params::default())?,
      )?;
      let lengths = branch_lengths
        .iter()
        .map(|(key, length)| (key, length.unwrap_or(0.0)))
        .collect();
      let sites_resolved_above_folding_threshold = graph
        .get_internal_nodes()
        .flat_map(|node| {
          let name = &names[node.key()];
          dense.node_states[&node.key()]
            .profile
            .dis
            .rows()
            .into_iter()
            .enumerate()
            .filter(|(_, row)| row.iter().any(|&p| p >= FOLDING_THRESHOLD))
            .map(|(pos, _)| (name.to_owned(), pos))
            .collect_vec()
        })
        .collect();
      Ok(DenseSummary {
        values: summary(&graph, &names, &lengths, ReconstructionRef::Dense(&dense), log_lh)?,
        sites_resolved_above_folding_threshold,
      })
    }

    fn summary(
      graph: &Graph,
      names: &NodeNames,
      lengths: &EdgeMap<f64>,
      reconstruction: ReconstructionRef<'_>,
      log_lh: f64,
    ) -> Result<BTreeMap<String, f64>, Report> {
      let changes = graph
        .get_edges()
        .map(|edge| edge.key())
        .cartesian_product(BRANCH_LENGTHS)
        .map(|(edge_key, length)| -> Result<_, Report> {
          let child = &names[graph.edge_endpoints(edge_key)?.1];
          let contribution = [reconstruction.create_edge_contribution(edge_key)?];
          let change = evaluate_partitions_log_lh_only(&contribution, length)?.value()
            - evaluate_partitions_log_lh_only(&contribution, lengths[&edge_key])?.value();
          Ok((format!("edge {child} t={length} objective change"), change))
        })
        .collect::<Result<Vec<_>, Report>>()?;
      Ok(once((String::from("root log_lh"), log_lh)).chain(changes).collect())
    }

    fn write_newick(graph: &Graph, names: &NodeNames, lengths: &EdgeMap<Option<f64>>) -> Result<String, Report> {
      nwk_write_str(
        &GraphSchedule::new(graph)?.into_tree()?,
        names,
        lengths,
        &NwkWriteOptions {
          weight_significant_digits: Some(17),
          ..NwkWriteOptions::default()
        },
        &NwkNodeComments::new(),
      )
    }
  }
}
