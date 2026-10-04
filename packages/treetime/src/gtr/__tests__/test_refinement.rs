#[cfg(test)]
mod tests {
  use crate::gtr::refinement::refine_gtr_model_and_rate;
  use crate::partition::marginal::shared::update::MarginalPasses;
  use crate::progress::NoopProgress;
  use eyre::Report;
  use treetime_utils::pretty_assert_abs_diff_eq;

  #[test]
  fn test_refinement_sparse_rate_search_ends_at_log_lh_maximum_in_mu() -> Result<(), Report> {
    let input = helpers::Input::new()?;
    let (partition, node_states) = input.sparse()?;
    let start = partition.marginal_update(&input.start_gtr, &input.graph, &input.branch_lengths, &node_states)?;

    let (gtr, refined) = refine_gtr_model_and_rate(
      &partition,
      input.start_gtr.clone(),
      start,
      0,
      None,
      1.0,
      None,
      &input.graph,
      &input.branch_lengths,
      &NoopProgress,
    )?;

    for factor in [0.99, 1.01] {
      let mut neighbour = gtr.clone();
      neighbour.mu *= factor;
      let at_neighbour = partition.marginal_update(&neighbour, &input.graph, &input.branch_lengths, &node_states)?;
      assert!(
        refined.log_lh.value() >= at_neighbour.log_lh.value(),
        "log likelihood at the optimized mu {} ({}) must not be lower than at mu {} ({})",
        gtr.mu,
        refined.log_lh.value(),
        neighbour.mu,
        at_neighbour.log_lh.value()
      );
    }
    Ok(())
  }

  #[test]
  #[ignore = "sparse transition counts differ from dense: kb/issues/M-gtr-sparse-transition-counts-diverge-from-dense.md"]
  fn test_refinement_sparse_gtr_and_rate_match_dense() -> Result<(), Report> {
    let input = helpers::Input::new()?;
    let (sparse, node_states) = input.sparse()?;
    let dense = input.dense()?;
    let sparse_start = sparse.marginal_update(&input.start_gtr, &input.graph, &input.branch_lengths, &node_states)?;
    let dense_start = dense.marginal_update(&input.start_gtr, &input.graph, &input.branch_lengths, &())?;

    let (sparse_gtr, sparse_refined) = refine_gtr_model_and_rate(
      &sparse,
      input.start_gtr.clone(),
      sparse_start,
      2,
      None,
      1.0,
      None,
      &input.graph,
      &input.branch_lengths,
      &NoopProgress,
    )?;
    let (dense_gtr, dense_refined) = refine_gtr_model_and_rate(
      &dense,
      input.start_gtr.clone(),
      dense_start,
      2,
      None,
      1.0,
      None,
      &input.graph,
      &input.branch_lengths,
      &NoopProgress,
    )?;

    pretty_assert_abs_diff_eq!(dense_gtr.mu, sparse_gtr.mu, epsilon = 1e-7);
    pretty_assert_abs_diff_eq!(dense_gtr.pi, sparse_gtr.pi, epsilon = 1e-7);
    pretty_assert_abs_diff_eq!(dense_gtr.W, sparse_gtr.W, epsilon = 1e-7);
    pretty_assert_abs_diff_eq!(
      dense_refined.log_lh.value(),
      sparse_refined.log_lh.value(),
      epsilon = 1e-7
    );
    Ok(())
  }

  mod helpers {
    use crate::alphabet::alphabet::{Alphabet, AlphabetName};
    use crate::branch_lengths::branch_lengths_or_zero;
    use crate::gtr::get_gtr::{JC69Params, jc69};
    use crate::gtr::gtr::GTR;
    use crate::partition::fitch::passes::create_fitch_partition;
    use crate::partition::marginal::dense::partition::PartitionMarginalDense;
    use crate::partition::marginal::sparse::partition::PartitionMarginalSparse;
    use crate::partition::storage::sparse::SparseNodeState;
    use crate::seq::alignment::{NodeSeqInput, node_seq_inputs};
    use eyre::Report;
    use indoc::indoc;
    use std::collections::BTreeMap;
    use treetime_graph::edge::GraphEdgeKey;
    use treetime_graph::graph::Graph;
    use treetime_graph::node::GraphNodeKey;
    use treetime_io::fasta::read_many_fasta_str;
    use treetime_io::nwk::nwk_read_str;
    use treetime_primitives::AlignmentRecord;

    pub(super) struct Input {
      pub(super) graph: Graph,
      pub(super) node_inputs: BTreeMap<GraphNodeKey, NodeSeqInput>,
      pub(super) branch_lengths: BTreeMap<GraphEdgeKey, f64>,
      pub(super) start_gtr: GTR,
    }

    impl Input {
      pub(super) fn new() -> Result<Self, Report> {
        let nwk_parsed = nwk_read_str("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?;
        let names = nwk_parsed.names();
        let aln: Vec<AlignmentRecord> = read_many_fasta_str(
          indoc! {r#"
          >A
          ACGTACGTACGTACGT
          >B
          ACGTACGAACGTACGA
          >C
          ACCTACGTTCGTACGG
          >D
          ACCTACGTTCGAACGC
        "#},
          &Alphabet::new(AlphabetName::Nuc)?,
        )?
        .into_iter()
        .map(AlignmentRecord::from)
        .collect();
        Ok(Self {
          node_inputs: node_seq_inputs(&nwk_parsed.graph, &names, aln),
          branch_lengths: branch_lengths_or_zero(&nwk_parsed.branch_lengths),
          graph: nwk_parsed.graph,
          start_gtr: jc69(JC69Params::default())?,
        })
      }

      pub(super) fn sparse(
        &self,
      ) -> Result<(PartitionMarginalSparse, BTreeMap<GraphNodeKey, SparseNodeState>), Report> {
        create_fitch_partition(
          &self.graph,
          0,
          Alphabet::new(AlphabetName::Nuc)?,
          self.node_inputs.clone(),
        )?
        .into_marginal_sparse(&self.graph)
      }

      pub(super) fn dense(&self) -> Result<PartitionMarginalDense, Report> {
        PartitionMarginalDense::new(0, Alphabet::new(AlphabetName::Nuc)?, &self.graph, &self.node_inputs)
      }
    }
  }
}
