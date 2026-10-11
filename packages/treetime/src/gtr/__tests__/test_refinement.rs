#[cfg(test)]
mod tests {
  use crate::gtr::refinement::refine_gtr_model_and_rate;
  use crate::partition::marginal::shared::update::MarginalPasses;
  use crate::progress::NoopProgress;
  use eyre::Report;

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

  mod helpers {
    use crate::alphabet::alphabet::{Alphabet, AlphabetName};
    use crate::branch_lengths::branch_lengths_or_zero;
    use crate::gtr::get_gtr::{JC69Params, jc69};
    use crate::gtr::gtr::GTR;
    use crate::partition::fitch::passes::create_fitch_partition;
    use crate::partition::marginal::sparse::partition::PartitionMarginalSparse;
    use crate::partition::storage::sparse::SparseNodeState;
    use crate::seq::alignment::NodeSeqInput;
    use crate::test_utils::leaf_seq_inputs;
    use eyre::Report;
    use indoc::indoc;
    use std::collections::BTreeMap;
    use treetime_graph::edge::GraphEdgeKey;
    use treetime_graph::graph::Graph;
    use treetime_graph::node::GraphNodeKey;
    use treetime_io::fasta::fasta_read;
    use treetime_io::nwk::nwk_read;
    use treetime_primitives::AlignmentRecord;

    pub(super) struct Input {
      pub(super) graph: Graph,
      pub(super) node_inputs: BTreeMap<GraphNodeKey, NodeSeqInput>,
      pub(super) branch_lengths: BTreeMap<GraphEdgeKey, f64>,
      pub(super) start_gtr: GTR,
    }

    impl Input {
      pub(super) fn new() -> Result<Self, Report> {
        let nwk_parsed = nwk_read(b"((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;".as_slice())?;
        let names = nwk_parsed.names();
        let aln: Vec<AlignmentRecord> = fasta_read(
          indoc! {r#"
          >A
          ACGTACGTACGTACGT
          >B
          ACGTACGAACGTACGA
          >C
          ACCTACGTTCGTACGG
          >D
          ACCTACGTTCGAACGC
        "#}
          .as_bytes(),
          &Alphabet::new(AlphabetName::Nuc)?,
        )?
        .into_iter()
        .map(AlignmentRecord::from)
        .collect();
        Ok(Self {
          node_inputs: leaf_seq_inputs(&nwk_parsed.graph, &names, aln),
          branch_lengths: branch_lengths_or_zero(&nwk_parsed.branch_lengths),
          graph: nwk_parsed.graph,
          start_gtr: jc69(JC69Params::default())?,
        })
      }

      pub(super) fn sparse(
        &self,
      ) -> Result<(PartitionMarginalSparse, BTreeMap<GraphNodeKey, SparseNodeState>), Report> {
        create_fitch_partition(&self.graph, Alphabet::new(AlphabetName::Nuc)?, self.node_inputs.clone())?
          .into_marginal_sparse(&self.graph)
      }
    }
  }
}
