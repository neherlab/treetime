#[cfg(test)]
mod tests {
  use crate::gtr::refinement::{RateRangeEdge, RefinementRate, optimize_gtr_rate, refine_gtr};
  use crate::partition::marginal::shared::update::MarginalPasses;
  use crate::progress::NoopProgress;
  use eyre::Report;
  use rstest::rstest;
  use std::collections::BTreeMap;
  use treetime_utils::{pretty_assert_abs_diff_eq, pretty_assert_eq, pretty_assert_ulps_eq};

  #[rustfmt::skip]
  #[rstest]
  #[case::keep_rate(        RefinementRate::Keep                           )]
  #[case::optimize_rate(    RefinementRate::Optimize { scale: None      } )]
  #[case::optimize_and_scale(RefinementRate::Optimize { scale: Some(0.5) } )]
  #[trace]
  fn test_refinement_returns_the_update_under_the_refined_model(#[case] rate: RefinementRate) -> Result<(), Report> {
    let input = helpers::Input::new()?;
    let (partition, node_states) = input.sparse()?;
    let start = partition.marginal_update(&input.start_gtr, &input.graph, &input.branch_lengths, &node_states)?;

    let (gtr, refined) = refine_gtr(
      &partition,
      input.start_gtr.clone(),
      start,
      &helpers::refinement(2, rate),
      &input.graph,
      &input.branch_lengths,
      &NoopProgress,
    )?;

    let fresh = partition.marginal_update(&gtr, &input.graph, &input.branch_lengths, &node_states)?;
    pretty_assert_abs_diff_eq!(fresh.log_lh.value(), refined.log_lh.value(), epsilon = 1e-10);
    Ok(())
  }

  #[test]
  fn test_refinement_scale_multiplies_the_optimized_rate_only() -> Result<(), Report> {
    let input = helpers::Input::new()?;
    let (partition, node_states) = input.sparse()?;
    let start = partition.marginal_update(&input.start_gtr, &input.graph, &input.branch_lengths, &node_states)?;
    let refine = |rate| {
      refine_gtr(
        &partition,
        input.start_gtr.clone(),
        start.clone(),
        &helpers::refinement(1, rate),
        &input.graph,
        &input.branch_lengths,
        &NoopProgress,
      )
    };

    let (unscaled, _) = refine(RefinementRate::Optimize { scale: None })?;
    let (scaled, _) = refine(RefinementRate::Optimize { scale: Some(0.5) })?;

    pretty_assert_ulps_eq!(0.5 * unscaled.mu(), scaled.mu(), max_ulps = 4);
    pretty_assert_eq!(unscaled.pi(), scaled.pi());
    pretty_assert_eq!(unscaled.W(), scaled.W());
    Ok(())
  }

  #[test]
  fn test_refinement_sparse_rate_search_ends_at_log_lh_maximum_in_mu() -> Result<(), Report> {
    let input = helpers::Input::new()?;
    let (partition, node_states) = input.sparse()?;
    let start = partition.marginal_update(&input.start_gtr, &input.graph, &input.branch_lengths, &node_states)?;

    let (gtr, refined) = refine_gtr(
      &partition,
      input.start_gtr.clone(),
      start,
      &helpers::refinement(0, RefinementRate::Optimize { scale: None }),
      &input.graph,
      &input.branch_lengths,
      &NoopProgress,
    )?;

    for factor in [0.99, 1.01] {
      let mut neighbour = gtr.clone();
      neighbour.set_mu(neighbour.mu() * factor)?;
      let at_neighbour = partition.marginal_update(&neighbour, &input.graph, &input.branch_lengths, &node_states)?;
      assert!(
        refined.log_lh.value() >= at_neighbour.log_lh.value(),
        "log likelihood at the optimized mu {} ({}) must not be lower than at mu {} ({})",
        gtr.mu(),
        refined.log_lh.value(),
        neighbour.mu(),
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

    let refinement = helpers::refinement(2, RefinementRate::Optimize { scale: None });

    let (sparse_gtr, sparse_refined) = refine_gtr(
      &sparse,
      input.start_gtr.clone(),
      sparse_start,
      &refinement,
      &input.graph,
      &input.branch_lengths,
      &NoopProgress,
    )?;
    let (dense_gtr, dense_refined) = refine_gtr(
      &dense,
      input.start_gtr.clone(),
      dense_start,
      &refinement,
      &input.graph,
      &input.branch_lengths,
      &NoopProgress,
    )?;

    pretty_assert_abs_diff_eq!(dense_gtr.mu(), sparse_gtr.mu(), epsilon = 1e-7);
    pretty_assert_abs_diff_eq!(dense_gtr.pi(), sparse_gtr.pi(), epsilon = 1e-7);
    pretty_assert_abs_diff_eq!(dense_gtr.W(), sparse_gtr.W(), epsilon = 1e-7);
    pretty_assert_abs_diff_eq!(
      dense_refined.log_lh.value(),
      sparse_refined.log_lh.value(),
      epsilon = 1e-7
    );
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::sparse(false)]
  #[case::dense( true )]
  #[trace]
  fn test_refinement_rate_step_on_an_alignment_without_change_takes_the_lower_range_end_with_its_messages(
    #[case] dense: bool,
  ) -> Result<(), Report> {
    let input = helpers::Input::without_change()?;
    let (gtr, log_lh, fresh_log_lh) = if dense {
      let partition = input.dense()?;
      let (gtr, nodes, _) = optimize_gtr_rate(
        &partition,
        input.start_gtr.clone(),
        &BTreeMap::new(),
        &input.graph,
        &input.branch_lengths,
        RateRangeEdge::LowerCostEnd,
        &NoopProgress,
      )?;
      let fresh = partition.marginal_update(&gtr, &input.graph, &input.branch_lengths, &())?;
      let log_lh = partition.root_log_lh(&input.graph, &nodes)?;
      (gtr, log_lh, fresh.log_lh)
    } else {
      let (partition, node_states) = input.sparse()?;
      let (gtr, nodes, _) = optimize_gtr_rate(
        &partition,
        input.start_gtr.clone(),
        &node_states,
        &input.graph,
        &input.branch_lengths,
        RateRangeEdge::LowerCostEnd,
        &NoopProgress,
      )?;
      let fresh = partition.marginal_update(&gtr, &input.graph, &input.branch_lengths, &node_states)?;
      let log_lh = partition.root_log_lh(&input.graph, &nodes)?;
      (gtr, log_lh, fresh.log_lh)
    };
    pretty_assert_ulps_eq!(1e-4 * input.start_gtr.mu(), gtr.mu(), max_ulps = 4);
    pretty_assert_abs_diff_eq!(fresh_log_lh.value(), log_lh.value(), epsilon = 1e-10);
    Ok(())
  }

  mod helpers {
    use crate::alphabet::alphabet::{Alphabet, AlphabetName};
    use crate::branch_lengths::branch_lengths_or_zero;
    use crate::gtr::get_gtr::{JC69Params, jc69};
    use crate::gtr::gtr::GTR;
    use crate::gtr::refinement::{GtrRefinement, RefinementRate};
    use crate::partition::fitch::passes::create_fitch_partition;
    use crate::partition::marginal::dense::partition::PartitionMarginalDense;
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

    pub(super) const fn refinement(iterations: usize, rate: RefinementRate) -> GtrRefinement {
      GtrRefinement {
        iterations,
        pc: 1.0,
        fixed_pi: None,
        rate,
      }
    }

    pub(super) struct Input {
      pub(super) graph: Graph,
      pub(super) node_inputs: BTreeMap<GraphNodeKey, NodeSeqInput>,
      pub(super) branch_lengths: BTreeMap<GraphEdgeKey, f64>,
      pub(super) start_gtr: GTR,
    }

    impl Input {
      pub(super) fn new() -> Result<Self, Report> {
        Self::from_fasta(indoc! {r#"
          >A
          ACGTACGTACGTACGT
          >B
          ACGTACGAACGTACGA
          >C
          ACCTACGTTCGTACGG
          >D
          ACCTACGTTCGAACGC
        "#})
      }

      pub(super) fn without_change() -> Result<Self, Report> {
        Self::from_fasta(indoc! {r#"
          >A
          ACGTACGTACGTACGT
          >B
          ACGTACGTACGTACGT
          >C
          ACGTACGTACGTACGT
          >D
          ACGTACGTACGTACGT
        "#})
      }

      fn from_fasta(fasta: &str) -> Result<Self, Report> {
        let nwk_parsed = nwk_read(b"((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;".as_slice())?;
        let names = nwk_parsed.names();
        let aln: Vec<AlignmentRecord> = fasta_read(fasta.as_bytes(), &Alphabet::new(AlphabetName::Nuc)?)?
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

      pub(super) fn dense(&self) -> Result<PartitionMarginalDense, Report> {
        PartitionMarginalDense::new(Alphabet::new(AlphabetName::Nuc)?, &self.graph, &self.node_inputs)
      }
    }
  }
}
