#[cfg(test)]
mod tests {
  use crate::partition::create::Representation;
  use crate::partition::set::PartitionSet;
  use crate::test_utils::{TestDataset, split_at_gap_free_column};
  use eyre::Report;
  use helpers::updated_set;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime_utils::pretty_assert_abs_diff_eq;

  #[rustfmt::skip]
  #[rstest]
  #[case::dense( Representation::Dense )]
  #[case::sparse(Representation::Sparse)]
  #[trace]
  fn test_partition_set_of_two_column_blocks_sums_to_the_whole_alignment(
    #[case] representation: Representation,
  ) -> Result<(), Report> {
    let whole = TestDataset::Flu20.alignment()?;
    let [first, second] = split_at_gap_free_column(&whole)?;
    let (graph, whole) = updated_set(representation, vec![whole])?;
    let (_, split) = updated_set(representation, vec![first, second])?;

    pretty_assert_abs_diff_eq!(
      whole.graph_log_lh(&graph)?.value(),
      split.graph_log_lh(&graph)?.value(),
      epsilon = 1e-8
    );
    let sub_counts = |set: &PartitionSet| -> Result<Vec<Option<usize>>, Report> {
      graph.get_edges().map(|edge| set.edge_sub_count(&graph, edge.key())).collect()
    };
    assert_eq!(
      (whole.total_length(), sub_counts(&whole)?),
      (split.total_length(), sub_counts(&split)?)
    );
    Ok(())
  }

  mod helpers {
    use crate::alphabet::alphabet::Alphabet;
    use crate::branch_lengths::branch_lengths_or_zero;
    use crate::gtr::get_gtr::GtrModelName;
    use crate::partition::create::Representation;
    use crate::partition::set::{PartitionSet, build_partition_set};
    use crate::progress::NoopProgress;
    use crate::test_utils::{TestDataset, leaf_seq_inputs, sequence_source};
    use eyre::Report;
    use treetime_graph::graph::Graph;
    use treetime_io::nwk::nwk_read;
    use treetime_primitives::AlignmentRecord;

    pub(super) fn updated_set(
      representation: Representation,
      files: Vec<Vec<AlignmentRecord>>,
    ) -> Result<(Graph, PartitionSet), Report> {
      let parse = nwk_read(TestDataset::Flu20.tree().as_bytes())?;
      let names = parse.names();
      let sources = files
        .into_iter()
        .map(|records| sequence_source(Alphabet::default(), leaf_seq_inputs(&parse.graph, &names, records)))
        .collect::<Result<Vec<_>, Report>>()?;
      let branch_lengths = branch_lengths_or_zero(&parse.branch_lengths);
      let (_, set) = build_partition_set(
        representation,
        GtrModelName::JC69,
        &parse.graph,
        sources,
        &branch_lengths,
        &NoopProgress,
      )?;
      let (set, _) = set.marginal_update(&parse.graph, &branch_lengths)?;
      Ok((parse.graph, set))
    }
  }
}
