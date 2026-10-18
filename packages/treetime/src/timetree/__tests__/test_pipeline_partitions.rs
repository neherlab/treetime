#[cfg(test)]
mod tests {
  use crate::optimize::partition_rates::PartitionRates;
  use crate::test_utils::{Columns, TestDataset, split_at_gap_free_column};
  use eyre::Report;
  use helpers::timetree;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime_utils::{pretty_assert_abs_diff_eq, pretty_assert_map_abs_diff_eq};

  #[rustfmt::skip]
  #[rstest]
  #[case::zika_dense(  TestDataset::Zika20,  true,  Columns::Leading(200))]
  #[case::zika_sparse( TestDataset::Zika20,  false, Columns::All          )]
  #[case::flu_dense(   TestDataset::Flu20,   true,  Columns::Leading(200))]
  #[case::flu_sparse(  TestDataset::Flu20,   false, Columns::All          )]
  #[case::ebola_dense( TestDataset::Ebola20, true,  Columns::Leading(200))]
  #[case::ebola_sparse(TestDataset::Ebola20, false, Columns::All          )]
  #[trace]
  fn test_timetree_two_column_blocks_give_the_dates_of_one_alignment(
    #[case] dataset: TestDataset,
    #[case] dense: bool,
    #[case] selection: Columns,
  ) -> Result<(), Report> {
    let records = dataset.alignment_columns(selection)?;
    let [first, second] = split_at_gap_free_column(&records)?;
    let whole = timetree(dataset, vec![records], dense, PartitionRates::Equal)?;
    let split = timetree(dataset, vec![first, second], dense, PartitionRates::Equal)?;

    assert_eq!(
      (whole.edges, whole.total_length, vec![1, 2]),
      (split.edges, split.total_length, split.numbers)
    );
    pretty_assert_map_abs_diff_eq!(whole.node_dates, split.node_dates, epsilon = 1e-8);
    pretty_assert_abs_diff_eq!(whole.log_lh, split.log_lh, epsilon = 1e-8);
    Ok(())
  }

  mod helpers {
    use crate::alphabet::alphabet::Alphabet;
    use crate::cancel::NoopCancel;
    use crate::error::OperationError;
    use crate::optimize::partition_rates::PartitionRates;
    use crate::progress::NoopProgress;
    use crate::test_utils::{TestDataset, dates_by_node, leaf_seq_inputs, marginal_timetree_params, sequence_source};
    use crate::timetree::params::TimetreeParams;
    use crate::timetree::pipeline::{TimetreeInput, run};
    use eyre::Report;
    use std::collections::{BTreeMap, BTreeSet};
    use treetime_graph::node::GraphNodeKey;
    use treetime_io::nwk::nwk_read;
    use treetime_primitives::{AlignmentRecord, LogLh};

    pub(super) struct Timetree {
      pub(super) edges: BTreeSet<(Option<String>, Option<String>)>,
      pub(super) node_dates: BTreeMap<GraphNodeKey, f64>,
      pub(super) total_length: usize,
      pub(super) numbers: Vec<usize>,
      pub(super) log_lh: f64,
    }

    pub(super) fn timetree(
      dataset: TestDataset,
      files: Vec<Vec<AlignmentRecord>>,
      dense: bool,
      partition_rates: PartitionRates,
    ) -> Result<Timetree, Report> {
      let parse = nwk_read(dataset.tree().as_bytes())?;
      let names = parse.names();
      let sources = files
        .into_iter()
        .map(|records| sequence_source(Alphabet::default(), leaf_seq_inputs(&parse.graph, &names, records)))
        .collect::<Result<Vec<_>, Report>>()?;
      let dates = dates_by_node(dataset.dates()?, &parse.graph, &names);
      let params = TimetreeParams {
        dense: Some(dense),
        partition_rates,
        keep_root: false,
        max_iter: 1,
        ..marginal_timetree_params()
      };
      let input = TimetreeInput {
        graph: parse.graph,
        names,
        sources,
        dates: Some(dates),
        branch_lengths: parse.branch_lengths,
      };
      let output = run(&params, input, None, None, &NoopCancel, &NoopProgress, &NoopProgress)
        .map_err(OperationError::into_report)?;
      let edges = output
        .graph
        .get_edges()
        .map(|edge| {
          (
            output.names[&edge.source()].clone(),
            output.names[&edge.target()].clone(),
          )
        })
        .collect();
      let node_dates = output
        .node_dates
        .into_iter()
        .map(|(key, date)| (key, date.expect("every node of a timetree has a date")))
        .collect();
      let log_lh = output
        .results
        .iter()
        .map(|result| result.log_lh.map(LogLh::value))
        .sum::<Option<f64>>()
        .expect("a marginal timetree reports the sequence log-likelihood");
      Ok(Timetree {
        edges,
        node_dates,
        total_length: output.results.iter().map(|result| result.length).sum(),
        numbers: output.results.iter().map(|result| result.number).collect(),
        log_lh,
      })
    }
  }
}
