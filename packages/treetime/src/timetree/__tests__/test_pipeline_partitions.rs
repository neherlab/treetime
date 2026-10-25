#[cfg(test)]
mod tests {
  use crate::optimize::partition_rates::PartitionRates;
  use crate::test_utils::{Columns, TestDataset, split_at_gap_free_column};
  use eyre::Report;
  use helpers::{Split, codon_case, timetree};
  use ndarray::Array1;
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

  #[rustfmt::skip]
  #[rstest]
  #[case::zika_dense(  TestDataset::Zika20,  true,  Columns::Leading(198))]
  #[case::flu_dense(   TestDataset::Flu20,   true,  Columns::Leading(198))]
  #[case::ebola_dense( TestDataset::Ebola20, true,  Columns::Leading(198))]
  #[case::zika_sparse( TestDataset::Zika20,  false, Columns::WholeCodons  )]
  #[case::flu_sparse(  TestDataset::Flu20,   false, Columns::WholeCodons  )]
  #[case::ebola_sparse(TestDataset::Ebola20, false, Columns::WholeCodons  )]
  #[trace]
  fn test_timetree_codon_positions_give_the_dates_and_mutations_of_one_alignment(
    #[case] dataset: TestDataset,
    #[case] dense: bool,
    #[case] selection: Columns,
  ) -> Result<(), Report> {
    let records = dataset.alignment_columns(selection)?;
    let whole = codon_case(dataset, records.clone(), dense, Split::Whole)?;
    let split = codon_case(dataset, records, dense, Split::CodonPositions)?;

    assert_eq!(
      (whole.edges, whole.total_length, vec![1], whole.sequences),
      (split.edges, split.total_length, split.numbers, split.sequences)
    );
    pretty_assert_abs_diff_eq!(
      Array1::from_iter(whole.indel_rates.iter().flatten().copied()),
      Array1::from_iter(split.indel_rates.iter().flatten().copied()),
      epsilon = 1e-8
    );
    pretty_assert_map_abs_diff_eq!(whole.node_dates, split.node_dates, epsilon = 1e-8);
    pretty_assert_abs_diff_eq!(whole.log_lh, split.log_lh, epsilon = 1e-8);
    Ok(())
  }

  mod helpers {
    use crate::alphabet::alphabet::Alphabet;
    use crate::cancel::NoopCancel;
    use crate::error::OperationError;
    use crate::gtr::get_gtr::{GtrModel, JC69Params};
    use crate::optimize::partition_rates::PartitionRates;
    use crate::progress::NoopProgress;
    use crate::reroot::params::RerootSpec;
    use crate::seq::alignment::SequenceSource;
    use crate::seq::mutation::Mutation;
    use crate::test_utils::{
      TestDataset, codon_positions_input, dates_by_node, jc69_timetree_sequences, leaf_seq_inputs,
      marginal_timetree_params, sequence_source,
    };
    use crate::timetree::params::TimetreeParams;
    use crate::timetree::pipeline::{TimetreeInput, TimetreeSequences, run};
    use eyre::Report;
    use std::collections::BTreeSet;
    use treetime_graph::key_map::{EdgeMap, NodeMap};

    use treetime_io::nwk::nwk_read;
    use treetime_primitives::{AlignmentRecord, LogLh, Seq};

    pub(super) struct Timetree {
      pub(super) edges: BTreeSet<(String, String)>,
      pub(super) node_dates: NodeMap<f64>,
      pub(super) total_length: usize,
      pub(super) numbers: Vec<usize>,
      pub(super) sequences: Vec<Option<(Seq, EdgeMap<Vec<Mutation>>)>>,
      pub(super) indel_rates: Vec<Option<f64>>,
      pub(super) log_lh: f64,
    }

    pub(super) fn timetree(
      dataset: TestDataset,
      files: Vec<Vec<AlignmentRecord>>,
      dense: bool,
      partition_rates: PartitionRates,
    ) -> Result<Timetree, Report> {
      let options = RunOptions {
        dense,
        partition_rates,
        resolve_polytomies: false,
        sequence_outputs_requested: false,
      };
      run_timetree(dataset, options, jc69_timetree_sequences, files)
    }

    #[derive(Clone, Copy)]
    pub(super) enum Split {
      Whole,
      CodonPositions,
    }

    pub(super) fn codon_case(
      dataset: TestDataset,
      records: Vec<AlignmentRecord>,
      dense: bool,
      split: Split,
    ) -> Result<Timetree, Report> {
      let options = RunOptions {
        dense,
        partition_rates: PartitionRates::Equal,
        resolve_polytomies: true,
        sequence_outputs_requested: true,
      };
      run_timetree(
        dataset,
        options,
        |sources| match split {
          Split::Whole => jc69_timetree_sequences(sources),
          Split::CodonPositions => Ok(TimetreeSequences::Alignments(codon_positions_input(
            sources,
            &GtrModel::JC69(JC69Params::default()),
          )?)),
        },
        vec![records],
      )
    }

    struct RunOptions {
      dense: bool,
      partition_rates: PartitionRates,
      resolve_polytomies: bool,
      sequence_outputs_requested: bool,
    }

    fn run_timetree(
      dataset: TestDataset,
      RunOptions {
        dense,
        partition_rates,
        resolve_polytomies,
        sequence_outputs_requested,
      }: RunOptions,
      sequences: impl FnOnce(Vec<SequenceSource>) -> Result<TimetreeSequences, Report>,
      files: Vec<Vec<AlignmentRecord>>,
    ) -> Result<Timetree, Report> {
      let parse = nwk_read(dataset.tree().as_bytes())?;
      let names = parse.names.clone();
      let sources = files
        .into_iter()
        .map(|records| sequence_source(Alphabet::default(), leaf_seq_inputs(&parse.graph, &names, records)))
        .collect::<Result<Vec<_>, Report>>()?;
      let dates = dates_by_node(dataset.dates()?, &parse.graph, &names);
      let params = TimetreeParams {
        dense: Some(dense),
        partition_rates,
        reroot_spec: Some(RerootSpec::default()),
        max_iter: 1,
        resolve_polytomies,
        sequence_outputs_requested,
        ..marginal_timetree_params()
      };
      let input = TimetreeInput {
        graph: parse.graph,
        names,
        sequences: sequences(sources)?,
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
            output.names[edge.source()].to_owned(),
            output.names[edge.target()].to_owned(),
          )
        })
        .collect();
      let node_dates = output
        .node_dates
        .into_iter()
        .map(|(key, date)| (key, date.expect("every node of a timetree has a date")))
        .collect();
      let log_lh = output
        .partitions
        .iter()
        .map(|result| result.log_lh.map(LogLh::value))
        .sum::<Option<f64>>()
        .expect("a marginal timetree reports the sequence log-likelihood");
      let sequences = output
        .alignments
        .iter()
        .map(|alignment| {
          alignment
            .sequences
            .as_ref()
            .map(|sequences| (sequences.root_sequence.clone(), sequences.edge_mutations.clone()))
        })
        .collect();
      Ok(Timetree {
        edges,
        node_dates,
        total_length: output.alignments.iter().map(|result| result.length).sum(),
        numbers: output.alignments.iter().map(|result| result.number).collect(),
        sequences,
        indel_rates: output.alignments.iter().map(|result| result.indel_rate).collect(),
        log_lh,
      })
    }
  }
}
