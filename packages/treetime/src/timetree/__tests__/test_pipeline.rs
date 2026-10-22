#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::Alphabet;
  use crate::cancel::NoopCancel;
  use crate::clock::divergence::root_to_node_divergences;
  use crate::error::OperationError;
  use crate::optimize::params::BranchLengthMode;
  use crate::partition::result::{PartitionResult, PartitionSequences};
  use crate::progress::{NoopProgress, RunWarning, RunWarningKind, StageSink};
  use crate::seq::sink::SeqSink;
  use crate::test_utils::{
    RecordingLog, RecordingSeqSink, TestDataset, dates_by_node, find_node_key_by_name, leaf_seq_inputs,
    marginal_timetree_params, sequence_source,
  };
  use crate::timetree::params::{TimeMarginalMode, TimetreeParams};
  use crate::timetree::pipeline::{self, TimetreeInput, TimetreeOutput};
  use eyre::Report;
  use itertools::Itertools;
  use maplit::btreemap;
  use parking_lot::Mutex;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::collections::{BTreeMap, BTreeSet};
  use treetime_graph::graph::Graph;
  use treetime_graph::node::GraphNodeKey;
  use treetime_graph::node_names::NodeNames;
  use treetime_io::dates_csv::DateConstraint;
  use treetime_io::nwk::nwk_read;
  use treetime_primitives::date::{DateRange, DateValue};
  use treetime_utils::{assert_error, o, pretty_assert_abs_diff_eq, pretty_assert_ulps_eq};

  const STAR_TREE: &str = "(A:0.5,B:1.0,C:1.5,D:2.0)root;";

  const STAR_CLOCK_RATE: f64 = 0.25;

  const CLOCK_FILTER_IQD: f64 = 3.0;

  const STEM_NAME: &str = "OUTER";

  const FORCED_OUTLIER_DATE: f64 = 2000.0;

  const ZIKA_CLOCK_RATE: f64 = 1e-3;

  const ZIKA_SEQUENCE_LENGTH: usize = 10_807;

  const ZIKA_COALESCENT_TC: f64 = 0.01;

  const TIP_PAST_OPTIMUM_TREE: &str = "(A:0.001,(B:0.002,C:0.003)N:0.001)R;";

  const TIP_PAST_OPTIMUM_RATE: f64 = 0.01;

  const TIP_PAST_OPTIMUM_SEQUENCE_LENGTH: usize = 1000;

  #[test]
  fn test_pipeline_marginal_mode_requires_an_alignment() -> Result<(), Report> {
    let params = TimetreeParams {
      clock_rate: Some(STAR_CLOCK_RATE),
      ..marginal_timetree_params()
    };
    let input = helpers::star_input(Some(helpers::star_dates()))?;

    assert_error!(
      helpers::run(&params, input, None),
      "Alignment required for marginal reconstruction"
    );
    Ok(())
  }

  #[test]
  fn test_pipeline_rejects_a_sequence_sink_in_input_mode() -> Result<(), Report> {
    let params = TimetreeParams {
      branch_length_mode: BranchLengthMode::Input,
      sequence_outputs_requested: true,
      ..marginal_timetree_params()
    };
    let input = helpers::star_input(Some(helpers::star_dates()))?;

    assert_error!(
      helpers::run(&params, input, Some(&mut RecordingSeqSink::default())),
      "Reconstructed sequence output requires ancestral reconstruction; incompatible with --branch-length-mode=input"
    );
    Ok(())
  }

  #[test]
  fn test_pipeline_rejects_a_sequence_sink_without_requested_sequence_outputs() -> Result<(), Report> {
    let params = TimetreeParams {
      sequence_outputs_requested: false,
      ..marginal_timetree_params()
    };
    let input = helpers::star_input(Some(helpers::star_dates()))?;

    assert_error!(
      helpers::run(&params, input, Some(&mut RecordingSeqSink::default())),
      "A sequence sink was passed, but the parameters request no sequence outputs"
    );
    Ok(())
  }

  #[test]
  fn test_pipeline_clock_filter_outlier_is_a_bad_branch_of_the_time_inference() -> Result<(), Report> {
    let mut dates = helpers::zika_dates()?;
    let outlier_name = dates.keys().next().expect("the metadata has samples").clone();
    dates.insert(outlier_name.clone(), Some(DateConstraint::exact(FORCED_OUTLIER_DATE)));
    let input = helpers::zika_input(&helpers::zika_newick()?, dates, true)?;
    let params = TimetreeParams {
      clock_filter: CLOCK_FILTER_IQD,
      clock_rate: Some(ZIKA_CLOCK_RATE),
      ..marginal_timetree_params()
    };

    let output = helpers::run(&params, input, None)?;

    let outlier = find_node_key_by_name(&output.graph, &output.names, &outlier_name).expect("outlier leaf must exist");
    assert!(
      output.outliers.contains(&outlier),
      "the leaf dated {FORCED_OUTLIER_DATE} must be a clock outlier"
    );
    assert!(
      output.bad_branches[&outlier],
      "a clock outlier must be a bad branch of the time inference"
    );
    let flagged_good: BTreeSet<GraphNodeKey> = output
      .outliers
      .iter()
      .filter(|key| !output.bad_branches[*key])
      .copied()
      .collect();
    assert_eq!(BTreeSet::new(), flagged_good);
    Ok(())
  }

  #[test]
  fn test_pipeline_warns_about_a_clock_outlier_dated_outside_its_input_date_range() -> Result<(), Report> {
    let mut dates = helpers::zika_dates()?;
    let outlier_name = dates.keys().next().expect("the metadata has samples").clone();
    dates.insert(
      outlier_name.clone(),
      Some(DateConstraint {
        raw: o!("forced outlier"),
        value: DateValue::Range(DateRange {
          start: FORCED_OUTLIER_DATE - 1.0,
          end: FORCED_OUTLIER_DATE,
        }),
      }),
    );
    let input = helpers::zika_input(&helpers::zika_newick()?, dates, true)?;
    let params = TimetreeParams {
      clock_filter: CLOCK_FILTER_IQD,
      clock_rate: Some(ZIKA_CLOCK_RATE),
      ..marginal_timetree_params()
    };
    let log = RecordingLog::default();

    let output = pipeline::run(&params, input, None, None, &NoopCancel, &NoopProgress, &log)
      .map_err(OperationError::into_report)?;

    let outlier = find_node_key_by_name(&output.graph, &output.names, &outlier_name).expect("outlier leaf must exist");
    let excess = output.node_dates[&outlier].expect("the outlier is dated") - FORCED_OUTLIER_DATE;
    let expected = vec![RunWarning {
      kind: RunWarningKind::DatesOutsideInputDates,
      message: format!(
        "Some final dates lie outside the input dates: 1 tip dated outside the input date range of the tip, by up to \
         {excess:.2} years ({outlier_name}). The clock alone dates these nodes, as it does for the tips that the clock \
         filter flagged and for their ancestors"
      ),
      names: vec![outlier_name],
    }];
    let actual = log
      .run_warnings()
      .into_iter()
      .filter(|warning| warning.kind == RunWarningKind::DatesOutsideInputDates)
      .collect_vec();
    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  fn test_pipeline_does_not_warn_about_dates_without_clock_outliers() -> Result<(), Report> {
    let input = helpers::zika_input(&helpers::zika_newick()?, helpers::zika_dates()?, true)?;
    let params = TimetreeParams {
      clock_rate: Some(ZIKA_CLOCK_RATE),
      ..marginal_timetree_params()
    };
    let log = RecordingLog::default();

    pipeline::run(&params, input, None, None, &NoopCancel, &NoopProgress, &log).map_err(OperationError::into_report)?;

    let kinds = log.run_warnings().into_iter().map(|warning| warning.kind).collect_vec();
    assert!(
      !kinds.contains(&RunWarningKind::DatesOutsideInputDates),
      "no date may lie outside the input dates, got warnings {kinds:?}"
    );
    Ok(())
  }

  #[test]
  fn test_pipeline_input_mode_undated_single_child_root_fails_with_point_division_known_issue_h_timetree_input_branch_lengths_abort_on_point_division()
  -> Result<(), Report> {
    let input = helpers::zika_input(&helpers::zika_newick_with_stem()?, helpers::zika_dates()?, false)?;
    let params = TimetreeParams {
      branch_length_mode: BranchLengthMode::Input,
      keep_root: false,
      clock_filter: CLOCK_FILTER_IQD,
      sequence_length: Some(ZIKA_SEQUENCE_LENGTH),
      ..marginal_timetree_params()
    };

    let Err(report) = helpers::run(&params, input, None) else {
      panic!("input mode must fail on this tree until the point division is defined");
    };

    assert_eq!(
      "Cannot divide point by point: operation not well-defined",
      report.root_cause().to_string()
    );
    Ok(())
  }

  #[test]
  fn test_pipeline_reroot_past_a_tip_dates_the_tip_at_its_input_date() -> Result<(), Report> {
    let params = TimetreeParams {
      branch_length_mode: BranchLengthMode::Input,
      keep_root: false,
      clock_rate: Some(TIP_PAST_OPTIMUM_RATE),
      sequence_length: Some(TIP_PAST_OPTIMUM_SEQUENCE_LENGTH),
      ..marginal_timetree_params()
    };

    let output = helpers::run(&params, helpers::tip_past_optimum_input()?, None)?;

    let a_key = find_node_key_by_name(&output.graph, &output.names, "A").expect("A exists");
    assert!(
      output.graph.is_leaf(a_key),
      "v0 roots a split at the tip end with Biopython root_with_outgroup(tip): a new root with a 0-length branch to \
       the tip (packages/legacy/treetime/treetime/treeregression.py#L443-L463)"
    );
    pretty_assert_ulps_eq!(
      2000.0,
      output.node_dates[&a_key].expect("the tip is dated"),
      max_ulps = 0
    );
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::regression_rate(None)]
  #[case::marginal_likelihood_rate(Some(ZIKA_COALESCENT_TC))]
  #[trace]
  fn test_pipeline_estimated_rate_dates_are_inferred_at_the_reported_rate(
    #[case] coalescent: Option<f64>,
  ) -> Result<(), Report> {
    let params = |time_marginal| TimetreeParams {
      coalescent,
      time_marginal,
      ..marginal_timetree_params()
    };

    let output = helpers::run(&params(TimeMarginalMode::Never), helpers::zika_input(&helpers::zika_newick()?, helpers::zika_dates()?, true)?, None)?;
    let at_reported_rate = helpers::run(
      &params(TimeMarginalMode::OnlyFinal),
      helpers::zika_input(&helpers::zika_newick()?, helpers::zika_dates()?, true)?,
      None,
    )?;

    pretty_assert_ulps_eq!(at_reported_rate.clock_model.clock_rate(), output.clock_model.clock_rate(), max_ulps = 0);
    for (key, date) in &output.node_dates {
      pretty_assert_abs_diff_eq!(
        at_reported_rate.node_dates[key].expect("every node is dated"),
        date.expect("every node is dated"),
        epsilon = 1e-9
      );
    }
    Ok(())
  }

  #[test]
  fn test_pipeline_refinement_progress_counts_iterations_from_one() -> Result<(), Report> {
    let params = TimetreeParams {
      clock_rate: Some(ZIKA_CLOCK_RATE),
      max_iter: 2,
      ..marginal_timetree_params()
    };
    let input = helpers::zika_input(&helpers::zika_newick()?, helpers::zika_dates()?, true)?;
    let stages = helpers::RecordingStages::default();

    pipeline::run(&params, input, None, None, &NoopCancel, &stages, &NoopProgress)
      .map_err(OperationError::into_report)?;

    let iterations: Vec<(String, String)> = stages
      .reports
      .into_inner()
      .into_iter()
      .filter(|(stage, _, message)| stage == "Refinement loop" && !message.is_empty())
      .map(|(_, fraction, message)| (format!("{fraction:.2}"), message))
      .collect();
    let expected = vec![(o!("0.30"), o!("iteration 1/2")), (o!("0.55"), o!("iteration 2/2"))];
    assert_eq!(expected, iterations);
    Ok(())
  }

  #[test]
  fn test_pipeline_date_branch_lengths_sum_to_node_dates() -> Result<(), Report> {
    let input = helpers::zika_input(&helpers::zika_newick()?, helpers::zika_dates()?, true)?;
    let params = TimetreeParams {
      clock_rate: Some(ZIKA_CLOCK_RATE),
      ..marginal_timetree_params()
    };

    let output = helpers::run(&params, input, None)?;

    let root_date = output.node_dates[&output.graph.root_key()?].expect("the root is dated");
    let summed = root_to_node_divergences(&output.graph, |edge_key| {
      output.date_branch_lengths[&edge_key].expect("every branch of a dated tree has a length")
    })?;
    for (key, date) in &output.node_dates {
      pretty_assert_abs_diff_eq!(
        date.expect("every node is dated"),
        root_date + summed[key],
        epsilon = 1e-9
      );
    }
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::never(                    (TimeMarginalMode::Never,     false, None),       false)]
  #[case::always(                   (TimeMarginalMode::Always,    false, None),       true)]
  #[case::only_final(               (TimeMarginalMode::OnlyFinal, false, None),       true)]
  #[case::confidence_with_rate_std( (TimeMarginalMode::Never,     true,  Some(1e-4)), true)]
  #[case::always_with_rate_std(     (TimeMarginalMode::Always,    true,  Some(1e-4)), true)]
  #[trace]
  fn test_pipeline_gives_a_confidence_interval_for_every_dated_node_when_marginal_times_are_requested(
    #[case] (time_marginal, confidence, clock_std_dev): (TimeMarginalMode, bool, Option<f64>),
    #[case] expected_intervals: bool,
  ) -> Result<(), Report> {
    let input = helpers::zika_input(&helpers::zika_newick()?, helpers::zika_dates()?, true)?;
    let params = TimetreeParams {
      clock_rate: Some(ZIKA_CLOCK_RATE),
      time_marginal,
      confidence,
      clock_std_dev,
      ..marginal_timetree_params()
    };

    let output = helpers::run(&params, input, None)?;

    let dated_nodes = output.node_dates.values().flatten().count();
    let interval_count = output.confidence_intervals.map(|intervals| intervals.len());
    assert_eq!(expected_intervals.then_some(dated_nodes), interval_count);
    Ok(())
  }

  #[test]
  fn test_pipeline_date_branch_lengths_need_both_end_dates() -> Result<(), Report> {
    let parsed = nwk_read(b"((A:1,B:1)X:1,C:1)root;".as_slice())?;
    let names = parsed.names.clone();
    let graph = parsed.graph;
    let key = |name: &str| find_node_key_by_name(&graph, &names, name).expect("node must exist");
    let node_dates = btreemap! {
      key("root") => Some(2000.0),
      key("X") => Some(2001.5),
      key("A") => Some(2004.0),
      key("B") => None,
      key("C") => Some(2003.25),
    };

    let actual: BTreeMap<String, Option<f64>> = pipeline::date_branch_lengths(&graph, &node_dates)
      .into_iter()
      .map(|(edge_key, length)| {
        let child = graph.get_edge(edge_key).expect("edge must exist").target();
        (names[child].to_owned(), length)
      })
      .collect();

    let expected = btreemap! {
      o!("X") => Some(1.5),
      o!("A") => Some(2.5),
      o!("B") => None,
      o!("C") => Some(3.25),
    };
    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  fn test_pipeline_sequence_outputs_do_not_depend_on_the_sequence_sink() -> Result<(), Report> {
    let params = TimetreeParams {
      clock_rate: Some(ZIKA_CLOCK_RATE),
      sequence_outputs_requested: true,
      ..marginal_timetree_params()
    };
    let input = helpers::zika_input(&helpers::zika_newick()?, helpers::zika_dates()?, true)?;
    let without_sink = helpers::run(&params, input, None)?;
    let input = helpers::zika_input(&helpers::zika_newick()?, helpers::zika_dates()?, true)?;
    let with_sink = helpers::run(&params, input, Some(&mut RecordingSeqSink::default()))?;

    let edge_count = without_sink.graph.get_edges().count();
    let without_sink = helpers::only_sequences(without_sink);
    let with_sink = helpers::only_sequences(with_sink);
    assert_eq!(edge_count, without_sink.edge_mutations.len());
    assert_eq!(without_sink.root_sequence, with_sink.root_sequence);
    assert_eq!(without_sink.edge_mutations, with_sink.edge_mutations);
    assert_eq!(without_sink.edge_mutation_counts, with_sink.edge_mutation_counts);
    Ok(())
  }

  #[test]
  fn test_pipeline_marginal_mode_removes_an_undated_single_child_root() -> Result<(), Report> {
    let newick = helpers::zika_newick_with_stem()?;
    let samples = nwk_read((helpers::zika_newick()?).as_bytes())?;
    let sample_names = helpers::leaf_names(&samples.graph, &samples.names);
    let input = helpers::zika_input(&newick, helpers::zika_dates()?, true)?;
    let params = TimetreeParams {
      keep_root: false,
      clock_filter: CLOCK_FILTER_IQD,
      ..marginal_timetree_params()
    };

    let output = helpers::run(&params, input, None)?;

    let output_leaves = helpers::leaf_names(&output.graph, &output.names);
    assert!(!output_leaves.contains(STEM_NAME));
    assert_eq!(sample_names, output_leaves);
    Ok(())
  }

  mod helpers {
    use super::*;

    pub(super) fn run(
      params: &TimetreeParams,
      input: TimetreeInput,
      seq_sink: Option<&mut dyn SeqSink>,
    ) -> Result<TimetreeOutput, Report> {
      pipeline::run(params, input, None, seq_sink, &NoopCancel, &NoopProgress, &NoopProgress)
        .map_err(OperationError::into_report)
    }

    #[derive(Default)]
    pub(super) struct RecordingStages {
      pub reports: Mutex<Vec<(String, f64, String)>>,
    }

    impl StageSink for RecordingStages {
      fn report(&self, stage: &str, fraction: f64, message: &str) {
        self
          .reports
          .lock()
          .push((stage.to_owned(), fraction, message.to_owned()));
      }
    }

    pub(super) fn star_dates() -> BTreeMap<String, Option<DateConstraint>> {
      btreemap! {
        "A".to_owned() => Some(DateConstraint::exact(2002.0)),
        "B".to_owned() => Some(DateConstraint::exact(2004.0)),
        "C".to_owned() => Some(DateConstraint::exact(2006.0)),
        "D".to_owned() => Some(DateConstraint::exact(2008.0)),
      }
    }

    pub(super) fn only_sequences(output: TimetreeOutput) -> PartitionSequences {
      let [result] = <[PartitionResult; 1]>::try_from(output.results).expect("one alignment gives one result");
      result.sequences.expect("sequence outputs were requested")
    }

    pub(super) fn star_input(dates: Option<BTreeMap<String, Option<DateConstraint>>>) -> Result<TimetreeInput, Report> {
      let nwk_parsed = nwk_read(STAR_TREE.as_bytes())?;
      let names = nwk_parsed.names.clone();
      let dates = dates.map(|dates| dates_by_node(dates, &nwk_parsed.graph, &names));
      let input = TimetreeInput {
        graph: nwk_parsed.graph,
        names,
        sources: vec![],
        dates,
        branch_lengths: nwk_parsed.branch_lengths,
      };
      Ok(input)
    }

    pub(super) fn tip_past_optimum_input() -> Result<TimetreeInput, Report> {
      let nwk_parsed = nwk_read(TIP_PAST_OPTIMUM_TREE.as_bytes())?;
      let names = nwk_parsed.names.clone();
      let dates = btreemap! {
        o!("A") => Some(DateConstraint::exact(2000.0)),
        o!("B") => Some(DateConstraint::exact(2001.0)),
        o!("C") => Some(DateConstraint::exact(2001.0)),
      };
      let dates = dates_by_node(dates, &nwk_parsed.graph, &names);
      Ok(TimetreeInput {
        graph: nwk_parsed.graph,
        names,
        sources: vec![],
        dates: Some(dates),
        branch_lengths: nwk_parsed.branch_lengths,
      })
    }

    pub(super) fn zika_newick() -> Result<String, Report> {
      Ok(TestDataset::Zika20.tree().trim().to_owned())
    }

    pub(super) fn zika_newick_with_stem() -> Result<String, Report> {
      let tree = zika_newick()?;
      let body = tree.strip_suffix(';').expect("a Newick tree ends with a semicolon");
      Ok(format!("({body}:0.001){STEM_NAME};"))
    }

    pub(super) fn zika_dates() -> Result<BTreeMap<String, Option<DateConstraint>>, Report> {
      TestDataset::Zika20.dates()
    }

    pub(super) fn zika_input(
      newick: &str,
      dates: BTreeMap<String, Option<DateConstraint>>,
      with_alignment: bool,
    ) -> Result<TimetreeInput, Report> {
      let alphabet = Alphabet::default();
      let nwk_parsed = nwk_read(newick.as_bytes())?;
      let names = nwk_parsed.names.clone();
      let source = if with_alignment {
        let records = TestDataset::Zika20.alignment()?;
        Some(sequence_source(
          alphabet,
          leaf_seq_inputs(&nwk_parsed.graph, &names, records),
        )?)
      } else {
        None
      };
      let dates = dates_by_node(dates, &nwk_parsed.graph, &names);
      let input = TimetreeInput {
        graph: nwk_parsed.graph,
        names,
        sources: source.into_iter().collect(),
        dates: Some(dates),
        branch_lengths: nwk_parsed.branch_lengths,
      };
      Ok(input)
    }

    pub(super) fn leaf_names(graph: &Graph, names: &NodeNames) -> BTreeSet<String> {
      graph.get_leaves().map(|leaf| names[leaf.key()].to_owned()).collect()
    }
  }
}
