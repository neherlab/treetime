#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::Alphabet;
  use crate::cancel::NoopCancel;
  use crate::clock::divergence::root_to_node_divergences;
  use crate::error::OperationError;
  use crate::optimize::params::BranchLengthMode;
  use crate::progress::NoopProgress;
  use crate::seq::sink::SeqSink;
  use crate::test_utils::{RecordingSeqSink, find_node_key_by_name, marginal_timetree_params};
  use crate::timetree::params::TimetreeParams;
  use crate::timetree::pipeline::{self, TimetreeInput, TimetreeOutput};
  use eyre::Report;
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use std::collections::{BTreeMap, BTreeSet};
  use std::path::PathBuf;
  use treetime_graph::graph::Graph;
  use treetime_graph::node::GraphNodeKey;
  use treetime_io::csv::default_name_candidates;
  use treetime_io::dates_csv::{DateConstraint, DatesMap, read_dates};
  use treetime_io::fasta::read_many_fasta_path;
  use treetime_io::nwk::nwk_read_str;
  use treetime_primitives::AlignmentRecord;
  use treetime_utils::{assert_error, pretty_assert_abs_diff_eq};

  const STAR_TREE: &str = "(A:0.5,B:1.0,C:1.5,D:2.0)root;";

  const STAR_CLOCK_RATE: f64 = 0.25;

  const CLOCK_FILTER_IQD: f64 = 3.0;

  const STEM_NAME: &str = "OUTER";

  const FORCED_OUTLIER_DATE: f64 = 2000.0;

  const ZIKA_CLOCK_RATE: f64 = 1e-3;

  #[test]
  fn test_pipeline_rejects_n_branches_posterior_before_any_inference() -> Result<(), Report> {
    let params = TimetreeParams {
      n_branches_posterior: Some(2),
      ..marginal_timetree_params()
    };
    let input = helpers::star_input(None)?;

    assert_error!(
      helpers::run(&params, input, None),
      "--n-branches-posterior is not yet implemented"
    );
    Ok(())
  }

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
  fn test_pipeline_input_mode_undated_single_child_root_fails_with_point_division_known_issue_h_timetree_input_branch_lengths_abort_on_point_division()
  -> Result<(), Report> {
    let input = helpers::zika_input(&helpers::zika_newick_with_stem()?, helpers::zika_dates()?, false)?;
    let params = TimetreeParams {
      branch_length_mode: BranchLengthMode::Input,
      keep_root: false,
      clock_filter: CLOCK_FILTER_IQD,
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

    assert!(without_sink.root_sequence.is_some());
    assert_eq!(
      without_sink.graph.get_edges().count(),
      without_sink.edge_mutations.len()
    );
    assert_eq!(without_sink.root_sequence, with_sink.root_sequence);
    assert_eq!(without_sink.edge_mutations, with_sink.edge_mutations);
    Ok(())
  }

  #[test]
  fn test_pipeline_marginal_mode_removes_an_undated_single_child_root() -> Result<(), Report> {
    let newick = helpers::zika_newick_with_stem()?;
    let samples = nwk_read_str(&helpers::zika_newick()?)?;
    let sample_names = helpers::leaf_names(&samples.graph, &samples.names());
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

    pub(super) fn star_dates() -> DatesMap {
      btreemap! {
        "A".to_owned() => Some(DateConstraint::exact(2002.0)),
        "B".to_owned() => Some(DateConstraint::exact(2004.0)),
        "C".to_owned() => Some(DateConstraint::exact(2006.0)),
        "D".to_owned() => Some(DateConstraint::exact(2008.0)),
      }
    }

    pub(super) fn star_input(dates: Option<DatesMap>) -> Result<TimetreeInput, Report> {
      let nwk_parsed = nwk_read_str(STAR_TREE)?;
      let names = nwk_parsed.names();
      let input = TimetreeInput {
        graph: nwk_parsed.graph,
        names,
        alphabet: Alphabet::default(),
        sequences: None,
        dates,
        branch_lengths: nwk_parsed.branch_lengths,
      };
      Ok(input)
    }

    pub(super) fn zika_newick() -> Result<String, Report> {
      Ok(std::fs::read_to_string(zika_path("tree.nwk"))?.trim().to_owned())
    }

    pub(super) fn zika_newick_with_stem() -> Result<String, Report> {
      let tree = zika_newick()?;
      let body = tree.strip_suffix(';').expect("a Newick tree ends with a semicolon");
      Ok(format!("({body}:0.001){STEM_NAME};"))
    }

    pub(super) fn zika_dates() -> Result<DatesMap, Report> {
      read_dates(
        zika_path("metadata.tsv"),
        &['\t'],
        &default_name_candidates(),
        &None,
        &None,
      )
    }

    pub(super) fn zika_input(newick: &str, dates: DatesMap, with_alignment: bool) -> Result<TimetreeInput, Report> {
      let alphabet = Alphabet::default();
      let sequences = if with_alignment {
        let records = read_many_fasta_path(&[zika_path("aln.fasta.xz")], &alphabet)?;
        Some(records.into_iter().map(AlignmentRecord::from).collect())
      } else {
        None
      };
      let nwk_parsed = nwk_read_str(newick)?;
      let names = nwk_parsed.names();
      let input = TimetreeInput {
        graph: nwk_parsed.graph,
        names,
        alphabet,
        sequences,
        dates: Some(dates),
        branch_lengths: nwk_parsed.branch_lengths,
      };
      Ok(input)
    }

    pub(super) fn leaf_names(graph: &Graph, names: &BTreeMap<GraphNodeKey, Option<String>>) -> BTreeSet<String> {
      graph
        .get_leaves()
        .map(|leaf| names[&leaf.key()].clone().expect("every leaf is named"))
        .collect()
    }

    fn zika_path(file: &str) -> PathBuf {
      PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../data/zika/20")
        .join(file)
    }
  }
}
