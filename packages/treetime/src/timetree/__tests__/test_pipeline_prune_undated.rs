#[cfg(test)]
mod tests {
  use crate::cancel::NoopCancel;
  use crate::clock::prune_tips::PruneTips;
  use crate::error::OperationError;
  use crate::gtr::get_gtr::{GtrModel, JC69Params};
  use crate::progress::{NoopProgress, RunWarning, RunWarningKind};
  use crate::test_utils::{RecordingLog, dates_by_node, find_node_key_by_name, marginal_timetree_params};
  use crate::timetree::params::{BranchLengthMode, TimetreeParams};
  use crate::timetree::pipeline::{self, SequenceLength, TimetreeInput, TimetreeSequences};
  use eyre::Report;
  use itertools::Itertools;
  use maplit::btreeset;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::collections::BTreeSet;
  use treetime_io::nwk::nwk_read;
  use treetime_primitives::date::DateConstraint;
  use treetime_utils::{assert_error, o, vec_of_owned};

  const TREE: &str = "((A:0.002,B:0.003):0.001,((C:0.002,D:0.004):0.002,(E:0.003,F:0.001):0.001):0.001);";

  const SEQUENCE_LENGTH: usize = 1000;

  const CLOCK_RATE: f64 = 1e-3;

  const CLOCK_FILTER_IQD: f64 = 4.0;

  #[rustfmt::skip]
  #[rstest]
  #[case::default(               (false, false, false), (vec_of_owned!["A", "B", "C", "D", "E", "F"], vec![]))]
  #[case::undated(               (false, true,  false), (vec_of_owned!["A", "B", "C", "D", "E"],      vec_of_owned!["F"]))]
  #[case::undated_keep_listed(   (false, true,  true),  (vec_of_owned!["A", "B", "C", "D", "E", "F"], vec![]))]
  #[case::filtered(              (true,  false, false), (vec_of_owned!["A", "B", "C", "D", "E", "F"], vec![]))]
  #[case::filtered_and_undated(  (true,  true,  false), (vec_of_owned!["A", "B", "C", "D", "E"],      vec_of_owned!["F"]))]
  #[trace]
  fn test_pipeline_prune_undated_tips_removes_the_tips_without_a_date(
    #[case] (filtered, undated, keep_f): (bool, bool, bool),
    #[case] (expected_leaves, expected_removed): (Vec<String>, Vec<String>),
  ) -> Result<(), Report> {
    let input = helpers::input()?;
    let keep = if keep_f {
      btreeset! {find_node_key_by_name(&input.graph, &input.names, "F").expect("F is a leaf")}
    } else {
      BTreeSet::new()
    };
    let params = TimetreeParams {
      prune_tips: PruneTips { filtered, undated, keep },
      ..helpers::params()
    };
    let log = RecordingLog::default();

    let output = pipeline::run(&params, input, None, None, &NoopCancel, &NoopProgress, &log)
      .map_err(OperationError::into_report)?;

    let expected_warnings = (!expected_removed.is_empty())
      .then(|| RunWarning {
        kind: RunWarningKind::UndatedTipsRemoved,
        message: format!("Removed 1 tip without a usable date from the tree ({})", expected_removed.join(", ")),
        names: expected_removed,
      })
      .into_iter()
      .collect_vec();
    let leaves = output
      .graph
      .get_leaves()
      .map(|leaf| output.names[leaf.key()].to_owned())
      .sorted()
      .collect_vec();
    let warnings = log
      .run_warnings()
      .into_iter()
      .filter(|warning| {
        matches!(
          warning.kind,
          RunWarningKind::UndatedTipsRemoved | RunWarningKind::FilteredTipsRemoved
        )
      })
      .collect_vec();
    let nodes_with_dates_in_graph = output
      .dates
      .as_ref()
      .expect("the run has dates")
      .keys()
      .all(|key| output.graph.get_node(key).is_some());
    assert_eq!(
      (expected_leaves, expected_warnings, BTreeSet::new(), true),
      (leaves, warnings, output.outliers, nodes_with_dates_in_graph)
    );
    Ok(())
  }

  #[test]
  fn test_pipeline_prune_undated_tips_requires_dates() -> Result<(), Report> {
    let input = TimetreeInput {
      dates: None,
      ..helpers::input()?
    };
    let params = TimetreeParams {
      prune_tips: PruneTips {
        undated: true,
        ..PruneTips::default()
      },
      ..helpers::params()
    };

    assert_error!(
      pipeline::run(&params, input, None, None, &NoopCancel, &NoopProgress, &NoopProgress)
        .map_err(OperationError::into_report),
      "--prune-undated-tips removes the tips without a usable date, but the run has no dates"
    );
    Ok(())
  }

  mod helpers {
    use super::*;

    pub(super) fn params() -> TimetreeParams {
      TimetreeParams {
        branch_length_mode: BranchLengthMode::Input,
        clock_rate: Some(CLOCK_RATE),
        clock_filter: CLOCK_FILTER_IQD,
        seed: 1,
        ..marginal_timetree_params()
      }
    }

    pub(super) fn input() -> Result<TimetreeInput, Report> {
      let nwk_parsed = nwk_read(TREE.as_bytes())?;
      let dates = [
        ("A", Some(2020.0)),
        ("B", Some(2021.0)),
        ("C", Some(2020.42)),
        ("D", Some(2022.0)),
        ("E", Some(2021.42)),
        ("F", None),
      ]
      .into_iter()
      .map(|(name, date)| (o!(name), date.map(DateConstraint::exact)));
      let names = nwk_parsed.names.clone();
      Ok(TimetreeInput {
        dates: Some(dates_by_node(dates, &nwk_parsed.graph, &names)),
        graph: nwk_parsed.graph,
        names,
        sequences: TimetreeSequences::Length(SequenceLength {
          length: Some(SEQUENCE_LENGTH),
          model: GtrModel::JC69(JC69Params::default()),
          alphabet: None,
        }),
        branch_lengths: nwk_parsed.branch_lengths,
      })
    }
  }
}
