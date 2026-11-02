#[cfg(test)]
mod tests {
  use crate::cancel::NoopCancel;
  use crate::clock::pipeline::{self, ClockInput, ClockOutput, ClockParams};
  use crate::clock::prune_tips::PruneTips;
  use crate::progress::{NoopProgress, RunWarning, RunWarningKind};
  use crate::reroot::params::RerootSpec;
  use crate::test_utils::{RecordingLog, dates_by_node, find_node_key_by_name};
  use eyre::Report;
  use itertools::Itertools;
  use maplit::btreeset;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::collections::BTreeSet;
  use treetime_primitives::date::DateConstraint;
  use treetime_utils::{assert_error, o, vec_of_owned};

  const CLOCK_FILTER_IQD: f64 = 3.0;

  const CLOCK_TREE: &str = "((A:0.104,B:0.106)AB:0.001,(C:0.126,D:0.124)CD:0.001,(E:0.144,F:0.146)EF:0.001,\
                            (G:0.166,H:0.164)GH:0.001,(I:0.184,J:0.186)IJ:0.001,(K:0.206,L:0.204)KL:0.001,\
                            (X:0.08,Y:0.11)XY:0.001)R;";

  #[rstest]
  #[case::reroot(Some(RerootSpec::default()))]
  #[case::keep_root(None)]
  #[trace]
  fn test_pipeline_prune_filtered_tips_removes_the_flagged_tips_that_the_keep_list_does_not_name(
    #[case] reroot_spec: Option<RerootSpec>,
  ) -> Result<(), Report> {
    let input = helpers::input(None)?;
    let kept = find_node_key_by_name(&input.graph, &input.names, "Y").expect("Y is a leaf");
    let params = ClockParams {
      clock_filter: CLOCK_FILTER_IQD,
      allow_negative_rate: false,
      reroot_spec,
      prune_tips: PruneTips {
        filtered: true,
        undated: false,
        keep: btreeset! {kept},
      },
    };
    let log = RecordingLog::default();

    let output = helpers::run(&params, input, &log)?;

    let expected = (
      vec_of_owned!["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "Y"],
      vec_of_owned!["Y"],
      vec![RunWarning {
        kind: RunWarningKind::FilteredTipsRemoved,
        message: o!("Removed 1 tip that the clock filter flagged from the tree (X)"),
        names: vec_of_owned!["X"],
      }],
      vec_of_owned!["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "Y"],
    );
    let rtt_tips = output
      .regression_results
      .iter()
      .filter(|result| {
        output
          .graph
          .get_leaves()
          .any(|leaf| output.names[leaf.key()] == result.name)
      })
      .map(|result| result.name.clone())
      .sorted()
      .collect_vec();
    assert_eq!(
      expected,
      (
        helpers::leaf_names(&output),
        helpers::names_of(&output.outliers, &output),
        log.run_warnings(),
        rtt_tips,
      )
    );
    Ok(())
  }

  #[test]
  fn test_pipeline_prune_filtered_tips_rejects_removing_the_outgroup() -> Result<(), Report> {
    let input = helpers::input(None)?;
    let outgroup_key = find_node_key_by_name(&input.graph, &input.names, "Y").expect("Y is a leaf");
    let params = ClockParams {
      clock_filter: CLOCK_FILTER_IQD,
      allow_negative_rate: false,
      reroot_spec: Some(RerootSpec::Tips(vec![outgroup_key])),
      prune_tips: PruneTips {
        filtered: true,
        ..PruneTips::default()
      },
    };

    assert_error!(
      helpers::run(&params, input, &RecordingLog::default()),
      "The clock filter flagged the outgroup tip Y, and --prune-filtered-tips would remove it. Name it in --keep-tips \
       to keep it in the tree"
    );
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::undated(              (0.0,              false, None),      (vec_of_owned!["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "X", "Y"],      vec![helpers::undated_warning()]))]
  #[case::undated_keep_listed(  (0.0,              false, Some("L")), (vec_of_owned!["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "X", "Y"], vec![]))]
  #[case::filtered_and_undated( (CLOCK_FILTER_IQD, true,  Some("Y")), (vec_of_owned!["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "Y"],           vec![helpers::filtered_warning(), helpers::undated_warning()]))]
  #[trace]
  fn test_pipeline_prune_undated_tips_removes_the_tips_without_a_date(
    #[case] (clock_filter, filtered, kept_name): (f64, bool, Option<&str>),
    #[case] (expected_leaves, expected_warnings): (Vec<String>, Vec<RunWarning>),
  ) -> Result<(), Report> {
    let input = helpers::input(Some("L"))?;
    let keep = kept_name
      .map(|name| find_node_key_by_name(&input.graph, &input.names, name).expect("the kept tip is a leaf"))
      .into_iter()
      .collect();
    let params = ClockParams {
      clock_filter,
      allow_negative_rate: false,
      reroot_spec: Some(RerootSpec::default()),
      prune_tips: PruneTips {
        filtered,
        undated: true,
        keep,
      },
    };
    let log = RecordingLog::default();

    let output = helpers::run(&params, input, &log)?;

    assert_eq!(
      (expected_leaves, expected_warnings),
      (helpers::leaf_names(&output), log.run_warnings())
    );
    Ok(())
  }

  #[test]
  fn test_pipeline_prune_undated_tips_rejects_removing_the_outgroup() -> Result<(), Report> {
    let input = helpers::input(Some("L"))?;
    let outgroup_key = find_node_key_by_name(&input.graph, &input.names, "L").expect("L is a leaf");
    let params = ClockParams {
      clock_filter: 0.0,
      allow_negative_rate: false,
      reroot_spec: Some(RerootSpec::Tips(vec![outgroup_key])),
      prune_tips: PruneTips {
        undated: true,
        ..PruneTips::default()
      },
    };

    assert_error!(
      helpers::run(&params, input, &RecordingLog::default()),
      "The outgroup tip L has no usable date, and --prune-undated-tips would remove it. Name it in --keep-tips to keep \
       it in the tree"
    );
    Ok(())
  }

  mod helpers {
    use super::*;
    use crate::error::OperationError;
    use crate::progress::LogSink;
    use treetime_graph::node::GraphNodeKey;
    use treetime_io::nwk::nwk_read;

    pub(super) fn run(params: &ClockParams, input: ClockInput, log: &dyn LogSink) -> Result<ClockOutput, Report> {
      pipeline::run(params, input, &NoopCancel, &NoopProgress, log).map_err(OperationError::into_report)
    }

    pub(super) fn input(undated: Option<&str>) -> Result<ClockInput, Report> {
      let nwk_parsed = nwk_read(CLOCK_TREE.as_bytes())?;
      let dates = [
        ("A", 2000.0),
        ("B", 2001.0),
        ("C", 2002.0),
        ("D", 2003.0),
        ("E", 2004.0),
        ("F", 2005.0),
        ("G", 2006.0),
        ("H", 2007.0),
        ("I", 2008.0),
        ("J", 2009.0),
        ("K", 2010.0),
        ("L", 2011.0),
        ("X", 2004.0),
        ("Y", 2007.0),
      ]
      .into_iter()
      .filter(|&(name, _)| Some(name) != undated)
      .map(|(name, date)| (name.to_owned(), Some(DateConstraint::exact(date))));
      let names = nwk_parsed.names.clone();
      Ok(ClockInput {
        dates: dates_by_node(dates, &nwk_parsed.graph, &names),
        graph: nwk_parsed.graph,
        names,
        branch_lengths: nwk_parsed.branch_lengths,
      })
    }

    pub(super) fn filtered_warning() -> RunWarning {
      RunWarning {
        kind: RunWarningKind::FilteredTipsRemoved,
        message: o!("Removed 1 tip that the clock filter flagged from the tree (X)"),
        names: vec_of_owned!["X"],
      }
    }

    pub(super) fn undated_warning() -> RunWarning {
      RunWarning {
        kind: RunWarningKind::UndatedTipsRemoved,
        message: o!("Removed 1 tip without a usable date from the tree (L)"),
        names: vec_of_owned!["L"],
      }
    }

    pub(super) fn leaf_names(output: &ClockOutput) -> Vec<String> {
      output
        .graph
        .get_leaves()
        .map(|leaf| output.names[leaf.key()].to_owned())
        .sorted()
        .collect()
    }

    pub(super) fn names_of(keys: &BTreeSet<GraphNodeKey>, output: &ClockOutput) -> Vec<String> {
      keys.iter().map(|&key| output.names[key].to_owned()).sorted().collect()
    }
  }
}
