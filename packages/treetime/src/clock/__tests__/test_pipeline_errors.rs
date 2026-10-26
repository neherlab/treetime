#[cfg(test)]
mod tests {
  use crate::cancel::NoopCancel;
  use crate::clock::pipeline::{self, ClockInput, ClockOutput, ClockParams};
  use crate::error::OperationError;
  use crate::progress::NoopProgress;
  use crate::reroot::params::RerootSpec;
  use crate::test_utils::dates_by_node;
  use eyre::Report;
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::collections::BTreeMap;
  use treetime_io::nwk::nwk_read;
  use treetime_primitives::date::DateConstraint;
  use treetime_utils::error::report_to_string;
  use treetime_utils::o;

  #[rustfmt::skip]
  #[rstest]
  #[case::same_dates(     btreemap! { o!("A") => 2000.0, o!("B") => 2000.0, o!("C") => 2000.0 }, "No variation in sampling dates! Please specify your clock rate explicitly.")]
  #[case::no_dates(       btreemap! {},                                                             "No valid date information found: no node of the tree has a usable date in the dates input"   )]
  #[trace]
  fn test_pipeline_errors_classify_data_problems_as_invalid_input(
    #[case] dates: BTreeMap<String, f64>,
    #[case] expected: &str,
  ) -> Result<(), Report> {
    let parsed = nwk_read(b"((A:0.1,B:0.2)X:0.1,C:0.3)root;".as_slice())?;
    let names = parsed.names.clone();
    let dates = dates_by_node(
      dates.into_iter().map(|(name, date)| (name, Some(DateConstraint::exact(date)))),
      &parsed.graph,
      &names,
    );
    let params = ClockParams {
      clock_filter: 3.0,
      allow_negative_rate: false,
      reroot_spec: None,
      prune_filtered_tips: None,
    };
    let input = ClockInput {
      graph: parsed.graph,
      names,
      dates,
      branch_lengths: parsed.branch_lengths,
    };

    let result = pipeline::run(&params, input, &NoopCancel, &NoopProgress, &NoopProgress);

    let Err(OperationError::InvalidInput(report)) = result else {
      panic!("a data problem must be classified as invalid input");
    };
    assert_eq!(expected, report_to_string(&report));
    Ok(())
  }

  #[test]
  fn test_pipeline_errors_prefilter_requires_a_positive_rate_also_when_negative_rates_are_allowed() -> Result<(), Report>
  {
    let result = helpers::run_negative_rate_tree(3.0);

    let Err(OperationError::InvalidInput(report)) = result else {
      panic!("a tree without a positive-rate root must fail the clock fit of the outlier filter as invalid input");
    };
    assert_eq!(
      "The clock fit of the outlier filter requires a positive clock rate, also with --allow-negative-rate. Run with \
       --clock-filter=0 to skip the outlier filter: Clock rate is negative for all root positions. The data may lack \
       temporal signal. Please specify --clock-rate explicitly.",
      report_to_string(&report)
    );
    Ok(())
  }

  #[test]
  fn test_pipeline_errors_negative_rate_is_allowed_without_the_outlier_filter() -> Result<(), Report> {
    let output = helpers::run_negative_rate_tree(0.0).map_err(OperationError::into_report)?;

    assert!(
      output.clock_model.clock_rate() < 0.0,
      "the tree has no positive-rate root, got rate {:.6e}",
      output.clock_model.clock_rate()
    );
    Ok(())
  }

  mod helpers {
    use super::*;

    const NEGATIVE_RATE_TREE: &str = "((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root;";

    pub(super) fn run_negative_rate_tree(clock_filter: f64) -> Result<ClockOutput, OperationError> {
      let parsed = nwk_read(NEGATIVE_RATE_TREE.as_bytes()).map_err(OperationError::InvalidInput)?;
      let names = parsed.names.clone();
      let dates = btreemap! {
        o!("A") => Some(DateConstraint::exact(2017.0)),
        o!("B") => Some(DateConstraint::exact(2005.0)),
        o!("C") => Some(DateConstraint::exact(2010.0)),
        o!("D") => Some(DateConstraint::exact(2022.0)),
      };
      let params = ClockParams {
        clock_filter,
        allow_negative_rate: true,
        reroot_spec: Some(RerootSpec::default()),
        prune_filtered_tips: None,
      };
      let input = ClockInput {
        dates: dates_by_node(dates, &parsed.graph, &names),
        names: names.clone(),
        graph: parsed.graph,
        branch_lengths: parsed.branch_lengths,
      };
      pipeline::run(&params, input, &NoopCancel, &NoopProgress, &NoopProgress)
    }
  }
}
