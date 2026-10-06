#[cfg(test)]
mod tests {
  use crate::cancel::NoopCancel;
  use crate::clock::clock_regression::ClockVarianceParams;
  use crate::clock::find_best_root::params::{BranchPointOptimizationParams, RerootSpec};
  use crate::clock::pipeline::{self, ClockInput, ClockParams};
  use crate::error::OperationError;
  use crate::progress::NoopProgress;
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
    let names = parsed.names();
    let dates = dates_by_node(
      dates.into_iter().map(|(name, date)| (name, Some(DateConstraint::exact(date)))),
      &parsed.graph,
      &names,
    );
    let params = ClockParams {
      clock_params: ClockVarianceParams::default(),
      clock_filter: 3.0,
      keep_root: true,
      allow_negative_rate: false,
      branch_params: BranchPointOptimizationParams::default(),
      reroot_spec: RerootSpec::default(),
    };
    let input = ClockInput {
      graph: parsed.graph,
      dates,
      branch_lengths: parsed.branch_lengths,
    };

    let result = pipeline::run(&params, input, &names, &NoopCancel, &NoopProgress, &NoopProgress);

    let Err(OperationError::InvalidInput(report)) = result else {
      panic!("a data problem must be classified as invalid input");
    };
    assert_eq!(expected, report_to_string(&report));
    Ok(())
  }
}
