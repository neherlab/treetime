#[cfg(test)]
mod tests {
  use crate::clock::find_best_root::params::{RerootMethod, RerootSpec};
  use crate::commands::optimize::args::TreetimeOptimizeArgs;
  use crate::commands::shared::reroot::RerootArgs;
  use crate::o;
  use pretty_assertions::assert_eq;
  use rstest::rstest;

  fn args_with(reroot: RerootArgs, keep_root: bool) -> TreetimeOptimizeArgs {
    TreetimeOptimizeArgs {
      reroot,
      keep_root,
      ..TreetimeOptimizeArgs::default()
    }
  }

  #[test]
  fn test_optimize_args_reroot_spec_default_keeps_root() {
    let args = TreetimeOptimizeArgs::default();
    assert_eq!(args.reroot_spec().unwrap(), None);
  }

  #[test]
  fn test_optimize_args_reroot_spec_keep_root_flag_keeps_root() {
    let args = args_with(RerootArgs::default(), true);
    assert_eq!(args.reroot_spec().unwrap(), None);
  }

  #[test]
  fn test_optimize_args_reroot_spec_min_dev() {
    let args = args_with(
      RerootArgs {
        reroot: Some(RerootMethod::MinDev),
        ..RerootArgs::default()
      },
      false,
    );
    assert_eq!(
      args.reroot_spec().unwrap(),
      Some(RerootSpec::Method(RerootMethod::MinDev))
    );
  }

  #[test]
  fn test_optimize_args_reroot_spec_tips() {
    let args = args_with(
      RerootArgs {
        reroot_tips: vec![o!("A"), o!("B")],
        ..RerootArgs::default()
      },
      false,
    );
    assert_eq!(
      args.reroot_spec().unwrap(),
      Some(RerootSpec::Tips(vec![o!("A"), o!("B")]))
    );
  }

  // Date-dependent methods cannot run in optimize (no sampling dates) and must be
  // rejected up front rather than producing a NaN clock objective.
  #[rstest]
  #[case::least_squares(RerootMethod::LeastSquares)]
  #[case::oldest(RerootMethod::Oldest)]
  #[case::clock_filter(RerootMethod::ClockFilter)]
  #[trace]
  fn test_optimize_args_reroot_spec_rejects_date_methods(#[case] method: RerootMethod) {
    let args = args_with(
      RerootArgs {
        reroot: Some(method),
        ..RerootArgs::default()
      },
      false,
    );
    let err = args
      .reroot_spec()
      .expect_err("date-dependent method should be rejected");
    assert!(err.to_string().contains("date-free"), "unexpected error: {err}");
  }
}
