#[cfg(test)]
mod tests {
  use crate::clock::find_best_root::params::{RerootMethod, RerootSpec};
  use crate::commands::optimize::args::{OptimizeRerootMethod, TreetimeOptimizeArgs};
  use crate::o;
  use pretty_assertions::assert_eq;

  fn args_with(
    reroot: Option<OptimizeRerootMethod>,
    reroot_tips: Vec<String>,
    keep_root: bool,
  ) -> TreetimeOptimizeArgs {
    TreetimeOptimizeArgs {
      reroot,
      reroot_tips,
      keep_root,
      ..TreetimeOptimizeArgs::default()
    }
  }

  #[test]
  fn test_optimize_args_reroot_spec_default_keeps_root() {
    let args = TreetimeOptimizeArgs::default();
    assert_eq!(None, args.reroot_spec());
  }

  #[test]
  fn test_optimize_args_reroot_spec_keep_root_flag_keeps_root() {
    let args = args_with(None, vec![], true);
    assert_eq!(None, args.reroot_spec());
  }

  #[test]
  fn test_optimize_args_reroot_spec_keep_root_overrides_method() {
    let args = args_with(Some(OptimizeRerootMethod::MinDev), vec![], true);
    assert_eq!(None, args.reroot_spec());
  }

  #[test]
  fn test_optimize_args_reroot_spec_min_dev() {
    let args = args_with(Some(OptimizeRerootMethod::MinDev), vec![], false);
    assert_eq!(Some(RerootSpec::Method(RerootMethod::MinDev)), args.reroot_spec());
  }

  #[test]
  fn test_optimize_args_reroot_spec_tips() {
    let args = args_with(None, vec![o!("A"), o!("B")], false);
    assert_eq!(Some(RerootSpec::Tips(vec![o!("A"), o!("B")])), args.reroot_spec());
  }

  #[test]
  fn test_optimize_reroot_method_converts_to_reroot_method() {
    let method: RerootMethod = OptimizeRerootMethod::MinDev.into();
    assert_eq!(RerootMethod::MinDev, method);
  }
}
