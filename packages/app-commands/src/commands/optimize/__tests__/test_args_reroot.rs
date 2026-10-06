#[cfg(test)]
mod tests {
  use crate::commands::optimize::args::{OptimizeRerootMethod, TreetimeOptimizeArgs, TreetimeOptimizeArgsRaw};
  use helpers::args_with;
  use pretty_assertions::assert_eq;
  use treetime::clock::find_best_root::params::{RerootMethod, RerootSpec};
  use treetime::o;
  use treetime_io::nwk::nwk_read;

  const TREE: &str = "(A:0.1,B:0.1,C:0.1)root;";

  #[test]
  fn test_optimize_args_reroot_spec_default_keeps_root() {
    let tree = nwk_read(TREE.as_bytes()).unwrap();
    let args = args_with(None, vec![], false);
    assert_eq!(None, args.reroot_spec(&tree.graph, &tree.names()).unwrap());
  }

  #[test]
  fn test_optimize_args_reroot_spec_keep_root_flag_keeps_root() {
    let tree = nwk_read(TREE.as_bytes()).unwrap();
    let args = args_with(None, vec![], true);
    assert_eq!(None, args.reroot_spec(&tree.graph, &tree.names()).unwrap());
  }

  #[test]
  fn test_optimize_args_reroot_spec_keep_root_overrides_method() {
    let tree = nwk_read(TREE.as_bytes()).unwrap();
    let args = args_with(Some(OptimizeRerootMethod::MinDev), vec![], true);
    assert_eq!(None, args.reroot_spec(&tree.graph, &tree.names()).unwrap());
  }

  #[test]
  fn test_optimize_args_reroot_spec_min_dev() {
    let tree = nwk_read(TREE.as_bytes()).unwrap();
    let args = args_with(Some(OptimizeRerootMethod::MinDev), vec![], false);
    assert_eq!(
      Some(RerootSpec::Method(RerootMethod::MinDev)),
      args.reroot_spec(&tree.graph, &tree.names()).unwrap()
    );
  }

  #[test]
  fn test_optimize_args_reroot_spec_tips() {
    let tree = nwk_read(TREE.as_bytes()).unwrap();
    let args = args_with(None, vec![o!("A"), o!("B")], false);
    assert_eq!(
      Some(RerootSpec::Tips(
        tree.graph.get_leaves().map(|leaf| leaf.key()).take(2).collect()
      )),
      args.reroot_spec(&tree.graph, &tree.names()).unwrap()
    );
  }

  #[test]
  fn test_optimize_reroot_method_converts_to_reroot_method() {
    let method: RerootMethod = OptimizeRerootMethod::MinDev.into();
    assert_eq!(RerootMethod::MinDev, method);
  }

  mod helpers {
    use super::*;

    pub(super) fn args_with(
      reroot: Option<OptimizeRerootMethod>,
      reroot_tips: Vec<String>,
      keep_root: bool,
    ) -> TreetimeOptimizeArgs {
      TreetimeOptimizeArgs::try_from(TreetimeOptimizeArgsRaw {
        tree: Some("tree.nwk".into()),
        reroot,
        reroot_tips,
        keep_root,
        ..Default::default()
      })
      .unwrap()
    }
  }
}
