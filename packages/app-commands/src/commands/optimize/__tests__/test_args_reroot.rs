#[cfg(test)]
mod tests {
  use crate::commands::optimize::args::{OptimizeRerootMethod, TreetimeOptimizeArgs, TreetimeOptimizeArgsRaw};
  use clap::Parser;
  use clap::error::ErrorKind;
  use helpers::args_with;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime::clock::find_best_root::params::{RerootMethod, RerootSpec};
  use treetime::o;
  use treetime_io::nwk::nwk_read;

  const TREE: &str = "(A:0.1,B:0.1,C:0.1)root;";

  #[test]
  fn test_optimize_args_reroot_spec_default_keeps_root() {
    let tree = nwk_read(TREE.as_bytes()).unwrap();
    let args = args_with(None, vec![], false);
    assert_eq!(None, helpers::reroot_spec(&args, &tree));
  }

  #[test]
  fn test_optimize_args_reroot_spec_keep_root_flag_keeps_root() {
    let tree = nwk_read(TREE.as_bytes()).unwrap();
    let args = args_with(None, vec![], true);
    assert_eq!(None, helpers::reroot_spec(&args, &tree));
  }

  #[test]
  fn test_optimize_args_reroot_spec_keep_root_overrides_method() {
    let tree = nwk_read(TREE.as_bytes()).unwrap();
    let args = args_with(Some(OptimizeRerootMethod::MinDev), vec![], true);
    assert_eq!(None, helpers::reroot_spec(&args, &tree));
  }

  #[test]
  fn test_optimize_args_reroot_spec_min_dev() {
    let tree = nwk_read(TREE.as_bytes()).unwrap();
    let args = args_with(Some(OptimizeRerootMethod::MinDev), vec![], false);
    assert_eq!(
      Some(RerootSpec::Method(RerootMethod::MinDev)),
      helpers::reroot_spec(&args, &tree)
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
      helpers::reroot_spec(&args, &tree)
    );
  }

  #[test]
  fn test_optimize_reroot_method_converts_to_reroot_method() {
    let method: RerootMethod = OptimizeRerootMethod::MinDev.into();
    assert_eq!(RerootMethod::MinDev, method);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::reroot_without_value(  &["--reroot"],                             ErrorKind::InvalidValue)]
  #[case::reroot_and_outgroup(   &["--reroot=min-dev", "--outgroup=A"],     ErrorKind::ArgumentConflict)]
  #[case::keep_root_and_outgroup(&["--keep-root", "--outgroup=A"],          ErrorKind::ArgumentConflict)]
  #[trace]
  fn test_optimize_args_rooting_rejects(#[case] flags: &[&str], #[case] expected: ErrorKind) {
    let args = ["optimize", "--tree=t.nwk", "--alignment=a.fasta"].into_iter().chain(flags.iter().copied());

    let result = TreetimeOptimizeArgsRaw::try_parse_from(args);

    assert_eq!(Some(expected), result.err().map(|err| err.kind()));
  }

  mod helpers {
    use super::*;
    use crate::commands::shared::alignment::AlignmentArgs;
    use crate::commands::shared::reroot::RootingArgs;
    use crate::commands::shared::tree_input::TreeArgs;
    use treetime_io::nwk::NwkParse;

    pub(super) fn reroot_spec(args: &TreetimeOptimizeArgs, tree: &NwkParse) -> Option<RerootSpec> {
      args
        .reroot
        .as_ref()
        .map(|reroot| reroot.spec(&tree.graph, &tree.names))
        .transpose()
        .unwrap()
    }

    pub(super) fn args_with(
      reroot: Option<OptimizeRerootMethod>,
      outgroup: Vec<String>,
      keep_root: bool,
    ) -> TreetimeOptimizeArgs {
      TreetimeOptimizeArgs::try_from(TreetimeOptimizeArgsRaw {
        tree: TreeArgs {
          tree: Some("tree.nwk".into()),
        },
        alignment: AlignmentArgs {
          alignment: vec!["a.fasta".into()],
        },
        reroot,
        rooting: RootingArgs { outgroup, keep_root },
        ..Default::default()
      })
      .unwrap()
    }
  }
}
