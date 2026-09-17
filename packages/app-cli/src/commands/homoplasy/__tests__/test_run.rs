#[cfg(test)]
mod tests {
  use crate::commands::ancestral::args::TreetimeAncestralArgsRaw;
  use crate::commands::homoplasy::args::{TreetimeHomoplasyArgs, TreetimeHomoplasyArgsRaw};
  use crate::commands::homoplasy::run::run_homoplasy;
  use pretty_assertions::assert_eq;
  use std::path::PathBuf;

  #[test]
  fn test_run_homoplasy_reports_not_implemented() {
    let raw = TreetimeHomoplasyArgsRaw {
      ancestral_args: TreetimeAncestralArgsRaw {
        tree: Some(PathBuf::from("tree.nwk")),
        ..TreetimeAncestralArgsRaw::default()
      },
      ..TreetimeHomoplasyArgsRaw::default()
    };
    let args = TreetimeHomoplasyArgs::try_from(raw).unwrap();

    let err = run_homoplasy(args).unwrap_err();
    assert_eq!("The homoplasy operation is not yet implemented in v1", err.to_string());
  }
}
