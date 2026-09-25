#[cfg(test)]
mod tests {
  use crate::command::AppCommand;
  use crate::config::source::{ConfigProblem, InvalidConfig};
  use indoc::indoc;
  use pretty_assertions::assert_eq;
  use rstest::rstest;

  #[rustfmt::skip]
  #[rstest]
  #[case::two_coalescent_priors( AppCommand::Timetree, "coalescent: 1.0\ncoalescent_opt: true\n",       "config::conflict",    "`coalescent` cannot be used together with `coalescent_opt`")]
  #[case::skyline_and_optimized( AppCommand::Timetree, "coalescent_opt: true\ncoalescent_skyline: true\n", "config::conflict",  "`coalescent_opt` cannot be used together with `coalescent_skyline`")]
  #[case::keep_root_and_reroot(  AppCommand::Clock,    "keep_root: true\nreroot: min-dev\n",             "config::conflict",    "`keep_root` cannot be used together with `reroot`")]
  #[case::optimize_reroot_tips(  AppCommand::Optimize, "reroot: min-dev\nreroot_tips: [A]\n",            "config::conflict",    "`reroot` cannot be used together with `reroot_tips`")]
  #[case::relax_with_one_value(  AppCommand::Timetree, "relax: [1.0]\n",                                 "config::value-count", "`relax` has 1 value, but `--relax` takes 2 values each time it is given")]
  #[case::relax_with_three_values(AppCommand::Timetree, "relax: [1.0, 0.0, 2.0]\n",                     "config::value-count", "`relax` has 3 values, but `--relax` takes 2 values each time it is given")]
  #[trace]
  fn test_cli_rules_reject_what_the_command_line_rejects(
    #[case] command: AppCommand,
    #[case] text: &str,
    #[case] code: &str,
    #[case] message: &str,
  ) {
    let report = command.prepare_text("config.yaml", text).err().unwrap();
    let invalid = report.downcast_ref::<InvalidConfig>().unwrap();
    assert_eq!(
      vec![(code, message)],
      invalid
        .problems
        .iter()
        .map(|problem: &ConfigProblem| (problem.code.as_str(), problem.message.as_str()))
        .collect::<Vec<_>>()
    );
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::one_coalescent_prior(     AppCommand::Timetree, "coalescent: 1.0\n")]
  #[case::keep_root_alone(          AppCommand::Timetree, "keep_root: true\n")]
  #[case::keep_root_with_unset_root(AppCommand::Clock,    "metadata: m.tsv\nkeep_root: true\nreroot: null\n")]
  #[case::relax_with_two_values(    AppCommand::Timetree, "relax: [1.0, 0.0]\n")]
  #[case::relax_twice(              AppCommand::Timetree, "relax: [1.0, 0.0, 2.0, 0.5]\n")]
  #[case::empty_relax(              AppCommand::Timetree, "relax: []\n")]
  #[case::list_setting(             AppCommand::Clock,    "metadata: m.tsv\nmetadata_id_columns: [strain, name, id]\n")]
  #[trace]
  fn test_cli_rules_accept_what_the_command_line_accepts(#[case] command: AppCommand, #[case] text: &str) {
    let error = command
      .prepare_text("config.yaml", text)
      .err()
      .map(|report| format!("{report:#}"));
    assert_eq!(None, error);
  }

  #[test]
  fn test_cli_rules_report_every_conflict_at_its_setting() {
    let text = indoc! {"
      coalescent: 1.0
      coalescent_opt: true
      coalescent_skyline: true
      keep_root: true
      reroot: min-dev
    "};
    let report = AppCommand::Timetree.prepare_text("config.yaml", text).err().unwrap();
    let invalid = report.downcast_ref::<InvalidConfig>().unwrap();
    assert_eq!(
      vec![
        ("`coalescent` cannot be used together with `coalescent_opt`", Some(0)),
        (
          "`coalescent` cannot be used together with `coalescent_skyline`",
          Some(0)
        ),
        (
          "`coalescent_opt` cannot be used together with `coalescent_skyline`",
          Some(16)
        ),
        ("`keep_root` cannot be used together with `reroot`", Some(16 + 21 + 25)),
      ],
      invalid
        .problems
        .iter()
        .map(|problem| (problem.message.as_str(), problem.span.map(|span| span.offset)))
        .collect::<Vec<_>>()
    );
  }
}
