#[cfg(test)]
mod tests {
  use crate::command::AppCommandExt;
  use app_commands::command::AppCommand;
  use app_commands::config::source::{ConfigProblem, InvalidConfig};
  use indoc::indoc;
  use pretty_assertions::assert_eq;
  use rstest::rstest;

  #[rustfmt::skip]
  #[rstest]
  #[case::keep_root_and_reroot(  AppCommand::Clock,    "keep_root: true\nreroot: min-dev\n",             "config::conflict",    "`keep_root` cannot be used together with `reroot`")]
  #[case::keep_root_and_oldest(  AppCommand::Clock,    "keep_root: true\nreroot: oldest\n",              "config::conflict",    "`keep_root` cannot be used together with `reroot`")]
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
  #[case::coalescent_number(        AppCommand::Timetree, "tree: t.nwk\nalignment: [a.fasta]\ncoalescent: 0.01\n")]
  #[case::coalescent_time(          AppCommand::Timetree, "tree: t.nwk\nalignment: [a.fasta]\ncoalescent: 10y\n")]
  #[case::coalescent_option(        AppCommand::Timetree, "tree: t.nwk\nalignment: [a.fasta]\ncoalescent: skyline\n")]
  #[case::keep_root_alone(          AppCommand::Timetree, "tree: t.nwk\nalignment: [a.fasta]\nkeep_root: true\n")]
  #[case::keep_root_with_default_root(AppCommand::Clock,  "tree: t.nwk\nmetadata: m.tsv\nkeep_root: true\n")]
  #[case::relax_with_two_values(    AppCommand::Timetree, "tree: t.nwk\nalignment: [a.fasta]\nrelax: [1.0, 0.0]\n")]
  #[case::relax_twice(              AppCommand::Timetree, "tree: t.nwk\nalignment: [a.fasta]\nrelax: [1.0, 0.0, 2.0, 0.5]\n")]
  #[case::empty_relax(              AppCommand::Timetree, "tree: t.nwk\nalignment: [a.fasta]\nrelax: []\n")]
  #[case::list_setting(             AppCommand::Clock,    "tree: t.nwk\nmetadata: m.tsv\nmetadata_id_columns: [strain, name, id]\n")]
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
      keep_root: true
      reroot: min-dev
      reroot_tips: [A]
    "};
    let report = AppCommand::Timetree.prepare_text("config.yaml", text).err().unwrap();
    let invalid = report.downcast_ref::<InvalidConfig>().unwrap();
    assert_eq!(
      vec![
        ("`keep_root` cannot be used together with `reroot`", Some(0)),
        ("`keep_root` cannot be used together with `reroot_tips`", Some(0)),
        ("`reroot` cannot be used together with `reroot_tips`", Some(16)),
      ],
      invalid
        .problems
        .iter()
        .map(|problem| (problem.message.as_str(), problem.span.map(|span| span.offset)))
        .collect::<Vec<_>>()
    );
  }
}
