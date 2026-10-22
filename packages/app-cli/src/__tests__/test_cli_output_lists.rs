#[cfg(test)]
mod tests {
  use crate::__tests__::test_cli_flag_families::tests::helpers::per_file_output_flags;
  use crate::__tests__::test_cli_support::tests::helpers::{
    AN, CL, HO, MU, OP, PR, TT, find_command, treetime_command,
  };
  use app_backend::results::outputs::result_outputs;
  use app_commands::command::AppCommand;
  use app_output::command_outputs::CommandKind;
  use app_output::output_conditions::OutputSettings;
  use app_output::output_plan::OutputSelection;
  use itertools::Itertools;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::collections::BTreeSet;

  #[rustfmt::skip]
  #[rstest]
  #[case::timetree( TT, CommandKind::Timetree)]
  #[case::ancestral(AN, CommandKind::Ancestral)]
  #[case::clock(    CL, CommandKind::Clock)]
  #[case::mugration(MU, CommandKind::Mugration)]
  #[case::homoplasy(HO, CommandKind::Homoplasy)]
  #[case::optimize( OP, CommandKind::Optimize)]
  #[case::prune(    PR, CommandKind::Prune)]
  #[trace]
  fn test_cli_output_selection_lists_the_outputs_of_the_command(#[case] name: &str, #[case] kind: CommandKind) {
    let root = treetime_command();
    let command = find_command(&root, name).expect("analysis command exists");
    let selection = command
      .get_arguments()
      .find(|arg| arg.get_long() == Some("output-selection"))
      .expect("analysis command has --output-selection");

    let actual: BTreeSet<String> = selection
      .get_possible_values()
      .iter()
      .map(|value| value.get_name().to_owned())
      .collect();

    let expected: BTreeSet<String> = kind
      .outputs()
      .iter()
      .map(|(selection, _)| selection.as_ref().to_owned())
      .chain([OutputSelection::All.as_ref().to_owned()])
      .collect();
    assert_eq!(expected, actual);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::timetree( TT, CommandKind::Timetree)]
  #[case::ancestral(AN, CommandKind::Ancestral)]
  #[case::clock(    CL, CommandKind::Clock)]
  #[case::mugration(MU, CommandKind::Mugration)]
  #[case::homoplasy(HO, CommandKind::Homoplasy)]
  #[case::optimize( OP, CommandKind::Optimize)]
  #[case::prune(    PR, CommandKind::Prune)]
  #[trace]
  fn test_cli_per_file_output_flags_match_the_outputs_of_the_command(#[case] name: &str, #[case] kind: CommandKind) {
    let root = treetime_command();
    let command = find_command(&root, name).expect("analysis command exists");

    let expected: BTreeSet<String> = kind
      .outputs()
      .iter()
      .map(|(selection, _)| selection.flag_name().to_owned())
      .collect();
    assert_eq!(expected, per_file_output_flags(command));
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::timetree( CommandKind::Timetree,  &["auspice", "clock-model", "dates", "model", "nexus", "node-data", "nwk", "sequences"])]
  #[case::ancestral(CommandKind::Ancestral, &["model", "nexus", "node-data", "nwk", "sequences"])]
  #[case::clock(    CommandKind::Clock,     &["clock-model", "nexus", "nwk", "rtt", "rtt-chart"])]
  #[case::mugration(CommandKind::Mugration, &["model", "nexus", "node-data", "nwk", "traits"])]
  #[case::homoplasy(CommandKind::Homoplasy, &["homoplasy-report", "homoplasy-stats", "nexus", "nwk"])]
  #[case::optimize( CommandKind::Optimize,  &["model", "nexus", "node-data", "nwk"])]
  #[case::prune(    CommandKind::Prune,     &["nexus", "nwk"])]
  #[trace]
  fn test_cli_default_outputs_of_each_command(#[case] kind: CommandKind, #[case] expected: &[&str]) {
    let actual = kind
      .default_outputs(OutputSettings::default())
      .into_iter()
      .map(|selection| selection.as_ref().to_owned())
      .sorted()
      .collect_vec();
    assert_eq!(expected, actual);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::timetree( AppCommand::Timetree,  CommandKind::Timetree)]
  #[case::ancestral(AppCommand::Ancestral, CommandKind::Ancestral)]
  #[case::clock(    AppCommand::Clock,     CommandKind::Clock)]
  #[case::mugration(AppCommand::Mugration, CommandKind::Mugration)]
  #[case::homoplasy(AppCommand::Homoplasy, CommandKind::Homoplasy)]
  #[case::optimize( AppCommand::Optimize,  CommandKind::Optimize)]
  #[case::prune(    AppCommand::Prune,     CommandKind::Prune)]
  #[trace]
  fn test_cli_app_reads_only_outputs_the_command_writes(#[case] command: AppCommand, #[case] kind: CommandKind) {
    let writable: BTreeSet<OutputSelection> = kind.outputs().iter().map(|&(selection, _)| selection).collect();
    let unwritable = result_outputs(command)
      .iter()
      .filter(|selection| !writable.contains(selection))
      .collect_vec();
    assert_eq!(Vec::<&OutputSelection>::new(), unwritable);
  }
}
