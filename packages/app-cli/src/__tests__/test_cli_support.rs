#[cfg(test)]
pub(crate) mod tests {
  pub(crate) mod helpers {
    use crate::cli::treetime_cli::TreetimeArgs;
    use clap::{Arg, Command, CommandFactory};
    use itertools::Itertools;

    pub(crate) const TT: &str = "timetree";
    pub(crate) const AN: &str = "ancestral";
    pub(crate) const CL: &str = "clock";
    pub(crate) const MU: &str = "mugration";
    pub(crate) const HO: &str = "homoplasy";
    pub(crate) const OP: &str = "optimize";
    pub(crate) const PR: &str = "prune";
    pub(crate) const PI: &str = "pipeline";
    pub(crate) const EXAMPLES_GET: &str = "examples get";
    pub(crate) const SCHEMA_WRITE: &str = "schema write";
    pub(crate) const COMPLETIONS: &str = "completions";
    pub(crate) const HELP_MARKDOWN: &str = "help-markdown";
    pub(crate) const ROOT: &str = "treetime";

    pub(crate) const ANALYSIS: &[&str] = &[TT, AN, CL, MU, HO, OP, PR];
    pub(crate) const RUN: &[&str] = &[TT, AN, CL, MU, HO, OP, PR, PI];
    pub(crate) const EVERY: &[&str] = &[
      TT,
      AN,
      CL,
      MU,
      HO,
      OP,
      PR,
      PI,
      EXAMPLES_GET,
      SCHEMA_WRITE,
      COMPLETIONS,
      HELP_MARKDOWN,
    ];

    #[derive(Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
    pub(crate) enum Pending {
      SettingsThatApplyOnlyWithAnother,
      ResolvePolytomiesByDefault,
      ConfidenceAsTheOnlyIntervalSwitch,
      RepresentationAndRootingGroup,
      NewickAnnotationsFlag,
      MugrationConfidenceFlag,
      FlagRenames,
      ShortHelp,
      RootCommand,
      DryRun,
    }

    pub(crate) fn treetime_command() -> Command {
      let mut command = TreetimeArgs::command();
      command.build();
      command
    }

    pub(crate) fn all_commands(root: &Command) -> Vec<(String, &Command)> {
      let mut commands = vec![(ROOT.to_owned(), root)];
      collect_subcommands_recursive(root, None, &mut commands);
      commands
    }

    pub(crate) fn leaf_commands(root: &Command) -> Vec<(String, &Command)> {
      all_commands(root)
        .into_iter()
        .filter(|(path, command)| path != ROOT && !command.has_subcommands() && !command.is_hide_set())
        .collect_vec()
    }

    pub(crate) fn find_command<'a>(root: &'a Command, path: &str) -> Option<&'a Command> {
      path
        .split(' ')
        .try_fold(root, |command, name| command.find_subcommand(name))
    }

    pub(crate) fn long_flags(command: &Command) -> Vec<String> {
      command
        .get_arguments()
        .filter_map(Arg::get_long)
        .map(|long| format!("--{long}"))
        .sorted()
        .collect_vec()
    }

    pub(crate) fn is_clap_flag(arg: &Arg) -> bool {
      matches!(arg.get_id().as_str(), "help" | "version")
    }

    fn collect_subcommands_recursive<'a>(
      command: &'a Command,
      prefix: Option<&str>,
      commands: &mut Vec<(String, &'a Command)>,
    ) {
      for sub in command.get_subcommands().filter(|sub| sub.get_name() != "help") {
        let path = prefix.map_or_else(
          || sub.get_name().to_owned(),
          |prefix| format!("{prefix} {}", sub.get_name()),
        );
        commands.push((path.clone(), sub));
        collect_subcommands_recursive(sub, Some(&path), commands);
      }
    }
  }
}
