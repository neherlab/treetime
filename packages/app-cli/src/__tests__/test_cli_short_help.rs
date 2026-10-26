#[cfg(test)]
mod tests {
  use crate::__tests__::test_cli_flag_families::tests::{FAMILIES, Family};
  use crate::__tests__::test_cli_support::tests::helpers::{
    AN, CL, HO, MU, OP, PR, ROOT, SCHEMA_WRITE, TT, all_commands, find_command, treetime_command,
  };
  use helpers::{Scope, ShortFlag, ShortFlagUse, ShortHelpDeviation, short_flag_deviations, short_help_deviation};
  use itertools::Itertools;
  use pretty_assertions::assert_eq;
  use std::collections::{BTreeMap, BTreeSet};

  const MAX_SHORT_HELP_LINES: usize = 30;

  #[rustfmt::skip]
  const SHORT_FLAGS: &[ShortFlag] = &[
    ShortFlag { short: 't', long: "tree",             scope: Scope::Every },
    ShortFlag { short: 'a', long: "alignment",        scope: Scope::Every },
    ShortFlag { short: 'm', long: "metadata",         scope: Scope::Every },
    ShortFlag { short: 'O', long: "output-all",       scope: Scope::Every },
    ShortFlag { short: 's', long: "output-selection", scope: Scope::Every },
    ShortFlag { short: 'j', long: "jobs",             scope: Scope::Every },
    ShortFlag { short: 'v', long: "verbose",          scope: Scope::Every },
    ShortFlag { short: 'q', long: "quiet",            scope: Scope::Every },
    ShortFlag { short: 'h', long: "help",             scope: Scope::Every },
    ShortFlag { short: 'V', long: "version",          scope: Scope::Only(&[ROOT]) },
    ShortFlag { short: 'n', long: "num-mut",          scope: Scope::Only(&[HO]) },
    ShortFlag { short: 'o', long: "output",           scope: Scope::Only(&[SCHEMA_WRITE]) },
  ];

  #[rustfmt::skip]
  const SHORT_HELP: &[(&str, &[&str])] = &[
    (TT, &["--tree", "--alignment", "--metadata", "--date-column", "--sequence-length", "--clock-rate", "--clock-std-dev",
           "--clock-filter", "--prune-filtered-tips", "--relax", "--reroot", "--outgroup", "--keep-root", "--branch-length-mode", "--confidence",
           "--keep-polytomies", "--max-iter", "--coalescent", "--n-skyline", "--model", "--output-all",
           "--output-selection", "--config", "--dry-run", "--seed", "--jobs", "--verbose", "--quiet", "--help"]),
    (AN, &["--tree", "--alignment", "--model", "--method-anc", "--include-leaves", "--impute-missing-data",
           "--reconstruct-tip-states", "--report-ambiguous", "--annotation", "--translations", "--cds-selection",
           "--aa-root-sequence", "--aa-model", "--output-all", "--output-selection", "--config", "--dry-run", "--seed",
           "--jobs", "--verbose", "--quiet", "--help"]),
    (CL, &["--tree", "--metadata", "--date-column", "--clock-filter", "--prune-filtered-tips", "--reroot", "--outgroup", "--keep-root",
           "--output-all", "--output-selection", "--config", "--dry-run", "--jobs", "--verbose", "--quiet", "--help"]),
    (MU, &["--tree", "--metadata", "--attribute", "--weights", "--missing-data", "--confidence", "--pc",
           "--sampling-bias-correction", "--output-all", "--output-selection", "--config", "--dry-run", "--jobs",
           "--verbose", "--quiet", "--help"]),
    (HO, &["--tree", "--alignment", "--model", "--method-anc", "--const", "--rescale", "--detailed", "--drms",
           "--num-mut", "--output-all", "--output-selection", "--config", "--dry-run", "--seed", "--jobs", "--verbose",
           "--quiet", "--help"]),
    (OP, &["--tree", "--alignment", "--model", "--max-iter", "--reroot", "--outgroup", "--keep-root", "--output-all",
           "--output-selection", "--config", "--dry-run", "--jobs", "--verbose", "--quiet", "--help"]),
    (PR, &["--tree", "--alignment", "--collapse-short", "--collapse-empty", "--merge-shared-mutations", "--remove",
           "--remove-file", "--output-all", "--output-selection", "--config", "--dry-run", "--jobs", "--verbose",
           "--quiet", "--help"]),
  ];

  const SHORT_HELP_ONLY_IN_ANCESTRAL: &[&str] = &[
    "--include-leaves",
    "--impute-missing-data",
    "--reconstruct-tip-states",
    "--report-ambiguous",
  ];

  #[test]
  fn test_cli_short_help_short_flags_come_from_one_list() {
    let root = treetime_command();
    let commands = all_commands(&root);
    let expected: BTreeSet<ShortFlagUse> = BTreeSet::new();

    let actual = short_flag_deviations(&commands, SHORT_FLAGS);

    assert_eq!(expected, actual);
  }

  #[test]
  fn test_cli_short_help_lists_the_main_flags() {
    let root = treetime_command();
    let expected: BTreeMap<&str, ShortHelpDeviation> = BTreeMap::new();

    let actual: BTreeMap<&str, ShortHelpDeviation> = SHORT_HELP
      .iter()
      .filter_map(|(path, listed)| {
        short_help_deviation(find_command(&root, path).unwrap(), listed, MAX_SHORT_HELP_LINES)
          .map(|deviation| (*path, deviation))
      })
      .collect();

    assert_eq!(expected, actual);
  }

  #[test]
  fn test_cli_short_help_shows_a_family_flag_in_every_command_of_the_family() {
    let expected: Vec<(&str, &str)> = vec![];

    let actual = FAMILIES
      .iter()
      .flat_map(|family: &Family| {
        let listed = family
          .commands
          .iter()
          .filter_map(|command| {
            SHORT_HELP
              .iter()
              .find(|(path, _)| path == command)
              .map(|(_, flags)| *flags)
          })
          .collect_vec();
        family.flags.iter().filter_map(move |flag| {
          let shown_in = listed.iter().filter(|flags| flags.contains(flag)).count();
          let partial = shown_in > 0 && shown_in < listed.len();
          (partial && !SHORT_HELP_ONLY_IN_ANCESTRAL.contains(flag)).then_some((family.name, *flag))
        })
      })
      .collect_vec();

    assert_eq!(expected, actual);
  }

  mod helpers {
    use clap::Command;
    use itertools::Itertools;
    use std::collections::BTreeSet;

    pub(super) struct ShortFlag {
      pub(super) short: char,
      pub(super) long: &'static str,
      pub(super) scope: Scope,
    }

    pub(super) enum Scope {
      Every,
      Only(&'static [&'static str]),
    }

    #[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
    pub(super) struct ShortFlagUse {
      pub(super) command: String,
      pub(super) long: String,
      pub(super) short: Option<char>,
    }

    #[derive(Clone, Debug, PartialEq, Eq)]
    pub(super) struct ShortHelpDeviation {
      pub(super) unlisted: BTreeSet<String>,
      pub(super) missing: BTreeSet<String>,
      pub(super) option_lines: usize,
    }

    pub(super) fn short_flag_deviations(
      commands: &[(String, &Command)],
      table: &[ShortFlag],
    ) -> BTreeSet<ShortFlagUse> {
      commands
        .iter()
        .flat_map(|(path, command)| {
          let in_scope = table
            .iter()
            .filter(|entry| match entry.scope {
              Scope::Every => true,
              Scope::Only(paths) => paths.contains(&path.as_str()),
            })
            .collect_vec();
          command.get_arguments().filter_map(move |arg| {
            let long = arg.get_long().unwrap_or_else(|| arg.get_id().as_str());
            let short = arg.get_short();
            let by_short = short.and_then(|short| in_scope.iter().find(|entry| entry.short == short));
            let by_long = in_scope.iter().find(|entry| entry.long == long);
            let unlisted_short = short.is_some() && by_short.is_none_or(|entry| entry.long != long);
            let missing_short = by_long.is_some_and(|entry| short != Some(entry.short));
            (unlisted_short || missing_short).then(|| ShortFlagUse {
              command: path.clone(),
              long: long.to_owned(),
              short,
            })
          })
        })
        .collect()
    }

    pub(super) fn short_help_deviation(
      command: &Command,
      listed: &[&str],
      max_option_lines: usize,
    ) -> Option<ShortHelpDeviation> {
      let listed: BTreeSet<String> = listed.iter().map(|flag| (*flag).to_owned()).collect();
      let shown = short_help_flags(command);
      let option_lines = option_line_count(command);
      (shown != listed || option_lines > max_option_lines).then(|| ShortHelpDeviation {
        unlisted: shown.difference(&listed).cloned().collect(),
        missing: listed.difference(&shown).cloned().collect(),
        option_lines,
      })
    }

    fn short_help_flags(command: &Command) -> BTreeSet<String> {
      command
        .get_arguments()
        .filter(|arg| !arg.is_hide_set() && !arg.is_hide_short_help_set())
        .map(|arg| format!("--{}", arg.get_long().unwrap_or_else(|| arg.get_id().as_str())))
        .collect()
    }

    fn option_line_count(command: &Command) -> usize {
      command
        .clone()
        .render_help()
        .to_string()
        .lines()
        .filter(|line| line.starts_with("  -") || line.starts_with("      --"))
        .count()
    }
  }
}
