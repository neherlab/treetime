#[cfg(test)]
mod tests {
  use crate::__tests__::test_cli_support::tests::helpers::{ANALYSIS, EVERY, find_command, treetime_command};
  use crate::cli::command_help::{HELP_HINT, HELP_HINT_OUTPUTS};
  use helpers::{examples, has_override_usage, headings};
  use itertools::Itertools;
  use pretty_assertions::assert_eq;

  const HEADING_ORDER: &[&str] = &[
    "Config",
    "Input data",
    "Molecular clock",
    "Rooting",
    "Branch lengths",
    "Dating",
    "Polytomies",
    "Coalescent prior",
    "Substitution model",
    "Ancestral reconstruction",
    "Homoplasy",
    "Pruning",
    "Output",
    "Tree ordering",
    "Reproducibility",
    "Parallelism",
    "Verbosity",
  ];

  #[test]
  fn test_cli_help_layout_headings_follow_one_order() {
    let root = treetime_command();
    let expected: Vec<(&str, Vec<&str>)> = vec![];

    let actual = ANALYSIS
      .iter()
      .filter_map(|path| {
        let headings = headings(find_command(&root, path).unwrap());
        let ordered = headings
          .iter()
          .map(|heading| HEADING_ORDER.iter().position(|known| known == heading))
          .collect_vec();
        let sorted = ordered.iter().all(Option::is_some) && ordered.is_sorted();
        (!sorted).then_some((*path, headings))
      })
      .collect_vec();

    assert_eq!(expected, actual);
  }

  #[test]
  fn test_cli_help_layout_every_command_has_usage_and_examples() {
    let root = treetime_command();
    let expected: Vec<(&str, &str)> = vec![];

    let actual = EVERY
      .iter()
      .flat_map(|path| {
        let command = find_command(&root, path).unwrap();
        let hint = if ANALYSIS.contains(path) {
          HELP_HINT_OUTPUTS
        } else {
          HELP_HINT
        };
        let after_help = command.get_after_help().map(ToString::to_string).unwrap_or_default();
        let after_long_help = command
          .get_after_long_help()
          .map(ToString::to_string)
          .unwrap_or_default();
        [
          (!has_override_usage(command)).then_some((*path, "override_usage")),
          (!(2..=3).contains(&examples(&after_long_help).len()))
            .then_some((*path, "two or three examples in after_long_help")),
          (after_help != format!("{after_long_help}\n\n{hint}")).then_some((*path, "after_help with the hint")),
        ]
      })
      .flatten()
      .collect_vec();

    assert_eq!(expected, actual);
  }

  #[test]
  fn test_cli_help_layout_examples_parse() {
    let root = treetime_command();
    let expected: Vec<(String, String)> = vec![];

    let actual = EVERY
      .iter()
      .flat_map(|path| {
        let after_long_help = find_command(&root, path)
          .unwrap()
          .get_after_long_help()
          .map(ToString::to_string)
          .unwrap_or_default();
        examples(&after_long_help)
      })
      .filter_map(|example| {
        let argv = example
          .split_whitespace()
          .take_while(|word| !matches!(*word, ">" | "|"))
          .map(|word| word.trim_matches('\''))
          .collect_vec();
        root
          .clone()
          .try_get_matches_from(argv)
          .err()
          .map(|error| (example, error.kind().to_string()))
      })
      .collect_vec();

    assert_eq!(expected, actual);
  }

  mod helpers {
    use clap::{Arg, Command};
    use itertools::Itertools;

    pub(super) fn headings(command: &Command) -> Vec<&str> {
      command
        .get_arguments()
        .filter_map(Arg::get_help_heading)
        .unique()
        .collect_vec()
    }

    pub(super) fn has_override_usage(command: &Command) -> bool {
      let usage = command.clone().render_usage().to_string();
      let generated = command.clone().override_usage(None::<&str>).render_usage().to_string();
      usage != generated
    }

    pub(super) fn examples(after_long_help: &str) -> Vec<String> {
      after_long_help
        .lines()
        .filter_map(|line| line.strip_prefix("  "))
        .filter(|line| line.starts_with("treetime "))
        .map(ToOwned::to_owned)
        .collect_vec()
    }
  }
}
