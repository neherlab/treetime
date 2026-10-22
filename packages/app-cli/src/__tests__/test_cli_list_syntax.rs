#[cfg(test)]
mod tests {
  use crate::__tests__::test_cli_support::tests::helpers::{all_commands, treetime_command};
  use helpers::{Arity, ListSyntaxDeviation, ListWithoutCommas, list_syntax_deviations};
  use itertools::Itertools;
  use pretty_assertions::assert_eq;

  #[rustfmt::skip]
  const LISTS_WITHOUT_COMMAS: &[ListWithoutCommas] = &[
    ListWithoutCommas { flag: "--alignment",           values: Arity::OnePerFlag, reason: "a path can contain a comma" },
    ListWithoutCommas { flag: "--model-params",        values: Arity::OneOrMore,  reason: "an item such as pis=0.25,0.25,0.25,0.25 contains commas" },
    ListWithoutCommas { flag: "--metadata-delimiters", values: Arity::OneOrMore,  reason: "a comma is one of its values" },
    ListWithoutCommas { flag: "--relax",               values: Arity::Exactly(2), reason: "it takes exactly two numbers" },
  ];

  #[test]
  fn test_cli_list_syntax_lists_take_repeated_space_and_comma_separated_values() {
    let root = treetime_command();
    let commands = all_commands(&root);
    let expected: Vec<ListSyntaxDeviation> = vec![];

    let actual = list_syntax_deviations(&commands, LISTS_WITHOUT_COMMAS);

    let rules = LISTS_WITHOUT_COMMAS
      .iter()
      .map(|list| format!("{} {:?}: {}", list.flag, list.values, list.reason))
      .join("\n");
    assert_eq!(expected, actual, "lists without comma-separated values:\n{rules}");
  }

  #[test]
  fn test_cli_list_syntax_lists_without_commas_name_existing_flags() {
    let root = treetime_command();
    let flags = all_commands(&root)
      .into_iter()
      .flat_map(|(_, command)| command.get_arguments().filter_map(|arg| arg.get_long()).collect_vec())
      .map(|long| format!("--{long}"))
      .collect_vec();
    let expected: Vec<&str> = vec![];

    let actual = LISTS_WITHOUT_COMMAS
      .iter()
      .map(|list| list.flag)
      .filter(|flag| !flags.iter().any(|existing| existing == flag))
      .collect_vec();

    assert_eq!(expected, actual);
  }

  mod helpers {
    use clap::builder::ValueRange;
    use clap::{ArgAction, Command};
    use itertools::Itertools;

    pub(super) struct ListWithoutCommas {
      pub(super) flag: &'static str,
      pub(super) values: Arity,
      pub(super) reason: &'static str,
    }

    #[derive(Copy, Clone, Debug)]
    pub(super) enum Arity {
      OnePerFlag,
      OneOrMore,
      Exactly(usize),
    }

    impl Arity {
      fn range(self) -> ValueRange {
        match self {
          Self::OnePerFlag => ValueRange::SINGLE,
          Self::OneOrMore => ValueRange::new(1..),
          Self::Exactly(count) => ValueRange::new(count),
        }
      }
    }

    #[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
    pub(super) struct ListSyntaxDeviation {
      pub(super) command: String,
      pub(super) flag: String,
      pub(super) syntax: String,
    }

    pub(super) fn list_syntax_deviations(
      commands: &[(String, &Command)],
      without_commas: &[ListWithoutCommas],
    ) -> Vec<ListSyntaxDeviation> {
      commands
        .iter()
        .flat_map(|(path, command)| {
          command
            .get_arguments()
            .filter(|arg| matches!(arg.get_action(), ArgAction::Append))
            .filter_map(move |arg| {
              let flag = format!("--{}", arg.get_long()?);
              let expected = without_commas
                .iter()
                .find(|list| list.flag == flag)
                .map_or_else(|| (Some(','), ValueRange::new(1..)), |list| (None, list.values.range()));
              let actual = (
                arg.get_value_delimiter(),
                arg.get_num_args().unwrap_or(ValueRange::SINGLE),
              );
              (expected != actual).then(|| ListSyntaxDeviation {
                command: path.clone(),
                flag,
                syntax: format!("{actual:?}"),
              })
            })
        })
        .sorted()
        .dedup()
        .collect_vec()
    }
  }
}
