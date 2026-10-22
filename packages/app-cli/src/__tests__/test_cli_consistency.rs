#[cfg(test)]
mod tests {
  use crate::__tests__::test_cli_support::tests::helpers::{Pending, is_clap_flag, leaf_commands, treetime_command};
  use helpers::{Property, differences};
  use itertools::Itertools;
  use pretty_assertions::assert_eq;
  use std::collections::BTreeSet;

  const ALLOWED_DIFFERENCES: &[(&str, &[Property], &str)] = &[
    (
      "--max-iter",
      &[Property::Defaults],
      "timetree refines in 2 rounds by default, optimize iterates up to 10 times",
    ),
    (
      "--reroot",
      &[
        Property::ValueType,
        Property::PossibleValues,
        Property::Help,
        Property::LongHelp,
      ],
      "optimize has no dates and offers only min-dev",
    ),
    (
      "--ambiguity-threshold",
      &[Property::LongHelp],
      "the help names the outputs of each command that carry the ambiguity report",
    ),
    (
      "--ambiguity-outputs",
      &[Property::ValueType, Property::PossibleValues],
      "each command lists its outputs that carry the ambiguity report",
    ),
    (
      "--output-selection",
      &[Property::ValueType, Property::PossibleValues],
      "each command lists the outputs it writes",
    ),
    (
      "--output-node-data",
      &[Property::LongHelp],
      "the node data of each command holds the per-node results of that command",
    ),
    (
      "--config",
      &[
        Property::Required,
        Property::Help,
        Property::LongHelp,
        Property::Heading,
      ],
      "pipeline requires its config, which names the pipeline file",
    ),
    (
      "--dry-run",
      &[Property::Help, Property::Heading],
      "pipeline prints the plan of every step instead of the settings of one command",
    ),
    (
      "--output-all",
      &[Property::Help, Property::LongHelp, Property::Heading],
      "pipeline writes the outputs of each step into a subfolder",
    ),
  ];

  #[rustfmt::skip]
  const KNOWN_EXCEPTIONS: &[(&str, &[Property], Pending)] = &[
    ("--confidence",    &[Property::ValueType, Property::Action, Property::NumArgs, Property::ValueNames, Property::Defaults,
                          Property::Aliases, Property::Heading, Property::Help, Property::LongHelp],
                                                                                                 Pending::MugrationConfidenceFlag),
  ];

  #[test]
  fn test_cli_consistency_shared_flags_agree_across_commands() {
    let root = treetime_command();
    let commands = leaf_commands(&root);
    let expected: BTreeSet<String> = ALLOWED_DIFFERENCES
      .iter()
      .map(|(flag, properties, _)| (*flag, *properties))
      .chain(
        KNOWN_EXCEPTIONS
          .iter()
          .map(|(flag, properties, _)| (*flag, *properties)),
      )
      .flat_map(|(flag, properties)| properties.iter().map(move |property| format!("{flag} {property:?}")))
      .collect();

    let actual = differences(&commands, is_clap_flag);

    let actual_keys: BTreeSet<String> = actual.keys().cloned().collect();
    let known = KNOWN_EXCEPTIONS
      .iter()
      .map(|(flag, properties, pending)| format!("{flag} {properties:?} until {pending:?}"))
      .join("\n");
    assert_eq!(
      expected,
      actual_keys,
      "{}\nknown exceptions:\n{known}",
      actual.values().join("\n")
    );
  }

  mod helpers {
    use clap::{Arg, Command};
    use itertools::Itertools;
    use std::collections::BTreeMap;

    #[derive(Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
    pub(super) enum Property {
      ValueType,
      Action,
      NumArgs,
      ValueNames,
      PossibleValues,
      Defaults,
      Required,
      Delimiter,
      Aliases,
      Heading,
      Help,
      LongHelp,
    }

    const PROPERTIES: [Property; 12] = [
      Property::ValueType,
      Property::Action,
      Property::NumArgs,
      Property::ValueNames,
      Property::PossibleValues,
      Property::Defaults,
      Property::Required,
      Property::Delimiter,
      Property::Aliases,
      Property::Heading,
      Property::Help,
      Property::LongHelp,
    ];

    pub(super) fn differences(commands: &[(String, &Command)], skip: fn(&Arg) -> bool) -> BTreeMap<String, String> {
      let by_flag = commands
        .iter()
        .flat_map(|(path, command)| {
          command
            .get_arguments()
            .filter(move |arg| !skip(arg))
            .flat_map(move |arg| {
              flag_names(arg)
                .into_iter()
                .map(move |name| (name, (path.as_str(), arg)))
            })
        })
        .into_group_map();
      by_flag
        .into_iter()
        .filter(|(_, uses)| uses.len() > 1)
        .flat_map(|(flag, uses)| {
          PROPERTIES.into_iter().filter_map(move |property| {
            let values = uses
              .iter()
              .map(|(path, arg)| (property_value(arg, property), *path))
              .into_group_map();
            (values.len() > 1).then(|| {
              let detail = values
                .into_iter()
                .sorted()
                .map(|(value, paths)| format!("  {}: {value}", paths.iter().join(", ")))
                .join("\n");
              (
                format!("{flag} {property:?}"),
                format!("{flag} {property:?}:\n{detail}"),
              )
            })
          })
        })
        .collect()
    }

    fn flag_names(arg: &Arg) -> Vec<String> {
      arg
        .get_long()
        .into_iter()
        .chain(arg.get_all_aliases().unwrap_or_default())
        .map(|name| format!("--{name}"))
        .collect_vec()
    }

    fn property_value(arg: &Arg, property: Property) -> String {
      match property {
        Property::ValueType => format!("{:?}", arg.get_value_parser().type_id()),
        Property::Action => format!("{:?}", arg.get_action()),
        Property::NumArgs => format!("{:?}", arg.get_num_args()),
        Property::ValueNames => format!("{:?}", arg.get_value_names()),
        Property::PossibleValues => format!(
          "{:?}",
          arg
            .get_possible_values()
            .iter()
            .map(|value| value.get_name().to_owned())
            .collect_vec()
        ),
        Property::Defaults => format!("{:?}", arg.get_default_values()),
        Property::Required => format!("{:?}", arg.is_required_set()),
        Property::Delimiter => format!("{:?}", arg.get_value_delimiter()),
        Property::Aliases => format!(
          "{:?}",
          arg
            .get_all_aliases()
            .unwrap_or_default()
            .into_iter()
            .sorted()
            .collect_vec()
        ),
        Property::Heading => format!("{:?}", arg.get_help_heading()),
        Property::Help => format!("{:?}", arg.get_help().map(ToString::to_string)),
        Property::LongHelp => format!("{:?}", arg.get_long_help().map(ToString::to_string)),
      }
    }
  }
}
