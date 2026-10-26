#[cfg(test)]
mod tests {
  use crate::__tests__::test_cli_support::tests::helpers::{find_command, treetime_command};
  use app_commands::command::AppCommand;
  use app_commands::config::properties::{PathRole, setting_properties};
  use itertools::Itertools;
  use pretty_assertions::assert_eq;
  use strum::IntoEnumIterator;
  use treetime_utils::o;

  const COMPRESSED_INPUT: &str = "A path that ends in `.gz`, `.bz2`, `.xz`, or `.zst` reads a compressed file.";

  #[test]
  fn test_cli_input_help_names_the_compressed_formats_on_every_input_flag() {
    let root = treetime_command();
    let expected: Vec<String> = vec![];

    let actual = AppCommand::iter()
      .flat_map(|command| {
        let inputs = setting_properties(command.config_schema().as_value())
          .expect("the config schema of every command lists its settings")
          .into_iter()
          .filter(|property| matches!(property.path_role, Some(PathRole::Input | PathRole::InputTemplate)))
          .map(|property| property.key)
          .chain([o!("config")]);
        let clap_command = find_command(&root, &command.to_string()).expect("every app command is a CLI command");
        inputs
          .filter(|key| {
            let arg = clap_command
              .get_arguments()
              .find(|arg| arg.get_id().as_str() == key)
              .expect("every input setting is a flag");
            let help = arg
              .get_long_help()
              .or_else(|| arg.get_help())
              .map(ToString::to_string)
              .unwrap_or_default();
            !help.contains(COMPRESSED_INPUT)
          })
          .map(|key| format!("{command} --{}", key.replace('_', "-")))
          .collect_vec()
      })
      .collect_vec();

    assert_eq!(expected, actual);
  }
}
