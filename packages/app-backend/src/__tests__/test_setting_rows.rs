#[cfg(test)]
mod tests {
  use crate::command::{AppCommandExt, offered_commands};
  use crate::config::catalog::command_settings;
  use crate::config::choices::ChoiceName;
  use crate::config::rows::RowControl;
  use app_commands::command::AppCommand;
  use itertools::Itertools;
  use pretty_assertions::assert_eq;

  #[test]
  fn test_setting_rows_hold_each_main_setting_once_in_order() {
    for command in offered_commands() {
      let keys = command_settings(command)
        .unwrap()
        .main_rows
        .into_iter()
        .flat_map(|row| row.settings.into_iter().map(|setting| setting.key))
        .collect_vec();
      assert_eq!(command.main_settings().to_vec(), keys, "{command}");
    }
  }

  #[test]
  fn test_setting_rows_timetree() {
    let rows = command_settings(AppCommand::Timetree)
      .unwrap()
      .main_rows
      .into_iter()
      .map(|row| (row.label, row.control))
      .collect_vec();
    let kinds = rows
      .iter()
      .map(|(label, control)| {
        let kind = match control {
          RowControl::Choice { choice } => format!("choice {choice:?}"),
          RowControl::Switches { switches } => {
            format!("switches {}", switches.iter().map(|switch| &switch.key).join(","))
          },
          RowControl::Coalescent { options } => {
            format!("coalescent {}", options.iter().map(|option| &option.value).join(","))
          },
          RowControl::ClockFilter { .. } => "clock-filter".to_owned(),
          RowControl::SubstitutionModel { .. } => "substitution-model".to_owned(),
          RowControl::Setting => "setting".to_owned(),
        };
        format!("{label}: {kind}")
      })
      .collect_vec();
    assert_eq!(
      vec![
        format!("Clock rate: choice {:?}", ChoiceName::ClockRate),
        "Confidence intervals: switches confidence".to_owned(),
        "Coalescent prior: coalescent none,fixed,opt,skyline".to_owned(),
        format!("Root: choice {:?}", ChoiceName::Root),
        "Clock filter: clock-filter".to_owned(),
        "Polytomies: switches keep_polytomies".to_owned(),
        "Substitution model: substitution-model".to_owned(),
        "Iterations: setting".to_owned(),
      ],
      kinds
    );
  }

  #[test]
  fn test_setting_rows_optimize_topology_switches_are_inverted() {
    let switches = command_settings(AppCommand::Optimize)
      .unwrap()
      .main_rows
      .into_iter()
      .find_map(|row| match row.control {
        RowControl::Switches { switches } => Some(switches),
        _ => None,
      })
      .unwrap()
      .into_iter()
      .map(|switch| (switch.key, switch.inverted))
      .collect_vec();
    assert_eq!(
      vec![
        ("no_merge_shared_mutations".to_owned(), true),
        ("no_collapse_short_branches".to_owned(), true),
        ("no_hoist_reversions".to_owned(), true),
      ],
      switches
    );
  }

  #[test]
  fn test_setting_rows_name_no_flags() {
    let flagged = offered_commands()
      .flat_map(|command| command_settings(command).unwrap().main_rows)
      .filter(|row| {
        row.label.contains("--")
          || row.description.as_deref().unwrap_or_default().contains("--")
          || row.settings.iter().any(|setting| setting.label.contains("--"))
      })
      .map(|row| row.label)
      .collect_vec();
    assert_eq!(Vec::<String>::new(), flagged);
  }
}
