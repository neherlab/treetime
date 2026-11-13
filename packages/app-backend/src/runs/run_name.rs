use crate::check_inputs::InputKind;
use crate::config::catalog::{CommandSettings, SettingPurpose};
use crate::config::effective::value_text;
use crate::config::rows::RowControl;
use crate::runs::record::RunInput;
use crate::runs::setting_differences::{DifferenceValues, SettingDifference};
use app_commands::json_value::JsonValue;
use itertools::Itertools;
use serde_json::{Map, Value};
use std::path::Path;

const NAMED_CHANGES: usize = 3;

const UNSET_VALUE: &str = "none";

pub fn run_name(dataset: Option<&str>, catalog: &CommandSettings, differences: &[SettingDifference]) -> String {
  let head = match dataset {
    Some(dataset) => format!("{dataset}, {}", lowercase_first(&catalog.title)),
    None => catalog.title.clone(),
  };
  let changes = differences
    .iter()
    .filter(|difference| difference.purpose == SettingPurpose::Analysis)
    .filter_map(|difference| match &difference.values {
      DifferenceValues::Setting { second, .. } => Some(change_text(catalog, difference, second.as_ref())),
      DifferenceValues::Input { .. } => None,
    })
    .collect_vec();
  if changes.is_empty() {
    return head;
  }
  let named = changes.iter().take(NAMED_CHANGES).join(", ");
  match changes.len().saturating_sub(NAMED_CHANGES) {
    0 => format!("{head}: {named}"),
    more => format!("{head}: {named} +{more}"),
  }
}

pub fn dataset_label(example: Option<&str>, inputs: &[RunInput]) -> Option<String> {
  let tree = inputs
    .iter()
    .find(|input| input.setting == InputKind::Tree.setting())
    .map(|input| input.path.as_path());
  label_of(example, tree)
}

pub fn draft_dataset_label(example: Option<&str>, config: Option<&Map<String, Value>>) -> Option<String> {
  let tree = config
    .and_then(|config| config.get(InputKind::Tree.setting()))
    .and_then(Value::as_str)
    .map(Path::new);
  label_of(example, tree)
}

fn label_of(example: Option<&str>, tree: Option<&Path>) -> Option<String> {
  example.map(ToOwned::to_owned).or_else(|| {
    tree
      .and_then(Path::file_name)
      .map(|name| name.to_string_lossy().into_owned())
  })
}

fn change_text(catalog: &CommandSettings, difference: &SettingDifference, value: Option<&JsonValue>) -> String {
  let label = lowercase_first(&difference.label);
  match value.map(|value| &value.0) {
    Some(Value::Bool(true)) => label,
    Some(Value::Bool(false)) => format!("{label} off"),
    value => {
      let text = value.map_or_else(|| UNSET_VALUE.to_owned(), |value| value_text(Some(value)));
      option_label(catalog, &difference.key, &text).map_or_else(|| format!("{label} {text}"), lowercase_first)
    },
  }
}

fn option_label<'a>(catalog: &'a CommandSettings, key: &str, value: &str) -> Option<&'a str> {
  catalog
    .main_rows
    .iter()
    .filter(|row| row.settings.first().is_some_and(|setting| setting.key == key))
    .find_map(|row| match &row.control {
      RowControl::Coalescent { options, .. } => options
        .iter()
        .find(|option| option.value == value)
        .map(|option| option.label.as_str()),
      _ => None,
    })
}

fn lowercase_first(text: &str) -> String {
  let mut chars = text.chars();
  match (chars.next(), chars.next()) {
    (Some(_), Some(second)) if second.is_uppercase() => text.to_owned(),
    (Some(first), _) => first.to_lowercase().chain(text.chars().skip(1)).collect(),
    (None, _) => String::new(),
  }
}
