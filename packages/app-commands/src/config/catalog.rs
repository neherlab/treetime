use crate::check_inputs::CommandInput;
use crate::command::AppCommand;
use crate::config::properties::{LeafProperty, PathRole, leaf_properties};
use crate::config::settings::setting_ref;
use crate::config::source::escape_pointer;
use clap::{Arg, Command};
use eyre::Report;
use itertools::Itertools;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use serde_json::Value;
use strum::IntoEnumIterator;
use treetime_utils::{make_error, make_report};

const DEFS_PREFIX: &str = "#/$defs/";

pub fn setting_catalog() -> Result<SettingCatalog, Report> {
  Ok(SettingCatalog {
    commands: AppCommand::iter().map(command_settings).try_collect()?,
  })
}

pub fn command_settings(command: AppCommand) -> Result<CommandSettings, Report> {
  let schema = command.config_schema().to_value();
  let cli = command.cli_command();
  let defaults = Value::Object(command.default_config()?);
  let leaves = leaf_properties(&schema)?;

  let mut settings = leaves
    .iter()
    .map(|leaf| {
      let arg = setting_arg(&cli, leaf)?;
      let order = cli
        .get_arguments()
        .position(|candidate| candidate.get_id() == arg.get_id())
        .unwrap_or(usize::MAX);
      let spec = setting_spec(&schema, leaf, arg, &defaults)?;
      Ok((order, arg.get_display_order(), spec))
    })
    .collect::<Result<Vec<_>, Report>>()?;

  let groups = settings
    .iter()
    .sorted_by_key(|(order, _, _)| *order)
    .map(|(_, _, spec)| spec.group.clone())
    .unique()
    .collect_vec();
  settings.sort_by_key(|(order, display_order, spec)| {
    (
      groups.iter().position(|group| *group == spec.group),
      *display_order,
      *order,
    )
  });

  Ok(CommandSettings {
    command,
    inputs: command.inputs().to_vec(),
    uses_dates: command.uses_dates(),
    groups,
    settings: settings.into_iter().map(|(_, _, spec)| spec).collect(),
  })
}

pub fn setting_arg<'a>(cli: &'a Command, leaf: &LeafProperty) -> Result<&'a Arg, Report> {
  cli
    .get_arguments()
    .find(|arg| arg.get_id() == leaf.key() && arg.get_long().is_some())
    .ok_or_else(|| make_report!("config key `{}` has no command-line flag", leaf.key_path.join(".")))
}

/// Settings of every command the app runs, as the settings form shows them.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
pub struct SettingCatalog {
  /// Settings of each command.
  pub commands: Vec<CommandSettings>,
}

/// Settings of one command.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
pub struct CommandSettings {
  /// The command.
  pub command: AppCommand,
  /// Input files the command reads, in the order the form asks for them.
  pub inputs: Vec<CommandInput>,
  /// Whether the command reads sampling dates from the metadata.
  pub uses_dates: bool,
  /// Setting groups, in the order `treetime <command> --help` lists their headings.
  pub groups: Vec<String>,
  /// Every setting of the command, by group, in the order `--help` lists them.
  pub settings: Vec<SettingSpec>,
}

/// One setting of a command.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
pub struct SettingSpec {
  /// Key path of the setting joined with `.`, for example `branch_split.method`.
  pub key: String,
  /// Key path of the setting in the configuration.
  pub path: Vec<String>,
  /// Command-line flag, for example `--clock-rate`.
  pub flag: String,
  /// Heading of the setting in `--help`.
  pub group: String,
  /// What the setting names: a file, or a value.
  pub role: SettingRole,
  /// Form of the value, which selects the control that edits it.
  pub kind: SettingKind,
  /// Whether the setting can be unset.
  pub nullable: bool,
  /// Allowed values of an `enum` or `enum-list` setting.
  pub options: Vec<SettingOption>,
  /// Type of the items of a `list` setting.
  pub item_kind: ListItemKind,
  /// Value of the setting when the configuration does not set it.
  pub default_value: Value,
  /// Smallest allowed value of a number, when there is one.
  pub minimum: Option<f64>,
  /// First paragraph of the setting's description.
  pub help: String,
  /// Remaining paragraphs of the setting's description.
  pub more: String,
}

/// What a setting names.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
#[serde(rename_all = "kebab-case")]
pub enum SettingRole {
  /// A value that controls the analysis.
  Setting,
  /// Path of an input file.
  Input,
  /// Path template of input files, with a placeholder.
  InputTemplate,
  /// Path of an output file.
  Output,
}

impl From<Option<PathRole>> for SettingRole {
  fn from(role: Option<PathRole>) -> Self {
    match role {
      None => Self::Setting,
      Some(PathRole::Input) => Self::Input,
      Some(PathRole::InputTemplate) => Self::InputTemplate,
      Some(PathRole::Output) => Self::Output,
    }
  }
}

/// Form of a setting's value.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
#[serde(rename_all = "kebab-case")]
pub enum SettingKind {
  /// On or off.
  Switch,
  /// On, off, or unset.
  Tristate,
  /// One of a set of values.
  Enum,
  /// A whole number.
  Integer,
  /// A number.
  Number,
  /// Text.
  Text,
  /// A list of text or numbers.
  List,
  /// A list of values from a set.
  EnumList,
}

/// Type of the items of a list setting.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
#[serde(rename_all = "kebab-case")]
pub enum ListItemKind {
  String,
  Number,
  Integer,
}

/// One allowed value of a setting.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
pub struct SettingOption {
  /// The value, as the configuration spells it.
  pub value: String,
  /// Description of the value, when it has one.
  pub help: String,
}

fn setting_spec(schema: &Value, leaf: &LeafProperty, arg: &Arg, defaults: &Value) -> Result<SettingSpec, Report> {
  let key = leaf.key_path.join(".");
  let property = schema
    .pointer(&leaf.schema_pointer)
    .ok_or_else(|| make_report!("schema has no property at `{}`", leaf.schema_pointer))?;
  let long = arg
    .get_long()
    .ok_or_else(|| make_report!("config key `{key}` has no command-line flag"))?;
  let group = arg
    .get_help_heading()
    .ok_or_else(|| make_report!("config key `{key}` has no help heading"))?;
  let description = Description::parse(property.get("description").and_then(Value::as_str).unwrap_or(""));
  let nullable = is_nullable(property);
  let base = non_null_branch(schema, property)?;
  let options = enum_options(schema, base)?;
  let types = type_list(base);

  let (kind, options, item_kind) = if !options.is_empty() {
    (SettingKind::Enum, options, ListItemKind::String)
  } else if types.contains(&"array") {
    let items = base.get("items").map(|items| deref(schema, items)).transpose()?;
    let item_options = items
      .map(|items| enum_options(schema, items))
      .transpose()?
      .unwrap_or_default();
    if item_options.is_empty() {
      let item_types = items.map(type_list).unwrap_or_default();
      (SettingKind::List, vec![], list_item_kind(&item_types, &key)?)
    } else {
      (SettingKind::EnumList, item_options, ListItemKind::String)
    }
  } else if types.contains(&"boolean") {
    let kind = if nullable {
      SettingKind::Tristate
    } else {
      SettingKind::Switch
    };
    (kind, vec![], ListItemKind::String)
  } else if types.contains(&"integer") {
    (SettingKind::Integer, vec![], ListItemKind::String)
  } else if types.contains(&"number") {
    (SettingKind::Number, vec![], ListItemKind::String)
  } else if types.contains(&"string") {
    (SettingKind::Text, vec![], ListItemKind::String)
  } else {
    return make_error!("setting `{key}` has a schema construct the settings form does not render");
  };

  Ok(SettingSpec {
    key,
    path: leaf.key_path.clone(),
    flag: format!("--{long}"),
    group: group.to_owned(),
    role: leaf.path_role.into(),
    kind,
    nullable,
    options,
    item_kind,
    default_value: setting_ref(
      defaults
        .as_object()
        .ok_or_else(|| make_report!("the defaults of a command must be a mapping of settings"))?,
      &leaf.key_path,
    )
    .cloned()
    .unwrap_or(Value::Null),
    minimum: base.get("minimum").and_then(Value::as_f64),
    help: description.help,
    more: description.more,
  })
}

fn list_item_kind(types: &[&str], key: &str) -> Result<ListItemKind, Report> {
  if types.contains(&"integer") {
    Ok(ListItemKind::Integer)
  } else if types.contains(&"number") {
    Ok(ListItemKind::Number)
  } else if types.contains(&"string") {
    Ok(ListItemKind::String)
  } else {
    make_error!("list setting `{key}` has items the settings form does not render")
  }
}

fn is_nullable(property: &Value) -> bool {
  type_list(property).contains(&"null") || alternatives(property).any(|branch| type_list(branch).contains(&"null"))
}

fn non_null_branch<'a>(schema: &'a Value, property: &'a Value) -> Result<&'a Value, Report> {
  let branches = property.get("anyOf").and_then(Value::as_array).map(|branches| {
    branches
      .iter()
      .filter(|branch| !type_list(branch).contains(&"null"))
      .collect_vec()
  });
  match branches.as_deref() {
    Some([branch]) => deref(schema, branch),
    _ => deref(schema, property),
  }
}

fn enum_options(schema: &Value, node: &Value) -> Result<Vec<SettingOption>, Report> {
  let target = deref(schema, node)?;
  let direct = target
    .get("enum")
    .and_then(Value::as_array)
    .into_iter()
    .flatten()
    .filter_map(Value::as_str)
    .map(|value| SettingOption {
      value: value.to_owned(),
      help: String::new(),
    });
  let constant = target.get("const").and_then(Value::as_str).map(|value| SettingOption {
    value: value.to_owned(),
    help: Description::parse(target.get("description").and_then(Value::as_str).unwrap_or("")).help,
  });
  let mut options = direct.chain(constant).collect_vec();
  for branch in target.get("oneOf").and_then(Value::as_array).into_iter().flatten() {
    options.extend(enum_options(schema, branch)?);
  }
  Ok(options)
}

fn deref<'a>(schema: &'a Value, node: &'a Value) -> Result<&'a Value, Report> {
  match node.get("$ref").and_then(Value::as_str) {
    None => Ok(node),
    Some(reference) => {
      let name = reference
        .strip_prefix(DEFS_PREFIX)
        .ok_or_else(|| make_report!("schema reference `{reference}` is not a local definition"))?;
      schema
        .pointer(&format!("/$defs/{}", escape_pointer(name)))
        .ok_or_else(|| make_report!("schema has no definition `{name}`"))
    },
  }
}

fn type_list(node: &Value) -> Vec<&str> {
  match node.get("type") {
    Some(Value::String(single)) => vec![single.as_str()],
    Some(Value::Array(types)) => types.iter().filter_map(Value::as_str).collect(),
    _ => vec![],
  }
}

fn alternatives(node: &Value) -> impl Iterator<Item = &Value> {
  ["anyOf", "oneOf"]
    .into_iter()
    .filter_map(|key| node.get(key).and_then(Value::as_array))
    .flatten()
}

struct Description {
  help: String,
  more: String,
}

impl Description {
  fn parse(description: &str) -> Self {
    let paragraphs = description
      .lines()
      .map(str::trim)
      .chunk_by(|line| line.is_empty())
      .into_iter()
      .filter(|(blank, _)| !blank)
      .map(|(_, mut lines)| lines.join(" "))
      .collect_vec();
    Self {
      help: paragraphs.first().cloned().unwrap_or_default(),
      more: paragraphs.iter().skip(1).join("\n\n"),
    }
  }
}
