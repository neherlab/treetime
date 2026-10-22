use crate::check_inputs::InputSlot;
use crate::command::AppCommandExt;
use crate::config::choices::{SettingChoice, setting_choices};
use crate::config::labels::setting_label;
use app_commands::command::AppCommand;
use app_commands::config::cli_flags::setting_arg;
use app_commands::config::cli_rules::setting_conflicts;
use app_commands::config::properties::{DEFS_PREFIX, PathRole, SettingProperty, def_pointer, setting_properties};
use app_commands::json_value::JsonValue;
use clap::Arg;
use deser::Serialize;
use eyre::Report;
use itertools::{Itertools, izip};
use schemars::JsonSchema;
use serde_json::{Map, Value};
use strum::IntoEnumIterator;
use treetime_schema::UNSET_KEY;
use treetime_schema::skip_serializing_optionals;
use treetime_utils::{make_error, make_report};

pub fn setting_catalog() -> Result<SettingCatalog, Report> {
  Ok(SettingCatalog {
    commands: AppCommand::iter().map(command_settings).try_collect()?,
  })
}

pub fn command_settings(command: AppCommand) -> Result<CommandSettings, Report> {
  let schema = command.config_schema().to_value();
  let cli = command.cli_command();
  let defaults = command.default_config()?;
  let properties = setting_properties(&schema)?;

  let args: Vec<&Arg> = properties
    .iter()
    .map(|property| setting_arg(&cli, property))
    .try_collect()?;
  let mut settings = izip!(&properties, &args)
    .map(|(property, arg)| {
      let order = cli
        .get_arguments()
        .position(|candidate| candidate.get_id() == arg.get_id())
        .unwrap_or(usize::MAX);
      let spec = setting_spec(
        &schema,
        property,
        arg,
        &defaults,
        setting_conflicts(&cli, &properties, &args, arg),
      )?;
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

  let specs = settings.iter().map(|(_, _, spec)| spec).collect_vec();
  let inputs = command
    .inputs()
    .iter()
    .map(|input| {
      let list = specs
        .iter()
        .any(|spec| spec.key == input.kind.setting() && spec.kind == SettingKind::List);
      InputSlot::new(*input, list)
    })
    .collect();
  let settings = settings.into_iter().map(|(_, _, spec)| spec).collect_vec();
  Ok(CommandSettings {
    command,
    title: command.title().to_owned(),
    description: command.description().to_owned(),
    inputs,
    uses_dates: command.uses_dates(),
    groups,
    main_settings: command.main_settings().iter().map(|key| (*key).to_owned()).collect(),
    choices: setting_choices(&settings)?,
    settings,
  })
}

/// Settings of every command the app runs, as the settings form shows them.
#[derive(Clone, Debug, PartialEq, JsonSchema, Serialize)]
pub struct SettingCatalog {
  /// Settings of each command.
  pub commands: Vec<CommandSettings>,
}

/// Settings of one command.
#[derive(Clone, Debug, PartialEq, JsonSchema, Serialize)]
pub struct CommandSettings {
  /// The command.
  pub command: AppCommand,
  /// Name of the analysis as the app shows it, for example `Time tree`.
  pub title: String,
  /// One sentence on what the analysis does.
  pub description: String,
  /// Input files the command reads, in the order the form asks for them.
  pub inputs: Vec<InputSlot>,
  /// Whether the command reads sampling dates from the metadata.
  pub uses_dates: bool,
  /// Setting groups, in the order `treetime <command> --help` lists their headings.
  pub groups: Vec<String>,
  /// Keys of the settings the form shows first, in the order it shows them.
  pub main_settings: Vec<String>,
  /// Settings that the form shows as one choice between options.
  pub choices: Vec<SettingChoice>,
  /// Every setting of the command, by group, in the order `--help` lists them.
  pub settings: Vec<SettingSpec>,
}

fn setting_spec(
  schema: &Value,
  setting: &SettingProperty,
  arg: &Arg,
  defaults: &Map<String, Value>,
  conflicts: Vec<String>,
) -> Result<SettingSpec, Report> {
  let key = setting.key.clone();
  let property = schema
    .pointer(&setting.schema_pointer)
    .ok_or_else(|| make_report!("schema has no property at `{}`", setting.schema_pointer))?;
  let long = arg
    .get_long()
    .ok_or_else(|| make_report!("config key `{key}` has no command-line flag"))?;
  let group = arg
    .get_help_heading()
    .ok_or_else(|| make_report!("config key `{key}` has no help heading"))?;
  let description = Description::parse(property.get("description").and_then(Value::as_str).unwrap_or(""));
  let nullable = is_nullable(property);
  let base = deref(schema, property)?;
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
    label: setting_label(&key),
    flag: format!("--{long}"),
    group: group.to_owned(),
    role: setting.path_role.into(),
    kind,
    nullable,
    options,
    item_kind,
    default_value: defaults.get(&key).cloned().map(JsonValue),
    minimum: base.get("minimum").and_then(Value::as_f64),
    examples: property
      .get("examples")
      .and_then(Value::as_array)
      .map(|examples| examples.iter().cloned().map(JsonValue).collect())
      .unwrap_or_default(),
    value_names: arg
      .get_value_names()
      .unwrap_or_default()
      .iter()
      .map(ToString::to_string)
      .collect(),
    conflicts,
    help: description.help,
    more: description.more,
    key,
  })
}

/// One setting of a command.
#[derive(Clone, Debug, PartialEq, JsonSchema, Serialize)]
#[deser(skip_serializing_optionals)]
#[schemars(transform = skip_serializing_optionals)]
pub struct SettingSpec {
  /// Key of the setting in the configuration, for example `clock_rate`.
  pub key: String,
  /// Name of the setting as the form shows it, for example `Clock rate std. dev.`.
  pub label: String,
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
  /// Value of the setting when the configuration does not set it; absent when the setting has no default.
  pub default_value: Option<JsonValue>,
  /// Smallest allowed value of a number, when there is one.
  pub minimum: Option<f64>,
  /// Typical values of the setting, as the configuration spells them; the first one is the value the form fills in
  /// when the setting is turned on.
  pub examples: Vec<JsonValue>,
  /// Names of the values the command-line flag takes, in order, for example `SLACK` and `COUPLING`.
  pub value_names: Vec<String>,
  /// Keys of the settings this setting cannot be used together with, as the command line rejects them.
  pub conflicts: Vec<String>,
  /// First paragraph of the setting's description.
  pub help: String,
  /// Remaining paragraphs of the setting's description.
  pub more: String,
}

/// What a setting names.
#[derive(Clone, Copy, Debug, PartialEq, Eq, JsonSchema, Serialize)]
#[schemars(rename_all = "kebab-case")]
#[deser(rename_all = "kebab-case")]
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
#[derive(Clone, Copy, Debug, PartialEq, Eq, JsonSchema, Serialize)]
#[schemars(rename_all = "kebab-case")]
#[deser(rename_all = "kebab-case")]
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

/// Type of the items of a list setting.
#[derive(Clone, Copy, Debug, PartialEq, Eq, JsonSchema, Serialize)]
#[schemars(rename_all = "kebab-case")]
#[deser(rename_all = "kebab-case")]
pub enum ListItemKind {
  String,
  Number,
  Integer,
}

fn is_nullable(property: &Value) -> bool {
  property.get(UNSET_KEY) == Some(&Value::Bool(true))
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

/// One allowed value of a setting.
#[derive(Clone, Debug, PartialEq, Eq, JsonSchema, Serialize)]
pub struct SettingOption {
  /// The value, as the configuration spells it.
  pub value: String,
  /// Description of the value, when it has one.
  pub help: String,
}

fn deref<'a>(schema: &'a Value, node: &'a Value) -> Result<&'a Value, Report> {
  match node.get("$ref").and_then(Value::as_str) {
    None => Ok(node),
    Some(reference) => {
      let name = reference
        .strip_prefix(DEFS_PREFIX)
        .ok_or_else(|| make_report!("schema reference `{reference}` is not a local definition"))?;
      schema
        .pointer(&def_pointer(name))
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
