use crate::cli::pipeline::suggest::suggestion_suffix;
use eyre::Report;
use itertools::Itertools;
use serde::Serialize;
use serde_json::{Map, Value};
use treetime::commands::ancestral::args::TreetimeAncestralArgs;
use treetime::commands::clock::args::TreetimeClockArgs;
use treetime::commands::mugration::args::TreetimeMugrationArgs;
use treetime::commands::optimize::args::TreetimeOptimizeArgs;
use treetime::commands::prune::args::TreetimePruneArgs;
use treetime::commands::shared::output::CommandKind;
use treetime::commands::timetree::args::TreetimeTimetreeArgs;
use treetime_utils::make_error;

/// Reserved top-level key that associates a JSON config with its schema in an editor.
///
/// Editors read it; the loader ignores it. It is neither a pipeline field nor a command tag, so it
/// must be filtered out before the "exactly one command tag" and unknown-key checks run.
pub const SCHEMA_KEY: &str = "$schema";

/// Command tags accepted inside a pipeline step, in the pipeline's declared order.
///
/// This is the analysis-only subset of the CLI's commands (D2a): the utility commands
/// (`completions`, `schema`, `help-markdown`, `debug`, `arg`) are deliberately absent because they
/// produce no analysis outputs and would validate but do nothing inside a pipeline.
pub const COMMAND_TAGS: [&str; 6] = ["timetree", "optimize", "prune", "ancestral", "clock", "mugration"];

/// A single analysis command invocation within a pipeline.
///
/// Externally tagged and kebab-cased so a step's payload is exactly the command's serialized args
/// object (the same shape the per-command `--config` accepts). `timetree` is boxed to match the
/// CLI enum and keep the variant sizes balanced.
#[derive(Debug, Serialize)]
#[serde(rename_all = "kebab-case")]
pub enum PipelineStepCommand {
  Timetree(Box<TreetimeTimetreeArgs>),
  Optimize(TreetimeOptimizeArgs),
  Prune(TreetimePruneArgs),
  Ancestral(TreetimeAncestralArgs),
  Clock(TreetimeClockArgs),
  Mugration(TreetimeMugrationArgs),
}

impl PipelineStepCommand {
  /// External tag as it appears in the config (kebab-case command name).
  pub fn tag(&self) -> &'static str {
    match self {
      Self::Timetree(_) => "timetree",
      Self::Optimize(_) => "optimize",
      Self::Prune(_) => "prune",
      Self::Ancestral(_) => "ancestral",
      Self::Clock(_) => "clock",
      Self::Mugration(_) => "mugration",
    }
  }

  /// Output taxonomy this command draws from, used to resolve `{{ steps.x.outputs.* }}` chaining.
  pub fn command_kind(&self) -> CommandKind {
    match self {
      Self::Timetree(_) => CommandKind::Timetree,
      Self::Optimize(_) => CommandKind::Optimize,
      Self::Prune(_) => CommandKind::Prune,
      Self::Ancestral(_) => CommandKind::Ancestral,
      Self::Clock(_) => CommandKind::Clock,
      Self::Mugration(_) => CommandKind::Mugration,
    }
  }

  /// Deserialize a command payload once its tag has been isolated.
  ///
  /// The tag is validated against the analysis set with a "did you mean?" hint, so a utility command
  /// or a misspelling is rejected here rather than silently ignored.
  pub fn from_tag_and_value(tag: &str, payload: Value) -> Result<Self, Report> {
    match tag {
      "timetree" => Ok(Self::Timetree(serde_json::from_value(payload)?)),
      "optimize" => Ok(Self::Optimize(serde_json::from_value(payload)?)),
      "prune" => Ok(Self::Prune(serde_json::from_value(payload)?)),
      "ancestral" => Ok(Self::Ancestral(serde_json::from_value(payload)?)),
      "clock" => Ok(Self::Clock(serde_json::from_value(payload)?)),
      "mugration" => Ok(Self::Mugration(serde_json::from_value(payload)?)),
      other => make_error!("unknown command `{other}`; {}", suggestion_suffix(other, &COMMAND_TAGS)),
    }
  }
}

/// One named step: a stable id plus exactly one command invocation.
///
/// The `name` is explicit (not the command name) because `--steps=` selection and
/// `{{ steps.<name>... }}` references need stable ids and must allow the same command twice.
#[derive(Debug, Serialize)]
pub struct PipelineStep {
  pub name: String,
  #[serde(flatten)]
  pub command: PipelineStepCommand,
}

impl PipelineStep {
  /// Parse one step object into a name and a typed command.
  ///
  /// Enforces the step shape directly, because `#[serde(deny_unknown_fields)]` is silently ignored
  /// on flattened structs: an object must carry a string `name` and exactly one command tag (a
  /// reserved `$schema` key is permitted and ignored). Zero tags means the command is missing; more
  /// than one means the step tries to run several commands.
  pub fn from_value(value: Value) -> Result<Self, Report> {
    let Value::Object(mut map) = value else {
      return make_error!("a pipeline step must be a mapping with a `name` and one command");
    };
    map.remove(SCHEMA_KEY);

    let name = match map.remove("name") {
      Some(Value::String(name)) => name,
      Some(_) => return make_error!("pipeline step `name` must be a string"),
      None => return make_error!("pipeline step is missing a `name`"),
    };

    let tags: Vec<String> = map.keys().cloned().collect();
    let (tag, payload) = match tags.as_slice() {
      [] => return make_error!("pipeline step `{name}` has no command; expected one of {}", commands_list()),
      [tag] => {
        let payload = map.remove(tag).unwrap_or(Value::Null);
        (tag.clone(), payload)
      },
      _ => {
        return make_error!(
          "pipeline step `{name}` has more than one command ({}); a step runs exactly one command",
          tags.iter().sorted().map(|tag| format!("`{tag}`")).join(", ")
        );
      },
    };

    let command = PipelineStepCommand::from_tag_and_value(&tag, payload)
      .map_err(|err| eyre::eyre!("in pipeline step `{name}`: {err}"))?;

    Ok(Self { name, command })
  }
}

/// The whole pipeline in typed form, used for schema generation and for dumping an example config.
///
/// The loader does not deserialize into this type directly: `vars`, `output_all`, and each step's
/// outputs are resolved in a staged, backward-only pass (interpolation depends on earlier results),
/// after which the typed steps are assembled. This type fixes the on-disk shape that the staged pass
/// and the generated schema must agree on.
#[derive(Debug, Serialize)]
pub struct Pipeline {
  #[serde(rename = "$schema", skip_serializing_if = "Option::is_none")]
  pub schema_ref: Option<String>,

  #[serde(default, skip_serializing_if = "Map::is_empty")]
  pub vars: Map<String, Value>,

  #[serde(skip_serializing_if = "Option::is_none")]
  pub output_all: Option<String>,

  pub steps: Vec<PipelineStep>,
}

/// Human-readable command list for error messages, in declared order.
fn commands_list() -> String {
  COMMAND_TAGS.iter().map(|tag| format!("`{tag}`")).collect::<Vec<_>>().join(", ")
}

#[cfg(test)]
mod tests {
  use super::*;
  use pretty_assertions::assert_eq;
  use serde_json::json;
  use treetime_utils::assert_error;

  fn method_anc(command: &PipelineStepCommand) -> Option<String> {
    match command {
      PipelineStepCommand::Ancestral(args) => Some(format!("{:?}", args.method_anc)),
      _ => None,
    }
  }

  // A well-formed ancestral step parses, snake_case fields land on the typed args, and the tag maps
  // to the right command kind.
  #[test]
  fn test_types_step_parses_ancestral_with_snake_case_fields() {
    let value = json!({
      "name": "anc",
      "ancestral": { "tree": "t.nwk", "method_anc": "marginal", "dense": true }
    });
    let step = PipelineStep::from_value(value).unwrap();
    assert_eq!("anc", step.name);
    assert_eq!("ancestral", step.command.tag());
    assert_eq!(CommandKind::Ancestral, step.command.command_kind());
    assert_eq!(Some("Marginal".to_owned()), method_anc(&step.command));
  }

  // The kebab command tag is what selects the command; the same command may appear twice under
  // different step names.
  #[test]
  fn test_types_step_parses_timetree_tag() {
    let step = PipelineStep::from_value(json!({ "name": "tt", "timetree": { "clock_rate": 0.003 } })).unwrap();
    assert_eq!("timetree", step.command.tag());
  }

  // A utility command (not in the analysis set) is rejected with a did-you-mean among analysis
  // commands, never silently accepted.
  #[test]
  fn test_types_step_rejects_utility_command_tag() {
    let result = PipelineStep::from_value(json!({ "name": "x", "debug": {} }));
    assert_error!(
      result,
      "in pipeline step `x`: unknown command `debug`; valid values: `ancestral`, `clock`, `mugration`, `optimize`, `prune`, `timetree`"
    );
  }

  // A misspelled command tag yields a "did you mean?" pointing at the closest analysis command.
  #[test]
  fn test_types_step_suggests_closest_command_for_typo() {
    let result = PipelineStep::from_value(json!({ "name": "x", "timtree": {} }));
    assert_error!(
      result,
      "in pipeline step `x`: unknown command `timtree`; did you mean `timetree`? Valid values: `ancestral`, `clock`, `mugration`, `optimize`, `prune`, `timetree`"
    );
  }

  // A step with two command tags is rejected: a step runs exactly one command.
  #[test]
  fn test_types_step_rejects_multiple_commands() {
    let result = PipelineStep::from_value(json!({ "name": "x", "timetree": {}, "clock": {} }));
    assert_error!(
      result,
      "pipeline step `x` has more than one command (`clock`, `timetree`); a step runs exactly one command"
    );
  }

  // A step with no command tag is rejected.
  #[test]
  fn test_types_step_rejects_missing_command() {
    let result = PipelineStep::from_value(json!({ "name": "x" }));
    assert_error!(
      result,
      "pipeline step `x` has no command; expected one of `timetree`, `optimize`, `prune`, `ancestral`, `clock`, `mugration`"
    );
  }

  // A step without a name is rejected.
  #[test]
  fn test_types_step_rejects_missing_name() {
    let result = PipelineStep::from_value(json!({ "timetree": {} }));
    assert_error!(result, "pipeline step is missing a `name`");
  }

  // A reserved `$schema` key is ignored, not treated as a command tag.
  #[test]
  fn test_types_step_ignores_schema_key() {
    let value = json!({ "$schema": "./pipeline.schema.json", "name": "tt", "timetree": {} });
    let step = PipelineStep::from_value(value).unwrap();
    assert_eq!("timetree", step.command.tag());
  }
}
