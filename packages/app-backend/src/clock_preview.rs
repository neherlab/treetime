use crate::command::config_from_request;
use crate::config::catalog::{SettingPurpose, command_settings};
use crate::results::clock::{ClockResults, clock_results_of_tips};
use crate::runs::errors::invalid;
use app_commands::command::AppCommand;
use app_commands::command_config::CommandConfig;
use app_commands::commands::clock::run::{ClockRegression, clock_regression};
use app_commands::json_value::SparseConfig;
use deser::{Deserialize, Serialize};
use eyre::Report;
use schemars::JsonSchema;
use serde_json::Map;
use std::collections::BTreeSet;
use treetime::cancel::NoopCancel;
use treetime::progress::NoopProgress;
use treetime_io::tree::tree_read_file;
use treetime_schema::skip_serializing_optionals;
use treetime_utils::make_internal_error;

pub const CLOCK_PREVIEW_TIP_LIMIT: usize = 500;

pub fn clock_preview(request: &ClockPreviewRequest) -> Result<ClockPreview, Report> {
  let config = clock_config(request.command, &request.config)?;
  let CommandConfig::Clock { config } = config_from_request(AppCommand::Clock, &config)? else {
    return make_internal_error!("a clock configuration reads as the configuration of another command");
  };
  let Some(tree) = &config.tree.tree else {
    return Err(invalid("a clock preview needs a tree"));
  };
  if config.metadata.metadata.is_none() {
    return Err(invalid("a clock preview needs the metadata with the sampling dates"));
  }
  let tips = tree_read_file(tree, config.newick_dialect.newick_dialect.dialect())?
    .graph
    .num_leaves();
  if tips > CLOCK_PREVIEW_TIP_LIMIT && !request.beyond_limit {
    return Ok(ClockPreview {
      tips,
      tip_limit: CLOCK_PREVIEW_TIP_LIMIT,
      clock: None,
    });
  }
  let ClockRegression {
    clock_model,
    regression_results,
    tip_names,
  } = clock_regression(*config, &NoopCancel, &NoopProgress, &NoopProgress)?;
  let tip_names: BTreeSet<&str> = tip_names.iter().map(String::as_str).collect();
  Ok(ClockPreview {
    tips,
    tip_limit: CLOCK_PREVIEW_TIP_LIMIT,
    clock: Some(clock_results_of_tips(
      Some(&tip_names),
      Some(&clock_model),
      Some(&regression_results),
    )),
  })
}

/// The settings of a draft to preview the root-to-tip regression of, as the `clock` command computes it.
#[derive(Clone, Debug, JsonSchema, Deserialize)]
#[schemars(deny_unknown_fields)]
#[deser(deny_unknown_fields)]
pub struct ClockPreviewRequest {
  /// Command of the draft.
  pub command: AppCommand,
  /// Settings of the draft. The preview takes the settings that the `clock` command shares with the draft's command
  /// and that concern the input data or the analysis, such as the input files, the dates, the rooting, and the clock
  /// filter, and ignores the rest.
  pub config: SparseConfig,
  /// Whether to compute the preview also for a tree with more tips than `tip_limit`.
  #[schemars(default)]
  #[deser(default)]
  pub beyond_limit: bool,
}

/// The root-to-tip regression of a draft, computed as a `clock` run computes it, without a run.
#[derive(Clone, Debug, PartialEq, JsonSchema, Serialize)]
#[deser(skip_serializing_optionals)]
#[schemars(transform = skip_serializing_optionals)]
pub struct ClockPreview {
  /// Number of tips of the tree.
  pub tips: usize,
  /// Largest tree, in tips, whose preview the back end computes without `beyond_limit`.
  pub tip_limit: usize,
  /// Estimates and root-to-tip regression; absent when the tree has more tips than `tip_limit` and the request does
  /// not set `beyond_limit`.
  pub clock: Option<ClockResults>,
}

fn clock_config(command: AppCommand, config: &SparseConfig) -> Result<SparseConfig, Report> {
  let shared: BTreeSet<String> = command_settings(AppCommand::Clock)?
    .settings
    .into_iter()
    .filter(|spec| matches!(spec.purpose, SettingPurpose::Data | SettingPurpose::Analysis))
    .map(|spec| spec.key)
    .collect();
  let draft_keys: BTreeSet<String> = command_settings(command)?
    .settings
    .into_iter()
    .map(|spec| spec.key)
    .collect();
  let settings: Map<_, _> = config
    .0
    .iter()
    .filter(|(key, _)| shared.contains(*key) && draft_keys.contains(*key))
    .map(|(key, value)| (key.clone(), value.clone()))
    .collect();
  Ok(SparseConfig(settings))
}
