use crate::command::{AppCommandExt, config_from_request};
use crate::config::catalog::command_settings;
use crate::results::auspice::{AuspiceDocument, display_auspice};
use crate::results::clock::{ClockResults, clock_results_of_tips};
use crate::runs::errors::invalid;
use app_commands::command::AppCommand;
use app_commands::command_config::CommandConfig;
use app_commands::commands::clock::run::{ClockRegression, clock_regression, input_tree_view};
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
use treetime_utils::error::report_to_string;
use treetime_utils::make_internal_error;

pub const CLOCK_PREVIEW_TIP_LIMIT: usize = 5100;

pub fn clock_preview(request: &ClockPreviewRequest) -> Result<ClockPreview, Report> {
  let config = clock_config(request.command, &request.config)?;
  let CommandConfig::Clock { config } = config_from_request(AppCommand::Clock, &config)? else {
    return make_internal_error!("a clock configuration reads as the configuration of another command");
  };
  let Some(tree) = &config.tree.tree else {
    return Err(invalid("a clock preview needs a tree"));
  };
  let parsed = tree_read_file(tree, config.newick_dialect.newick_dialect.dialect())?;
  let tips = parsed.graph.num_leaves();
  let preview = |clock, tree, problem| ClockPreview {
    tips,
    tip_limit: CLOCK_PREVIEW_TIP_LIMIT,
    clock,
    tree,
    problem,
  };
  if tips > CLOCK_PREVIEW_TIP_LIMIT && !request.beyond_limit {
    return Ok(preview(None, None, None));
  }
  let input_tree = || -> Result<PreviewTree, Report> {
    Ok(PreviewTree {
      root: PreviewRoot::Input,
      auspice: display_auspice(input_tree_view(&parsed)?)?,
    })
  };
  if !request.command.uses_dates() || config.metadata.metadata.is_none() {
    return Ok(preview(None, Some(input_tree()?), None));
  }
  match clock_regression(*config, &NoopCancel, &NoopProgress, &NoopProgress) {
    Ok(ClockRegression {
      clock_model,
      regression_results,
      tip_names,
      tree,
    }) => {
      let tip_names: BTreeSet<&str> = tip_names.iter().map(String::as_str).collect();
      let clock = clock_results_of_tips(Some(&tip_names), Some(&clock_model), Some(&regression_results));
      let tree = PreviewTree {
        root: PreviewRoot::Regression,
        auspice: display_auspice(tree)?,
      };
      Ok(preview(Some(clock), Some(tree), None))
    },
    Err(report) => Ok(preview(None, Some(input_tree()?), Some(report_to_string(&report)))),
  }
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

/// The root-to-tip regression of a draft, computed as a `clock` run computes it, and its tree, without a run.
#[derive(Clone, Debug, JsonSchema, Serialize)]
#[deser(skip_serializing_optionals)]
#[schemars(transform = skip_serializing_optionals)]
pub struct ClockPreview {
  /// Number of tips of the tree.
  pub tips: usize,
  /// Largest tree, in tips, whose preview the back end computes without `beyond_limit`.
  pub tip_limit: usize,
  /// Estimates and root-to-tip regression; absent when the tree has more tips than `tip_limit` and the request does
  /// not set `beyond_limit`, when the draft has no metadata, when the draft's command uses no dates, and when the
  /// regression fails.
  pub clock: Option<ClockResults>,
  /// The tree in divergence units, without dates; absent when the tree has more tips than `tip_limit` and the request
  /// does not set `beyond_limit`.
  pub tree: Option<PreviewTree>,
  /// Why the root-to-tip regression failed, as the CLI prints the error; the tree then has the root of the input tree
  /// file.
  pub problem: Option<String>,
}

/// The tree of a draft, as the preview draws it.
#[derive(Clone, Debug, JsonSchema, Serialize)]
pub struct PreviewTree {
  /// Where the root of the tree comes from.
  pub root: PreviewRoot,
  /// The tree as an Auspice document, with the divergence of every node and no dates.
  pub auspice: AuspiceDocument,
}

/// Where the root of a preview tree comes from.
#[derive(Clone, Copy, Debug, PartialEq, Eq, JsonSchema, Serialize)]
#[schemars(rename_all = "kebab-case")]
#[deser(rename_all = "kebab-case")]
pub enum PreviewRoot {
  /// The root that the root-to-tip regression places, as in `clock`.
  Regression,
  /// The root of the input tree file, because the preview has no regression: the draft has no dates, or the regression
  /// failed.
  Input,
}

fn clock_config(command: AppCommand, config: &SparseConfig) -> Result<SparseConfig, Report> {
  let read: BTreeSet<String> = command_settings(command)?.preview_settings.into_iter().collect();
  let settings: Map<_, _> = config
    .0
    .iter()
    .filter(|(key, _)| read.contains(*key))
    .map(|(key, value)| (key.clone(), value.clone()))
    .collect();
  Ok(SparseConfig(settings))
}
