use crate::commands::shared::dependent_settings::check_applies;
#[cfg(feature = "clap")]
use clap::ValueHint;
use deser::{Deserialize, Serialize};
use eyre::Report;
use itertools::{Either, Itertools};
use schemars::JsonSchema;
use smart_default::SmartDefault;
use std::collections::BTreeSet;
use std::path::PathBuf;
use treetime::clock::prune_filtered::PruneFilteredTips;
use treetime::progress::LogSink;
use treetime::progress_info;
use treetime_graph::graph::Graph;
use treetime_graph::node_names::NodeNames;
use treetime_graph::pair_by_name::key_by_name;
use treetime_io::name_list::name_list_read_file;
use treetime_schema::schema_defaults;
use treetime_utils::fmt::string::counted;
use treetime_utils::o;

#[derive(Debug, Clone, SmartDefault, JsonSchema, Serialize, Deserialize)]
#[schemars(default, deny_unknown_fields)]
#[schemars(transform = schema_defaults::<Self>)]
#[deser(default, deny_unknown_fields)]
#[cfg_attr(feature = "clap", derive(clap::Args))]
pub struct ClockFilterArgs {
  /// Ignore tips that do not follow a loose clock: the number of interquartile ranges from the root-to-tip regression
  /// beyond which a tip is an outlier.
  ///
  /// 0 switches the filter off.
  #[default = 4.0]
  #[cfg_attr(
    feature = "clap",
    clap(long, default_value_t = ClockFilterArgs::default().clock_filter, help_heading = "Molecular clock")
  )]
  pub clock_filter: f64,

  /// Allow negative clock rates.
  ///
  /// By default, rates are forced to be positive. For trees with little temporal signal, removing this restriction
  /// gives essentially mid-point rooting.
  #[cfg_attr(
    feature = "clap",
    clap(long, hide_short_help = true, help_heading = "Molecular clock")
  )]
  pub allow_negative_rate: bool,

  /// Remove the tips that the clock filter flags from the tree.
  ///
  /// The clock filter flags a tip and leaves it in the tree by default. With this flag, the run removes the flagged tips
  /// after the filter: timetree before the time inference, as augur refine does, and clock before the final clock fit.
  /// Every output leaves the removed tips out. The tips that --keep-tips and --keep-tips-file name stay in the tree,
  /// flagged. Requires the clock filter (--clock-filter above 0).
  #[cfg_attr(feature = "clap", clap(long, help_heading = "Molecular clock"))]
  pub prune_filtered_tips: bool,

  /// Names of tips that stay in the tree when --prune-filtered-tips removes the flagged tips
  ///
  /// Accepts several names, separated by commas. A kept tip stays flagged: the clock fit leaves it out, and the time
  /// inference dates it from the tree. Names that match no tip are ignored.
  #[cfg_attr(
    feature = "clap",
    clap(
      long,
      num_args = 1..,
      value_delimiter = ',',
      value_name = "NAME",
      hide_short_help = true,
      help_heading = "Molecular clock"
    )
  )]
  pub keep_tips: Vec<String>,

  /// File with names of tips that stay in the tree when --prune-filtered-tips removes the flagged tips, one name per
  /// line
  ///
  /// Surrounding whitespace is trimmed and empty lines are ignored. Names that match no tip are ignored, so the file
  /// can list samples of several datasets. Use '-' to read from standard input (stdin).
  #[cfg_attr(
    feature = "clap",
    clap(
      long,
      value_hint = ValueHint::FilePath,
      value_name = "FILE",
      hide_short_help = true,
      help_heading = "Molecular clock"
    )
  )]
  #[schemars(extend("x-path" = "input"))]
  pub keep_tips_file: Option<PathBuf>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FilteredTipPruning {
  keep_tips: Vec<String>,
  keep_tips_file: Option<PathBuf>,
}

impl FilteredTipPruning {
  pub fn new(
    clock_filter: f64,
    prune_filtered_tips: bool,
    keep_tips: Vec<String>,
    keep_tips_file: Option<PathBuf>,
  ) -> Result<Option<Self>, Report> {
    let keep_rule = || o!("applies only with --prune-filtered-tips: it names flagged tips that the run keeps");
    check_applies(
      "keep_tips",
      !keep_tips.is_empty(),
      "prune_filtered_tips",
      prune_filtered_tips,
      keep_rule,
    )?;
    check_applies(
      "keep_tips_file",
      keep_tips_file.is_some(),
      "prune_filtered_tips",
      prune_filtered_tips,
      keep_rule,
    )?;
    check_applies(
      "prune_filtered_tips",
      prune_filtered_tips,
      "clock_filter",
      clock_filter > 0.0,
      || o!("removes the tips that the clock filter flags; set --clock-filter above 0"),
    )?;
    Ok(prune_filtered_tips.then_some(Self {
      keep_tips,
      keep_tips_file,
    }))
  }

  pub fn resolve(&self, graph: &Graph, names: &NodeNames, log: &dyn LogSink) -> Result<PruneFilteredTips, Report> {
    let from_file = self.keep_tips_file.as_deref().map(name_list_read_file).transpose()?;
    let requested: BTreeSet<String> = self
      .keep_tips
      .iter()
      .cloned()
      .chain(from_file.into_iter().flatten())
      .collect();
    let leaf_of = key_by_name(graph.get_leaves().map(|leaf| leaf.key()), names);
    let (keep, unmatched): (BTreeSet<_>, Vec<_>) =
      requested.iter().partition_map(|name| match leaf_of.get(name.as_str()) {
        Some(&leaf) => Either::Left(leaf),
        None => Either::Right(name),
      });
    if !unmatched.is_empty() {
      progress_info!(
        log,
        "Ignoring {} of the keep list of --prune-filtered-tips without a matching tip in the tree",
        counted(unmatched.len(), "name")
      );
    }
    Ok(PruneFilteredTips { keep })
  }
}
