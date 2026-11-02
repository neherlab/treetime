use crate::commands::shared::argument_error::argument_error;
use crate::commands::shared::dependent_settings::{check_applies, setting_flag};
#[cfg(feature = "clap")]
use clap::ValueHint;
use deser::{Deserialize, Serialize};
use eyre::Report;
use itertools::{Either, Itertools};
use schemars::JsonSchema;
use smart_default::SmartDefault;
use std::collections::BTreeSet;
use std::path::PathBuf;
use treetime::clock::prune_tips::PruneTips;
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
  /// after the filter: timetree before the time inference, and clock before the final clock fit. Every output leaves the
  /// removed tips out, and a run warning names them. The tips that --keep-tips and --keep-tips-file name stay in the
  /// tree, flagged. Requires the clock filter (--clock-filter above 0).
  ///
  /// augur refine with a clock filter (--clock-filter-iqd) removes the flagged tips and also the tips without a usable
  /// date. To remove the same tips, add --prune-undated-tips.
  #[cfg_attr(feature = "clap", clap(long, help_heading = "Molecular clock"))]
  pub prune_filtered_tips: bool,

  /// Remove the tips without a usable date from the tree.
  ///
  /// A tip has no usable date when the metadata gives no date for it, or a value that is not a date, a numeric date, or
  /// a date range. Such a tip stays in the tree by default: timetree infers its date from the tree, and the clock fit
  /// leaves it out. With this flag, the run removes these tips: timetree before the time inference, and clock before
  /// the final clock fit. Every output leaves the removed tips out, and a run warning names them. The tips that
  /// --keep-tips and --keep-tips-file name stay in the tree.
  ///
  /// Together with --prune-filtered-tips, the run removes the same tips as augur refine with a clock filter
  /// (--clock-filter-iqd).
  #[cfg_attr(feature = "clap", clap(long, help_heading = "Molecular clock"))]
  pub prune_undated_tips: bool,

  /// Names of tips that stay in the tree when --prune-filtered-tips or --prune-undated-tips removes tips
  ///
  /// Accepts several names, separated by commas or spaces. A kept flagged tip stays flagged: the clock fit leaves it
  /// out, and timetree dates it from the tree, as it dates a kept tip without a usable date. Names that match no tip are
  /// ignored.
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

  /// File with names of tips that stay in the tree when --prune-filtered-tips or --prune-undated-tips removes tips, one
  /// name per line
  ///
  /// Surrounding whitespace is trimmed and empty lines are ignored. Names that match no tip are ignored, so the file
  /// can list samples of several datasets. A path that ends in `.gz`, `.bz2`, `.xz`, or `.zst` reads a compressed file.
  /// Use `-` to read uncompressed input from standard input (stdin).
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
pub struct TipPruning {
  filtered: bool,
  undated: bool,
  keep_tips: Vec<String>,
  keep_tips_file: Option<PathBuf>,
}

impl TipPruning {
  pub fn new(
    clock_filter: f64,
    prune_filtered_tips: bool,
    prune_undated_tips: bool,
    keep_tips: Vec<String>,
    keep_tips_file: Option<PathBuf>,
  ) -> Result<Self, Report> {
    let pruning = prune_filtered_tips || prune_undated_tips;
    for (setting, given) in [
      ("keep_tips", !keep_tips.is_empty()),
      ("keep_tips_file", keep_tips_file.is_some()),
    ] {
      if given && !pruning {
        return Err(argument_error(
          &[setting, "prune_filtered_tips", "prune_undated_tips"],
          format!(
            "{} applies only with --prune-filtered-tips or --prune-undated-tips: it names tips that the run keeps",
            setting_flag(setting)
          ),
        ));
      }
    }
    check_applies(
      "prune_filtered_tips",
      prune_filtered_tips,
      "clock_filter",
      clock_filter > 0.0,
      || o!("removes the tips that the clock filter flags; set --clock-filter above 0"),
    )?;
    Ok(Self {
      filtered: prune_filtered_tips,
      undated: prune_undated_tips,
      keep_tips,
      keep_tips_file,
    })
  }

  pub fn resolve(&self, graph: &Graph, names: &NodeNames, log: &dyn LogSink) -> Result<PruneTips, Report> {
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
        "Ignoring {} of the keep list (--keep-tips, --keep-tips-file) without a matching tip in the tree",
        counted(unmatched.len(), "name")
      );
    }
    Ok(PruneTips {
      filtered: self.filtered,
      undated: self.undated,
      keep,
    })
  }
}
