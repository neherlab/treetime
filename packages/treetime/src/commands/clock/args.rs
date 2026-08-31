use crate::ancestral::params::MethodAncestral;
use crate::clock::clock_regression::ClockParams;
use crate::clock::find_best_root::params::{BrentParams, GoldenSectionParams, GridSearchParams, OptimizationMethod};
use crate::commands::shared::alignment::AlignmentArgs;
use crate::commands::shared::config::ConfigArgs;
use crate::commands::shared::metadata::{DateColumnArgs, MetadataIdArgs};
use crate::commands::shared::model::ModelArgs;
use crate::commands::shared::output::{ClockOutputSelection, OutputCoreArgs, TopologyOrderArgs};
use crate::commands::shared::required::missing_required_args;
use crate::commands::shared::reroot::RerootArgs;
use crate::optimize::params::BranchLengthMode;
#[cfg(feature = "clap")]
use clap::ValueHint;
use eyre::Report;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use smart_default::SmartDefault;
use std::fmt::Debug;
use std::path::{Path, PathBuf};

#[derive(Debug, Clone, SmartDefault, Serialize, Deserialize, JsonSchema)]
#[serde(default, deny_unknown_fields)]
#[cfg_attr(feature = "clap", derive(clap::Parser))]
pub struct TreetimeClockArgsRaw {
  #[cfg_attr(feature = "clap", clap(flatten))]
  #[serde(skip)]
  pub config_args: ConfigArgs,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[serde(flatten)]
  pub alignment: AlignmentArgs,

  /// Name of file containing the tree in newick, nexus, or phylip format.
  ///
  /// If none is provided, treetime will attempt to build a tree from the alignment using fasttree, iqtree, or raxml (assuming they are installed)
  #[cfg_attr(feature = "clap", clap(long, short = 't'))]
  #[cfg_attr(feature = "clap", clap(value_hint = ValueHint::FilePath))]
  pub tree: Option<PathBuf>,

  /// Only for vcf input: fasta file of the sequence the VCF was mapped to.
  #[cfg_attr(feature = "clap", clap(long, short = 'r'))]
  #[cfg_attr(feature = "clap", clap(value_hint = ValueHint::FilePath))]
  pub vcf_reference: Option<PathBuf>,

  /// CSV/TSV file with metadata including sampling dates
  #[cfg_attr(feature = "clap", clap(long = "metadata", visible_alias = "dates", short = 'd'))]
  #[cfg_attr(feature = "clap", clap(value_hint = ValueHint::FilePath))]
  pub metadata: Option<PathBuf>,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[serde(flatten)]
  pub metadata_id: MetadataIdArgs,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[serde(flatten)]
  pub date_column: DateColumnArgs,

  /// Length of the sequence, used to calculate expected variation in branch length. Not required if alignment is provided.
  #[cfg_attr(feature = "clap", clap(long))]
  pub sequence_length: Option<usize>,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[serde(flatten)]
  pub model_args: ModelArgs,

  /// If set to 'input', the provided branch length will be used without modification. Note that branch lengths optimized by treetime are only accurate at short evolutionary distances.
  #[cfg_attr(feature = "clap", clap(long, value_enum, default_value_t = BranchLengthMode::default()))]
  pub branch_length_mode: BranchLengthMode,

  /// Method used for reconstructing ancestral sequences
  #[cfg_attr(feature = "clap", clap(long, value_enum, default_value_t = MethodAncestral::default()))]
  pub method_anc: MethodAncestral,

  /// ignore tips that don't follow a loose clock, 'clock-filter=number of interquartile ranges from regression'. Default=3.0, set to 0 to switch off.
  #[cfg_attr(feature = "clap", clap(long, default_value = "3.0"))]
  #[default = 3.0]
  pub clock_filter: f64,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[serde(flatten)]
  pub reroot: RerootArgs,

  /// don't reroot the tree. Otherwise, reroot to minimize the residual of the regression of
  /// root-to-tip distance and sampling time
  #[cfg_attr(feature = "clap", clap(long, conflicts_with_all = ["reroot", "reroot_tips"]))]
  pub keep_root: bool,

  #[cfg_attr(feature = "clap", clap(long))]
  pub prune_short: bool,

  /// excess variance associated with terminal nodes accounting for overdispersion of the molecular
  /// clock
  #[cfg_attr(feature = "clap", clap(long))]
  pub tip_slack: Option<f64>,

  /// Account for covariation when estimating rates or rerooting using root-to-tip regression
  #[cfg_attr(feature = "clap", clap(long))]
  pub covariation: bool,

  /// By default, rates are forced to be positive. For trees with little temporal signal it is advisable to remove this restriction to achieve essentially mid-point rooting.
  #[cfg_attr(feature = "clap", clap(long))]
  pub allow_negative_rate: bool,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[serde(flatten)]
  pub output: OutputCoreArgs,

  /// Path to output clock model JSON.
  ///
  /// Takes precedence over paths configured with `--output-all` and `--output-selection`.
  #[cfg_attr(feature = "clap", clap(long, value_hint = ValueHint::FilePath, help_heading = "Output"))]
  pub output_clock_model: Option<PathBuf>,

  /// Path to output clock regression CSV.
  ///
  /// Takes precedence over paths configured with `--output-all` and `--output-selection`.
  #[cfg_attr(feature = "clap", clap(long, value_hint = ValueHint::FilePath, help_heading = "Output"))]
  pub output_clock_csv: Option<PathBuf>,

  /// Comma-separated list of outputs to produce with `--output-all`.
  ///
  /// Restricts which outputs `--output-all` writes. Special value `all` expands to every output
  /// available for this command. Requires `--output-all`. Per-file flags are always honored
  /// regardless of this selection.
  #[cfg_attr(
    feature = "clap",
    clap(long, value_delimiter = ',', requires = "output_all", help_heading = "Output")
  )]
  pub output_selection: Vec<ClockOutputSelection>,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[serde(flatten)]
  pub topology_order: TopologyOrderArgs,

  /// Random seed
  #[cfg_attr(feature = "clap", clap(long, visible_alias = "rng-seed"))]
  pub seed: Option<u64>,

  /// Method for clock filter outlier detection (not yet implemented)
  #[cfg_attr(feature = "clap", clap(long, hide = true))]
  pub clock_filter_method: Option<String>,

  /// Filename to save root-to-tip regression plot (not yet implemented)
  #[cfg_attr(feature = "clap", clap(long, hide = true))]
  #[cfg_attr(feature = "clap", clap(value_hint = ValueHint::FilePath))]
  pub plot_rtt: Option<PathBuf>,

  /// Prune clock outlier tips from the tree (not yet implemented)
  #[cfg_attr(feature = "clap", clap(long, hide = true))]
  pub prune_outliers: bool,

  /// Branch split optimization parameters
  #[cfg_attr(feature = "clap", clap(flatten, next_help_heading = "Branch split optimization"))]
  pub branch_split: BranchSplitArgs,

  /// Clock regression model parameters
  #[cfg_attr(feature = "clap", clap(flatten, next_help_heading = "Clock regression"))]
  pub clock_regression: ClockRegressionArgs,
}

/// Clock arguments with required inputs proven present.
///
/// Produced from [`TreetimeClockArgsRaw`] by [`TryFrom`] once the `--config` overlay has run, so the
/// run code reads `metadata` without an `Option`.
#[derive(Debug, Clone)]
pub struct TreetimeClockArgs {
  pub alignment: AlignmentArgs,
  pub tree: Option<PathBuf>,
  pub vcf_reference: Option<PathBuf>,
  pub metadata: PathBuf,
  pub metadata_id: MetadataIdArgs,
  pub date_column: DateColumnArgs,
  pub sequence_length: Option<usize>,
  pub model_args: ModelArgs,
  pub branch_length_mode: BranchLengthMode,
  pub method_anc: MethodAncestral,
  pub clock_filter: f64,
  pub reroot: RerootArgs,
  pub keep_root: bool,
  pub prune_short: bool,
  pub tip_slack: Option<f64>,
  pub covariation: bool,
  pub allow_negative_rate: bool,
  pub output: OutputCoreArgs,
  pub output_clock_model: Option<PathBuf>,
  pub output_clock_csv: Option<PathBuf>,
  pub output_selection: Vec<ClockOutputSelection>,
  pub topology_order: TopologyOrderArgs,
  pub seed: Option<u64>,
  pub clock_filter_method: Option<String>,
  pub plot_rtt: Option<PathBuf>,
  pub prune_outliers: bool,
  pub branch_split: BranchSplitArgs,
  pub clock_regression: ClockRegressionArgs,
}

impl TreetimeClockArgs {
  /// Metadata (dates) path.
  pub fn metadata(&self) -> &Path {
    &self.metadata
  }
}

impl TryFrom<TreetimeClockArgsRaw> for TreetimeClockArgs {
  type Error = Report;

  fn try_from(raw: TreetimeClockArgsRaw) -> Result<Self, Report> {
    let metadata = raw
      .metadata
      .ok_or_else(|| missing_required_args::<TreetimeClockArgsRaw>(&["metadata"]))?;
    Ok(Self {
      alignment: raw.alignment,
      tree: raw.tree,
      vcf_reference: raw.vcf_reference,
      metadata,
      metadata_id: raw.metadata_id,
      date_column: raw.date_column,
      sequence_length: raw.sequence_length,
      model_args: raw.model_args,
      branch_length_mode: raw.branch_length_mode,
      method_anc: raw.method_anc,
      clock_filter: raw.clock_filter,
      reroot: raw.reroot,
      keep_root: raw.keep_root,
      prune_short: raw.prune_short,
      tip_slack: raw.tip_slack,
      covariation: raw.covariation,
      allow_negative_rate: raw.allow_negative_rate,
      output: raw.output,
      output_clock_model: raw.output_clock_model,
      output_clock_csv: raw.output_clock_csv,
      output_selection: raw.output_selection,
      topology_order: raw.topology_order,
      seed: raw.seed,
      clock_filter_method: raw.clock_filter_method,
      plot_rtt: raw.plot_rtt,
      prune_outliers: raw.prune_outliers,
      branch_split: raw.branch_split,
      clock_regression: raw.clock_regression,
    })
  }
}

/// Branch split optimization parameters
#[derive(Debug, Clone, SmartDefault, Serialize, Deserialize, JsonSchema)]
#[serde(default, deny_unknown_fields)]
#[cfg_attr(feature = "clap", derive(clap::Args))]
pub struct BranchSplitArgs {
  /// Optimization method to use for finding the best root position
  #[cfg_attr(feature = "clap", clap(long = "branch-split-method", value_enum, default_value_t = OptimizationMethod::default()))]
  #[default(OptimizationMethod::default())]
  pub method: OptimizationMethod,

  /// Grid search parameters
  #[cfg_attr(feature = "clap", clap(flatten))]
  #[serde(flatten)]
  pub grid_params: GridSearchParams,

  /// Brent's method parameters
  #[cfg_attr(feature = "clap", clap(flatten))]
  #[serde(flatten)]
  pub brent_params: BrentParams,

  /// Golden section search parameters
  #[cfg_attr(feature = "clap", clap(flatten))]
  #[serde(flatten)]
  pub golden_params: GoldenSectionParams,
}

/// Clock regression model parameters
#[derive(Debug, Clone, SmartDefault, Serialize, Deserialize, JsonSchema)]
#[serde(default, deny_unknown_fields)]
#[cfg_attr(feature = "clap", derive(clap::Args))]
pub struct ClockRegressionArgs {
  /// Clock regression model parameters
  #[cfg_attr(feature = "clap", clap(flatten))]
  #[serde(flatten)]
  pub clock_params: ClockParams,
}
