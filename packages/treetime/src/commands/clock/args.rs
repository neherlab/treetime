use crate::ancestral::params::MethodAncestral;
use crate::clock::clock_regression::ClockParams;
use crate::clock::find_best_root::params::{BrentParams, GoldenSectionParams, GridSearchParams, OptimizationMethod};
use crate::commands::shared::alignment::AlignmentArgs;
use crate::commands::shared::metadata::{DateColumnArgs, MetadataIdArgs};
use crate::commands::shared::model::ModelArgs;
use crate::commands::shared::output::OutputCoreArgs;
use crate::commands::shared::reroot::RerootArgs;
use crate::optimize::params::BranchLengthMode;
#[cfg(feature = "clap")]
use clap::ValueHint;
use serde::{Deserialize, Serialize};
use smart_default::SmartDefault;
use std::fmt::Debug;
use std::path::PathBuf;

#[derive(Debug, SmartDefault, Serialize, Deserialize)]
#[serde(default)]
#[cfg_attr(feature = "clap", derive(clap::Parser))]
pub struct TreetimeClockArgs {
  #[cfg_attr(feature = "clap", clap(flatten))]
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
  pub metadata: PathBuf,

  #[cfg_attr(feature = "clap", clap(flatten))]
  pub metadata_id: MetadataIdArgs,

  #[cfg_attr(feature = "clap", clap(flatten))]
  pub date_column: DateColumnArgs,

  /// Length of the sequence, used to calculate expected variation in branch length. Not required if alignment is provided.
  #[cfg_attr(feature = "clap", clap(long))]
  pub sequence_length: Option<usize>,

  #[cfg_attr(feature = "clap", clap(flatten))]
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

/// Branch split optimization parameters
#[derive(Debug, Clone, SmartDefault, Serialize, Deserialize)]
#[serde(default)]
#[cfg_attr(feature = "clap", derive(clap::Args))]
pub struct BranchSplitArgs {
  /// Optimization method to use for finding the best root position
  #[cfg_attr(feature = "clap", clap(long = "branch-split-method", value_enum, default_value_t = OptimizationMethod::default()))]
  #[default(OptimizationMethod::default())]
  pub method: OptimizationMethod,

  /// Grid search parameters
  #[cfg_attr(feature = "clap", clap(flatten))]
  pub grid_params: GridSearchParams,

  /// Brent's method parameters
  #[cfg_attr(feature = "clap", clap(flatten))]
  pub brent_params: BrentParams,

  /// Golden section search parameters
  #[cfg_attr(feature = "clap", clap(flatten))]
  pub golden_params: GoldenSectionParams,
}

/// Clock regression model parameters
#[derive(Debug, Clone, SmartDefault, Serialize, Deserialize)]
#[serde(default)]
#[cfg_attr(feature = "clap", derive(clap::Args))]
pub struct ClockRegressionArgs {
  /// Clock regression model parameters
  #[cfg_attr(feature = "clap", clap(flatten))]
  pub clock_params: ClockParams,
}
