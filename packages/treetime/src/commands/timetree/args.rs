use crate::ancestral::params::MethodAncestral;
use crate::commands::shared::alignment::AlignmentArgs;
use crate::commands::shared::alphabet::AlphabetArgs;
use crate::commands::shared::gap_fill::GapFillArgs;
use crate::commands::shared::metadata::{DateColumnArgs, MetadataIdArgs};
use crate::commands::shared::model::ModelArgs;
use crate::commands::shared::output::{DivergenceUnits, OutputCoreArgs, TimetreeOutputSelection, TopologyOrderArgs};
use crate::commands::shared::reroot::RerootArgs;
use crate::optimize::params::BranchLengthMode;
#[cfg(feature = "clap")]
use clap::ValueHint;
use serde::{Deserialize, Serialize};
use smart_default::SmartDefault;
use std::fmt::Debug;
use std::path::PathBuf;

#[cfg(feature = "clap")]
fn parse_n_skyline(s: &str) -> Result<usize, String> {
  let n: usize = s.parse().map_err(|_err| format!("'{s}' is not a valid number"))?;
  if n < 2 {
    return Err("n-skyline must be at least 2".to_owned());
  }
  Ok(n)
}

pub use crate::timetree::params::TimeMarginalMode;

#[derive(Debug, SmartDefault, Serialize, Deserialize)]
#[serde(default)]
#[cfg_attr(feature = "clap", derive(clap::Parser))]
pub struct TreetimeTimetreeArgs {
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
  pub metadata: Option<PathBuf>,

  #[cfg_attr(feature = "clap", clap(flatten))]
  pub metadata_id: MetadataIdArgs,

  #[cfg_attr(feature = "clap", clap(flatten))]
  pub date_column_args: DateColumnArgs,

  /// Length of the sequence, used to calculate expected variation in branch length. Not required if alignment is provided.
  #[cfg_attr(feature = "clap", clap(long))]
  pub sequence_length: Option<usize>,

  /// If specified, the rate of the molecular clock won't be optimized.
  #[cfg_attr(feature = "clap", clap(long))]
  pub clock_rate: Option<f64>,

  /// Standard deviation of the provided clock rate estimate
  #[cfg_attr(feature = "clap", clap(long))]
  pub clock_std_dev: Option<f64>,

  /// If set to 'input', the provided branch length will be used without modification. Branch lengths optimized by treetime are only accurate at short evolutionary distances.
  #[cfg_attr(feature = "clap", clap(long, value_enum, default_value_t = BranchLengthMode::default()))]
  pub branch_length_mode: BranchLengthMode,

  /// Control when marginal time distributions are used for output.
  ///
  /// All modes use marginal inference during optimization. The mode controls whether
  /// confidence intervals are extracted from the resulting distributions:
  ///
  /// - `never`: no confidence interval output (default)
  /// - `always`: write confidence intervals from distributions computed during optimization
  /// - `only-final`: run one extra inference pass after optimization, then write confidence intervals
  #[cfg_attr(feature = "clap", clap(long, value_enum, default_value_t = TimeMarginalMode::default()))]
  pub time_marginal: TimeMarginalMode,

  /// Add rate-uncertainty to confidence intervals.
  ///
  /// `--time-marginal=always` and `only-final` already write mutation-stochasticity CIs.
  /// This flag adds rate-uncertainty CIs (re-runs inference at rate +/- sigma), combined
  /// via quadrature sum. Requires `--covariation` or `--clock-std-dev`.
  /// When set with `--time-marginal=never` (default), automatically promotes to `only-final`.
  #[cfg_attr(feature = "clap", clap(long))]
  pub confidence: bool,

  /// Don't resolve polytomies using temporal information.
  #[cfg_attr(feature = "clap", clap(long))]
  pub keep_polytomies: bool,

  /// Resolve polytomies using temporal information
  #[cfg_attr(feature = "clap", clap(long))]
  pub resolve_polytomies: bool,

  /// use an autocorrelated molecular clock. Strength of the gaussian priors on branch specific rate
  /// deviation and the coupling of parent and offspring rates can be specified e.g. as --relax 1.0
  /// 0.5. Values around 1.0 correspond to weak priors, larger values constrain rate deviations more
  /// strongly. Coupling 0 (--relax 1.0 0) corresponds to an un-correlated clock.
  #[cfg_attr(feature = "clap", clap(long, num_args = 2, value_names = ["SLACK", "COUPLING"]))]
  pub relax: Vec<f64>,

  /// maximal number of iterations the inference cycle is run. For polytomy resolution and
  /// coalescence models max_iter should be at least 2
  #[default = 2]
  #[cfg_attr(feature = "clap", clap(long, default_value_t = TreetimeTimetreeArgs::default().max_iter))]
  pub max_iter: usize,

  /// Coalescent time scale in years.
  ///
  /// Sensible values are on the order of the time from the root to the tips and are given in units of time.
  #[cfg_attr(feature = "clap", clap(long))]
  pub coalescent: Option<f64>,

  /// Optimize coalescent time scale Tc to maximize coalescent likelihood.
  ///
  /// When set, TreeTime will find the optimal Tc using Brent's method. This is similar
  /// to Python v0's `--coalescent=opt`, but uses a closed-form solution.
  #[cfg_attr(
    feature = "clap",
    clap(long, conflicts_with = "coalescent", conflicts_with = "coalescent_skyline")
  )]
  pub coalescent_opt: bool,

  /// Use skyline coalescent model instead of constant Tc.
  ///
  /// Estimates a piecewise linear coalescent rate history. Requires --n-skyline to specify
  /// the number of grid points.
  #[cfg_attr(feature = "clap", clap(long))]
  #[cfg_attr(
    feature = "clap",
    clap(conflicts_with = "coalescent", conflicts_with = "coalescent_opt")
  )]
  pub coalescent_skyline: bool,

  /// Number of grid points in skyline coalescent model.
  ///
  /// Only used when --coalescent-skyline is set. Defines how many piecewise linear segments
  /// are used to model Tc(t) over time. Must be at least 2.
  #[default = 10]
  #[cfg_attr(feature = "clap", clap(long, default_value_t = TreetimeTimetreeArgs::default().n_skyline))]
  #[cfg_attr(feature = "clap", clap(value_parser = parse_n_skyline))]
  pub n_skyline: usize,

  /// Smoothing stiffness for the skyline coalescent.
  ///
  /// Penalizes changes in the inverse coalescent time scale 1/Tc between adjacent
  /// skyline segments: the objective adds `stiffness * Σ (1/Tc_{i+1} - 1/Tc_i)^2`.
  /// Because 1/Tc has units of 1/time, the stiffness has units of time^2. Larger
  /// values enforce a smoother Tc(t). Only used when --coalescent-skyline is set.
  #[default = 2.0]
  #[cfg_attr(feature = "clap", clap(long, default_value_t = TreetimeTimetreeArgs::default().skyline_stiffness))]
  pub skyline_stiffness: f64,

  /// add posterior LH to coalescent model: use the posterior probability distributions of
  /// divergence times for estimating the number of branches when calculating the coalescent
  /// mergerrate or use inferred time before present (default).
  #[cfg_attr(feature = "clap", clap(long, hide = true))]
  pub n_branches_posterior: Option<usize>,

  /// filename to save the plot to. Suffix will determine format (choices pdf, png, svg,
  /// default=pdf)
  #[cfg_attr(feature = "clap", clap(long, hide = true))]
  #[cfg_attr(feature = "clap", clap(value_hint = ValueHint::FilePath))]
  pub plot_tree: Option<PathBuf>,

  /// filename to save the plot to. Suffix will determine format (choices pdf, png, svg,
  /// default=pdf)
  #[cfg_attr(feature = "clap", clap(long, hide = true))]
  #[cfg_attr(feature = "clap", clap(value_hint = ValueHint::FilePath))]
  pub plot_rtt: Option<PathBuf>,

  /// add tip labels (default for small trees with <30 leaves)
  #[cfg_attr(feature = "clap", clap(long))]
  pub tip_labels: bool,

  /// don't show tip labels (default for trees with >=30 leaves)
  #[cfg_attr(feature = "clap", clap(long))]
  pub no_tip_labels: bool,

  /// ignore tips that don't follow a loose clock, 'clock-filter=number of inter-quartile ranges from
  /// regression'. Default=3.0, set to 0 to switch off.
  #[cfg_attr(feature = "clap", clap(long, default_value = "3.0"))]
  pub clock_filter: f64,

  /// Number of IQD (interquartile distance) for clock filter outlier detection
  #[cfg_attr(feature = "clap", clap(long))]
  pub n_iqd: Option<f64>,

  #[cfg_attr(feature = "clap", clap(flatten))]
  pub reroot: RerootArgs,

  /// don't reroot the tree. Otherwise, reroot to minimize the residual of the regression of
  /// root-to-tip distance and sampling time
  #[cfg_attr(feature = "clap", clap(long, conflicts_with_all = ["reroot", "reroot_tips"]))]
  pub keep_root: bool,

  /// By default, rates are forced to be positive. For trees with little temporal signal it is advisable to remove this restriction to achieve essentially mid-point rooting.
  #[cfg_attr(feature = "clap", clap(long))]
  pub allow_negative_rate: bool,

  /// excess variance associated with terminal nodes accounting for overdispersion of the molecular
  /// clock
  #[cfg_attr(feature = "clap", clap(long))]
  pub tip_slack: Option<f64>,

  /// Account for covariation when estimating rates or rerooting using root-to-tip regression
  #[cfg_attr(feature = "clap", clap(long))]
  pub covariation: bool,

  #[cfg_attr(feature = "clap", clap(flatten))]
  pub model_args: ModelArgs,

  /// Method used for reconstructing ancestral sequences
  #[cfg_attr(feature = "clap", clap(long, value_enum, default_value_t = MethodAncestral::default()))]
  pub method_anc: MethodAncestral,

  #[cfg_attr(feature = "clap", clap(flatten))]
  pub alphabet_args: AlphabetArgs,

  /// Use dense representation for sequences (store full probability distributions)
  #[cfg_attr(feature = "clap", clap(long))]
  pub dense: Option<bool>,

  #[cfg_attr(feature = "clap", clap(flatten))]
  pub gap_fill_args: GapFillArgs,

  /// Zero-based mutation indexing
  #[cfg_attr(feature = "clap", clap(long))]
  pub zero_based: bool,

  /// Overwrite ambiguous states on tips with the most likely inferred state
  #[cfg_attr(feature = "clap", clap(long))]
  pub reconstruct_tip_states: bool,

  /// Include transitions involving ambiguous states
  #[cfg_attr(feature = "clap", clap(long))]
  pub report_ambiguous: bool,

  /// Disable indel (insertion/deletion) contributions to branch-length
  /// optimization and branch-length distributions.
  ///
  /// When set, branch-length optimization uses substitution-only likelihood
  /// and timetree branch distributions exclude the Poisson indel term.
  /// Matches standard phylogenetic tools (RAxML, IQ-TREE, PhyML, BEAST)
  /// and enables v0 parity testing. Default: indels enabled.
  #[cfg_attr(feature = "clap", clap(long))]
  pub no_indels: bool,

  /// Units for divergence values in augur node data JSON and auspice output.
  ///
  /// `mutations-per-site` (default): branch divergence as substitutions per site.
  /// `mutations`: absolute count of reconstructed substitutions per branch,
  /// excluding ambiguous and gap positions. Requires ancestral reconstruction
  /// (incompatible with `--branch-length-mode=input`).
  #[cfg_attr(feature = "clap", clap(long, value_enum, default_value_t = DivergenceUnits::default()))]
  pub divergence_units: DivergenceUnits,

  /// Path to output augur-compatible node data JSON.
  ///
  /// Contains per-node dates, branch lengths, clock model parameters, confidence
  /// intervals, and divergence metrics. The output is compatible with augur
  /// export v2 --node-data for Nextstrain pipeline integration.
  ///
  /// Takes precedence over paths configured with `--output-all` and `--output-selection`.
  #[cfg_attr(feature = "clap", clap(long, value_hint = ValueHint::FilePath, help_heading = "Output"))]
  pub output_augur_node_data: Option<PathBuf>,

  /// Path to output GTR model JSON.
  ///
  /// Takes precedence over paths configured with `--output-all` and `--output-selection`.
  #[cfg_attr(feature = "clap", clap(long, value_hint = ValueHint::FilePath, help_heading = "Output"))]
  pub output_gtr: Option<PathBuf>,

  /// Path to output clock model JSON.
  ///
  /// Takes precedence over paths configured with `--output-all` and `--output-selection`.
  #[cfg_attr(feature = "clap", clap(long, value_hint = ValueHint::FilePath, help_heading = "Output"))]
  pub output_clock_model: Option<PathBuf>,

  /// Path to output date-confidence-interval TSV.
  ///
  /// Takes precedence over paths configured with `--output-all` and `--output-selection`.
  #[cfg_attr(feature = "clap", clap(long, value_hint = ValueHint::FilePath, help_heading = "Output"))]
  pub output_confidence_tsv: Option<PathBuf>,

  /// Path to output iteration-statistics tracelog CSV (monitors convergence).
  ///
  /// Takes precedence over paths configured with `--output-all` and `--output-selection`.
  #[cfg_attr(feature = "clap", clap(long, visible_alias = "tracelog", value_hint = ValueHint::FilePath, help_heading = "Output"))]
  pub output_tracelog: Option<PathBuf>,

  #[cfg_attr(feature = "clap", clap(flatten))]
  pub output: OutputCoreArgs,

  /// Comma-separated list of outputs to produce with `--output-all`.
  ///
  /// Restricts which outputs `--output-all` writes. Special value `all` expands to every output
  /// available for this command. Requires `--output-all`. Per-file flags are always honored
  /// regardless of this selection.
  #[cfg_attr(
    feature = "clap",
    clap(long, value_delimiter = ',', requires = "output_all", help_heading = "Output")
  )]
  pub output_selection: Vec<TimetreeOutputSelection>,

  #[cfg_attr(feature = "clap", clap(flatten))]
  pub topology_order: TopologyOrderArgs,

  /// Random seed
  #[cfg_attr(feature = "clap", clap(long, visible_alias = "rng-seed"))]
  pub seed: Option<u64>,

  /// Use amino-acid alphabet (v0 compat, equivalent to `--alphabet=aa`)
  #[cfg_attr(feature = "clap", clap(long, hide = true))]
  pub aa: bool,

  /// Load a custom GTR model from file (not yet implemented)
  #[cfg_attr(feature = "clap", clap(long, hide = true))]
  #[cfg_attr(feature = "clap", clap(value_hint = ValueHint::FilePath))]
  pub custom_gtr: Option<PathBuf>,

  /// Method for clock filter outlier detection (not yet implemented)
  #[cfg_attr(feature = "clap", clap(long, hide = true))]
  pub clock_filter_method: Option<String>,

  /// Generations per year for coalescent model (not yet implemented)
  #[cfg_attr(feature = "clap", clap(long, hide = true))]
  pub gen_per_year: Option<f64>,

  /// Use greedy polytomy resolution (not yet implemented)
  #[cfg_attr(feature = "clap", clap(long, hide = true))]
  pub greedy_resolve: bool,

  /// Use stochastic polytomy resolution (not yet implemented)
  #[cfg_attr(feature = "clap", clap(long, hide = true))]
  pub stochastic_resolve: bool,
}
