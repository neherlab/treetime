use crate::commands::shared::alignment::AlignmentArgs;
use crate::commands::shared::alphabet::AlphabetArgs;
use crate::commands::shared::ambiguity::TimetreeAmbiguityArgs;
use crate::commands::shared::branch_length_mode::BranchLengthModeCli;
use crate::commands::shared::branch_lengths::BranchLengthArgs;
use crate::commands::shared::clock_filter::ClockFilterArgs;
use crate::commands::shared::config::ConfigArgs;
use crate::commands::shared::gap_fill::GapFillArgs;
use crate::commands::shared::gtr_output::GtrOutputArgs;
use crate::commands::shared::metadata::{DateColumnArgs, MetadataIdArgs};
use crate::commands::shared::method_anc::MethodAncArgs;
use crate::commands::shared::mirror_enum::impl_from_mirror;
use crate::commands::shared::model::{ModelArgs, PartitionRatesArgs};
use crate::commands::shared::output_args::{ClockOutputArgs, OutputCoreArgs, TimetreeOutputSelectionArgs};
use crate::commands::shared::reconstructed_fasta::ReconstructedFastaOutputArgs;
use crate::commands::shared::representation::RepresentationArgs;
use crate::commands::shared::required::{missing_required_args, missing_settings};
use crate::commands::shared::reroot::{RerootArgs, RootingArgs};
use crate::commands::shared::seed::SeedArgs;
use crate::commands::shared::tip_sequences::{ImputeMissingDataArgs, TipSequenceArgs};
use crate::commands::shared::topology_order_args::TopologyOrderArgs;
use crate::commands::shared::tree_input::{TreeArgs, TreeDialectArgs};
use crate::commands::shared::zero_based::ZeroBasedArgs;
use app_output::ambiguity_outputs::AmbiguityOutputs;
#[cfg(feature = "clap")]
use clap::ValueHint;
use derive_more::{Deref, DerefMut};
use deser::{Deserialize, Serialize};
use eyre::Report;
use schemars::JsonSchema;
use smart_default::SmartDefault;
use std::path::PathBuf;
use treetime::optimize::indel::IndelPolicy;
use treetime::optimize::params::BranchLengthMode;
use treetime::timetree::params::TimeMarginalMode;
use treetime_grid::MaxGridPoints;
use treetime_schema::{schema_defaults, skip_serializing_optionals};
use treetime_utils::make_error;

#[cfg(feature = "clap")]
fn parse_skyline_n_points(s: &str) -> Result<usize, String> {
  let n: usize = s.parse().map_err(|err| format!("'{s}' is not a valid number: {err}"))?;
  if n < 2 {
    return Err("skyline-n-points must be at least 2".to_owned());
  }
  Ok(n)
}

#[derive(Debug, Clone, Deref, DerefMut)]
pub struct TreetimeTimetreeArgs {
  #[deref]
  #[deref_mut]
  pub(crate) raw: TreetimeTimetreeArgsRaw,
  pub(crate) tree: PathBuf,
  pub(crate) branch_length_mode: BranchLengthMode,
  pub(crate) indel_policy: IndelPolicy,
  pub(crate) time_marginal: TimeMarginalMode,
  pub(crate) max_grid_points: MaxGridPoints,
  pub(crate) include_leaves: bool,
  pub(crate) impute_missing_data: bool,
  pub(crate) ambiguity_outputs: AmbiguityOutputs,
}

impl TryFrom<TreetimeTimetreeArgsRaw> for TreetimeTimetreeArgs {
  type Error = Report;

  fn try_from(raw: TreetimeTimetreeArgsRaw) -> Result<Self, Report> {
    if raw.plot_rtt.is_some() {
      return make_error!("--plot-rtt is not yet implemented");
    }
    if raw.plot_tree.is_some() {
      return make_error!("--plot-tree is not yet implemented");
    }
    let tree = raw
      .tree
      .tree
      .clone()
      .ok_or_else(|| missing_required_args::<TreetimeTimetreeArgsRaw>(&["tree"]))?;
    check_sequence_source(&raw)?;
    Ok(Self {
      tree,
      branch_length_mode: raw.branch_length_mode.into(),
      indel_policy: raw.branch_lengths.indel_policy(),
      time_marginal: raw.time_marginal.into(),
      max_grid_points: raw.max_grid_points.unwrap_or_default(),
      include_leaves: raw.tip_sequences.include_leaves || raw.tip_sequences.reconstruct_tip_states,
      impute_missing_data: raw.impute_missing_data.impute_missing_data || raw.tip_sequences.reconstruct_tip_states,
      ambiguity_outputs: raw.ambiguity.resolve()?,
      raw,
    })
  }
}

fn check_sequence_source(raw: &TreetimeTimetreeArgsRaw) -> Result<(), Report> {
  if !raw.alignment.alignment.is_empty() {
    return Ok(());
  }
  if raw.branch_length_mode == BranchLengthModeCli::Marginal {
    return Err(missing_settings(
      &["alignment", "branch_length_mode"],
      "--branch-length-mode=marginal requires an alignment: add --alignment or use --branch-length-mode=input",
    ));
  }
  if raw.sequence_length.is_none_or(|length| length == 0) {
    return Err(missing_settings(
      &["alignment", "sequence_length"],
      "--sequence-length is required without an alignment: the time inference needs the sequence length",
    ));
  }
  Ok(())
}

#[derive(Debug, Clone, SmartDefault, JsonSchema, Serialize, Deserialize)]
#[deser(skip_serializing_optionals)]
#[schemars(transform = skip_serializing_optionals)]
#[schemars(default, deny_unknown_fields)]
#[schemars(transform = schema_defaults::<Self>)]
#[deser(default, deny_unknown_fields)]
#[cfg_attr(feature = "clap", derive(clap::Parser))]
#[schemars(rename = "TimetreeConfig")]
pub struct TreetimeTimetreeArgsRaw {
  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(skip)]
  #[deser(skip)]
  pub config_args: ConfigArgs,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  #[schemars(extend("x-path" = "input"))]
  pub alignment: AlignmentArgs,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub tree: TreeArgs,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub tree_dialect: TreeDialectArgs,

  /// Only for vcf input: fasta file of the sequence the VCF was mapped to.
  #[cfg_attr(feature = "clap", clap(long, short = 'r', help_heading = "Input data"))]
  #[cfg_attr(feature = "clap", clap(value_hint = ValueHint::FilePath))]
  #[schemars(extend("x-path" = "input"))]
  pub vcf_reference: Option<PathBuf>,

  /// Metadata table (CSV or TSV) with one row per sample
  ///
  /// `--metadata-id-columns` names the column that links a row to a tree tip. The other columns hold the data of the
  /// sample, such as its sampling date or its traits.
  #[cfg_attr(
    feature = "clap",
    clap(long = "metadata", visible_alias = "dates", short = 'd', help_heading = "Input data")
  )]
  #[cfg_attr(feature = "clap", clap(value_hint = ValueHint::FilePath))]
  #[schemars(extend("x-path" = "input"))]
  pub metadata: Option<PathBuf>,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub metadata_id: MetadataIdArgs,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub date_column: DateColumnArgs,

  /// Length of the sequence, used to calculate expected variation in branch length. Not required if alignment is provided.
  #[cfg_attr(feature = "clap", clap(long, help_heading = "Input data"))]
  pub sequence_length: Option<usize>,

  /// If specified, the rate of the molecular clock won't be optimized.
  #[schemars(example = 0.001)]
  #[cfg_attr(feature = "clap", clap(long, help_heading = "Molecular clock"))]
  pub clock_rate: Option<f64>,

  /// Standard deviation of the provided clock rate estimate
  #[cfg_attr(feature = "clap", clap(long, help_heading = "Molecular clock"))]
  pub clock_std_dev: Option<f64>,

  /// If set to 'input', the provided branch length will be used without modification. Branch lengths optimized by treetime are only accurate at short evolutionary distances.
  #[cfg_attr(feature = "clap", clap(long, value_enum, default_value_t = BranchLengthModeCli::default(), help_heading = "Branch lengths"))]
  pub branch_length_mode: BranchLengthModeCli,

  /// Control when marginal time distributions are used for output.
  ///
  /// All modes use marginal inference during optimization. The mode controls whether
  /// confidence intervals are extracted from the resulting distributions:
  ///
  /// - `never`: no confidence interval output (default)
  /// - `always`: write confidence intervals from distributions computed during optimization
  /// - `only-final`: run one extra inference pass after optimization, then write confidence intervals
  #[cfg_attr(feature = "clap", clap(long, value_enum, default_value_t = TimeMarginalModeCli::default(), help_heading = "Dating"))]
  pub time_marginal: TimeMarginalModeCli,

  /// Add rate-uncertainty to confidence intervals.
  ///
  /// `--time-marginal=always` and `only-final` already write mutation-stochasticity CIs.
  /// This flag adds rate-uncertainty CIs (re-runs inference at rate +/- sigma), combined
  /// via quadrature sum. The rate std. dev. is `--clock-std-dev` if given, else estimated by a Gaussian
  /// approximation of the timetree, which needs the alignment or `--sequence-length`.
  /// When set with `--time-marginal=never` (default), automatically promotes to `only-final`.
  #[cfg_attr(feature = "clap", clap(long, help_heading = "Dating"))]
  pub confidence: bool,

  /// Largest number of points of one probability grid during time inference.
  ///
  /// A run stops with an error when a grid would need more points. The value bounds single grids, not the total
  /// memory of a run. When unset, the run takes the limit of the app or server that runs it, otherwise 1000000.
  #[schemars(example = 1000000)]
  #[cfg_attr(feature = "clap", clap(long, help_heading = "Dating"))]
  pub max_grid_points: Option<MaxGridPoints>,

  /// Don't resolve polytomies using temporal information.
  #[cfg_attr(feature = "clap", clap(long, help_heading = "Polytomies"))]
  pub keep_polytomies: bool,

  /// Resolve polytomies using temporal information
  #[cfg_attr(feature = "clap", clap(long, help_heading = "Polytomies"))]
  pub resolve_polytomies: bool,

  /// use an autocorrelated molecular clock. Strength of the gaussian priors on branch specific rate
  /// deviation and the coupling of parent and offspring rates can be specified e.g. as --relax 1.0
  /// 0.5. Values around 1.0 correspond to weak priors, larger values constrain rate deviations more
  /// strongly. Coupling 0 (--relax 1.0 0) corresponds to an un-correlated clock.
  #[schemars(example = [1.0, 0.0])]
  #[cfg_attr(feature = "clap", clap(long, num_args = 2, value_names = ["SLACK", "COUPLING"], help_heading = "Molecular clock"))]
  pub relax: Vec<f64>,

  /// Maximum number of iterations of the inference loop
  ///
  /// Polytomy resolution and coalescent models of the time inference need at least 2 iterations.
  #[default = 2]
  #[cfg_attr(feature = "clap", clap(long, default_value_t = TreetimeTimetreeArgsRaw::default().max_iter, help_heading = "Branch lengths"))]
  pub max_iter: usize,

  /// Coalescent time scale in years.
  ///
  /// Sensible values are on the order of the time from the root to the tips and are given in units of time.
  #[schemars(example = 1.0)]
  #[cfg_attr(feature = "clap", clap(long, help_heading = "Coalescent prior"))]
  pub coalescent: Option<f64>,

  /// Optimize coalescent time scale Tc to maximize coalescent likelihood.
  ///
  /// When set, TreeTime finds the optimal constant Tc analytically (closed-form maximum
  /// of the coalescent likelihood). This is similar to Python v0's `--coalescent=opt`,
  /// which used a numerical search.
  #[cfg_attr(
    feature = "clap",
    clap(
      long,
      conflicts_with = "coalescent",
      conflicts_with = "coalescent_skyline",
      help_heading = "Coalescent prior"
    )
  )]
  pub coalescent_opt: bool,

  /// Use skyline coalescent model instead of constant Tc.
  ///
  /// Estimates a piecewise linear coalescent rate history. Requires --skyline-n-points to specify
  /// the number of grid points.
  #[cfg_attr(feature = "clap", clap(long, help_heading = "Coalescent prior"))]
  #[cfg_attr(
    feature = "clap",
    clap(conflicts_with = "coalescent", conflicts_with = "coalescent_opt")
  )]
  pub coalescent_skyline: bool,

  /// Number of grid points in skyline coalescent model.
  ///
  /// Only used when --coalescent-skyline is set. Defines how many piecewise linear segments
  /// are used to model Tc(t) over time. Must be at least 2. Matches Python v0's default.
  #[default = 20]
  #[cfg_attr(feature = "clap", clap(long, default_value_t = TreetimeTimetreeArgsRaw::default().skyline_n_points, help_heading = "Coalescent prior"))]
  #[cfg_attr(feature = "clap", clap(value_parser = parse_skyline_n_points))]
  pub skyline_n_points: usize,

  /// Smoothing stiffness for the skyline coalescent.
  ///
  /// Penalizes log-fold-changes of the coalescent time scale Tc between adjacent
  /// skyline segments: with z = ln Tc, the objective adds
  /// `(stiffness/2) * Σ (ln(Tc_{i+1}/Tc_i))^2`. Because it acts on log Tc, the
  /// stiffness is dimensionless and scale-independent. Larger values enforce a
  /// smoother Tc(t). Only used when --coalescent-skyline is set.
  #[default = 2.0]
  #[cfg_attr(feature = "clap", clap(long, default_value_t = TreetimeTimetreeArgsRaw::default().skyline_stiffness, help_heading = "Coalescent prior"))]
  pub skyline_stiffness: f64,

  /// Confidence level for coalescent time scale (Tc) bands, in standard deviations.
  ///
  /// Applies to every inferred coalescent mode (constant, --coalescent-opt, and
  /// --coalescent-skyline). The band spans `Tc * exp(±confidence * σ)`, where `σ` is
  /// the standard deviation of `ln Tc` from the coalescent likelihood curvature. A fixed
  /// --coalescent value is not inferred and therefore has no band.
  #[default = 2.0]
  #[cfg_attr(feature = "clap", clap(long, default_value_t = TreetimeTimetreeArgsRaw::default().coalescent_confidence, help_heading = "Coalescent prior"))]
  pub coalescent_confidence: f64,

  /// add posterior LH to coalescent model: use the posterior probability distributions of
  /// divergence times for estimating the number of branches when calculating the coalescent
  /// mergerrate or use inferred time before present (default).
  #[cfg_attr(feature = "clap", clap(long, hide = true, help_heading = "Dating"))]
  pub n_branches_posterior: Option<usize>,

  /// filename to save the plot to. Suffix will determine format (choices pdf, png, svg,
  /// default=pdf)
  #[cfg_attr(feature = "clap", clap(long, hide = true, help_heading = "Plots"))]
  #[cfg_attr(feature = "clap", clap(value_hint = ValueHint::FilePath))]
  #[schemars(extend("x-path" = "output"))]
  pub plot_tree: Option<PathBuf>,

  /// filename to save the plot to. Suffix will determine format (choices pdf, png, svg,
  /// default=pdf)
  #[cfg_attr(feature = "clap", clap(long, hide = true, help_heading = "Plots"))]
  #[cfg_attr(feature = "clap", clap(value_hint = ValueHint::FilePath))]
  #[schemars(extend("x-path" = "output"))]
  pub plot_rtt: Option<PathBuf>,

  /// add tip labels (default for small trees with <30 leaves)
  #[cfg_attr(feature = "clap", clap(long, help_heading = "Plots"))]
  pub tip_labels: bool,

  /// don't show tip labels (default for trees with >=30 leaves)
  #[cfg_attr(feature = "clap", clap(long, help_heading = "Plots"))]
  pub no_tip_labels: bool,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub clock_filter: ClockFilterArgs,

  /// Number of IQD (interquartile distance) for clock filter outlier detection
  #[cfg_attr(feature = "clap", clap(long, help_heading = "Molecular clock"))]
  pub n_iqd: Option<f64>,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub reroot: RerootArgs,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub rooting: RootingArgs,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub model_args: ModelArgs,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub partition_rates_args: PartitionRatesArgs,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub method_anc: MethodAncArgs,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub alphabet_args: AlphabetArgs,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub representation: RepresentationArgs,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub gap_fill_args: GapFillArgs,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub zero_based: ZeroBasedArgs,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub tip_sequences: TipSequenceArgs,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub impute_missing_data: ImputeMissingDataArgs,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub branch_lengths: BranchLengthArgs,

  /// Path to output augur-compatible node data JSON.
  ///
  /// Contains per-node dates, branch lengths, clock model parameters, confidence
  /// intervals, and divergence metrics. The output is compatible with augur
  /// export v2 --node-data for Nextstrain pipeline integration.
  ///
  /// Takes precedence over paths configured with `--output-all` and `--output-selection`.
  #[cfg_attr(feature = "clap", clap(long, value_hint = ValueHint::FilePath, help_heading = "Output"))]
  #[schemars(extend("x-path" = "output"))]
  pub output_augur_node_data: Option<PathBuf>,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub gtr_output: GtrOutputArgs,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub reconstructed_fasta_output: ReconstructedFastaOutputArgs,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub clock_output: ClockOutputArgs,

  /// Path to output date-confidence-interval TSV.
  ///
  /// Takes precedence over paths configured with `--output-all` and `--output-selection`.
  #[cfg_attr(feature = "clap", clap(long, value_hint = ValueHint::FilePath, help_heading = "Output"))]
  #[schemars(extend("x-path" = "output"))]
  pub output_confidence_tsv: Option<PathBuf>,

  /// Path to output iteration-statistics tracelog CSV (monitors convergence).
  ///
  /// Takes precedence over paths configured with `--output-all` and `--output-selection`.
  #[cfg_attr(feature = "clap", clap(long, visible_alias = "tracelog", value_hint = ValueHint::FilePath, help_heading = "Output"))]
  #[schemars(extend("x-path" = "output"))]
  pub output_tracelog: Option<PathBuf>,

  /// Path to output the coalescent time scale as a flat TSV (one row per skyline segment).
  ///
  /// Written when a coalescent model is set (`--coalescent`, `--coalescent-opt`, or
  /// `--coalescent-skyline`). A fixed `--coalescent` writes one band-less segment over the tree
  /// span. Takes precedence over paths configured with `--output-all` and `--output-selection`.
  #[cfg_attr(feature = "clap", clap(long, value_hint = ValueHint::FilePath, help_heading = "Output"))]
  #[schemars(extend("x-path" = "output"))]
  pub output_coalescent_tsv: Option<PathBuf>,

  /// Path to output the coalescent time scale as a flat CSV (one row per skyline segment).
  ///
  /// Written when a coalescent model is set. Takes precedence over paths configured with
  /// `--output-all` and `--output-selection`.
  #[cfg_attr(feature = "clap", clap(long, value_hint = ValueHint::FilePath, help_heading = "Output"))]
  #[schemars(extend("x-path" = "output"))]
  pub output_coalescent_csv: Option<PathBuf>,

  /// Path to output the coalescent time scale as a rich JSON document (inputs + segments).
  ///
  /// Written when a coalescent model is set. Takes precedence over paths configured with
  /// `--output-all` and `--output-selection`.
  #[cfg_attr(feature = "clap", clap(long, value_hint = ValueHint::FilePath, help_heading = "Output"))]
  #[schemars(extend("x-path" = "output"))]
  pub output_coalescent_json: Option<PathBuf>,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub output: OutputCoreArgs,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub output_selection: TimetreeOutputSelectionArgs,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub ambiguity: TimetreeAmbiguityArgs,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub topology_order: TopologyOrderArgs,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub seed_args: SeedArgs,

  /// Use amino-acid alphabet (v0 compat, equivalent to `--alphabet=aa`)
  #[cfg_attr(feature = "clap", clap(long, hide = true, help_heading = "Input data"))]
  pub aa: bool,

  /// Load a custom GTR model from file (not yet implemented)
  #[cfg_attr(feature = "clap", clap(long, hide = true, help_heading = "Substitution model"))]
  #[cfg_attr(feature = "clap", clap(value_hint = ValueHint::FilePath))]
  #[schemars(extend("x-path" = "input"))]
  pub custom_gtr: Option<PathBuf>,

  /// Method for clock filter outlier detection (not yet implemented)
  #[cfg_attr(feature = "clap", clap(long, hide = true, help_heading = "Molecular clock"))]
  pub clock_filter_method: Option<String>,

  /// Generations per year for converting the coalescent time scale Tc into an effective
  /// population size.
  ///
  /// The coalescent output reports an effective population size `N_e = Tc * gen_per_year`. Tc is
  /// already in calendar years, so this factor rescales it into generation units, the standard
  /// axis of a skyline plot. Only affects the reported `N_e`; it does not enter the inference.
  #[default = 50.0]
  #[cfg_attr(
    feature = "clap",
    clap(long, default_value_t = TreetimeTimetreeArgsRaw::default().gen_per_year, help_heading = "Coalescent prior")
  )]
  pub gen_per_year: f64,
}

#[derive(Copy, Debug, Clone, PartialEq, Eq, PartialOrd, Ord, SmartDefault, JsonSchema, Serialize, Deserialize)]
#[cfg_attr(feature = "clap", derive(clap::ValueEnum))]
#[schemars(rename_all = "kebab-case")]
#[deser(rename_all = "kebab-case")]
#[schemars(rename = "TimeMarginalMode")]
pub enum TimeMarginalModeCli {
  #[default]
  Never,
  Always,
  OnlyFinal,
}

impl_from_mirror!(TimeMarginalModeCli => TimeMarginalMode { Never, Always, OnlyFinal });
