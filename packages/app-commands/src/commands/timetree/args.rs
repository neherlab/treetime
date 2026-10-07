use crate::commands::shared::alignment::AlignmentArgs;
use crate::commands::shared::alphabet::AlphabetArgs;
use crate::commands::shared::branch_length_mode::BranchLengthModeCli;
use crate::commands::shared::config::ConfigArgs;
use crate::commands::shared::gap_fill::GapFillArgs;
use crate::commands::shared::metadata::{DateColumnArgs, MetadataIdArgs};
use crate::commands::shared::method_anc::MethodAncestralCli;
use crate::commands::shared::model::ModelArgs;
use crate::commands::shared::output_args::{DivergenceUnits, OutputCoreArgs, TimetreeOutputSelection};
use crate::commands::shared::required::missing_required_args;
use crate::commands::shared::reroot::RerootArgs;
use crate::commands::shared::seed::SeedArgs;
use crate::commands::shared::topology_order_args::TopologyOrderArgs;
use crate::commands::shared::tree_input::TreeDialectArgs;
#[cfg(feature = "clap")]
use clap::ValueHint;
use deser::{Deserialize, Serialize};
use eyre::Report;
use schemars::JsonSchema;
use smart_default::SmartDefault;
use std::path::PathBuf;
use treetime::ancestral::params::MethodAncestral;
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

#[derive(Debug, Clone)]
pub struct TreetimeTimetreeArgs {
  pub(crate) alignment: AlignmentArgs,
  pub(crate) tree: PathBuf,
  pub(crate) tree_dialect: TreeDialectArgs,
  #[expect(
    dead_code,
    reason = "VCF input is not implemented, see kb/issues/M-io-vcf-input-output-unimplemented.md"
  )]
  pub(crate) vcf_reference: Option<PathBuf>,
  pub(crate) metadata: Option<PathBuf>,
  pub(crate) metadata_id: MetadataIdArgs,
  pub(crate) date_column: DateColumnArgs,
  pub(crate) sequence_length: Option<usize>,
  pub(crate) clock_rate: Option<f64>,
  pub(crate) clock_std_dev: Option<f64>,
  pub(crate) branch_length_mode: BranchLengthMode,
  pub(crate) time_marginal: TimeMarginalMode,
  pub(crate) confidence: bool,
  pub(crate) max_grid_points: MaxGridPoints,
  #[expect(
    dead_code,
    reason = "parsed but not implemented, see kb/issues/M-cli-flags-parsed-but-ignored.md"
  )]
  pub(crate) keep_polytomies: bool,
  pub(crate) resolve_polytomies: bool,
  pub(crate) relax: Vec<f64>,
  pub(crate) max_iter: usize,
  pub(crate) coalescent: Option<f64>,
  pub(crate) coalescent_opt: bool,
  pub(crate) coalescent_skyline: bool,
  pub(crate) skyline_n_points: usize,
  pub(crate) skyline_stiffness: f64,
  pub(crate) coalescent_confidence: f64,
  pub(crate) n_branches_posterior: Option<usize>,
  #[expect(
    dead_code,
    reason = "parsed but not implemented, see kb/issues/M-cli-flags-parsed-but-ignored.md"
  )]
  pub(crate) tip_labels: bool,
  #[expect(
    dead_code,
    reason = "parsed but not implemented, see kb/issues/M-cli-flags-parsed-but-ignored.md"
  )]
  pub(crate) no_tip_labels: bool,
  pub(crate) clock_filter: f64,
  #[expect(
    dead_code,
    reason = "parsed but not implemented, see kb/issues/M-cli-flags-parsed-but-ignored.md"
  )]
  pub(crate) n_iqd: Option<f64>,
  pub(crate) reroot: RerootArgs,
  pub(crate) keep_root: bool,
  pub(crate) allow_negative_rate: bool,
  pub(crate) tip_slack: Option<f64>,
  pub(crate) covariation: bool,
  pub(crate) model_args: ModelArgs,
  #[expect(dead_code, reason = "see kb/issues/M-timetree-method-anc-ignored.md")]
  pub(crate) method_anc: MethodAncestral,
  pub(crate) alphabet_args: AlphabetArgs,
  pub(crate) dense: Option<bool>,
  pub(crate) gap_fill_args: GapFillArgs,
  #[expect(
    dead_code,
    reason = "parsed but not implemented, see kb/issues/M-cli-flags-parsed-but-ignored.md"
  )]
  pub(crate) zero_based: bool,
  pub(crate) include_leaves: bool,
  pub(crate) impute_missing_data: bool,
  pub(crate) report_ambiguous: bool,
  pub(crate) no_indels: bool,
  pub(crate) divergence_units: DivergenceUnits,
  pub(crate) output_augur_node_data: Option<PathBuf>,
  pub(crate) output_gtr: Option<PathBuf>,
  pub(crate) output_reconstructed_nuc_fasta: Option<PathBuf>,
  pub(crate) output_clock_model: Option<PathBuf>,
  pub(crate) output_clock_csv: Option<PathBuf>,
  pub(crate) output_confidence_tsv: Option<PathBuf>,
  pub(crate) output_tracelog: Option<PathBuf>,
  pub(crate) output_coalescent_tsv: Option<PathBuf>,
  pub(crate) output_coalescent_csv: Option<PathBuf>,
  pub(crate) output_coalescent_json: Option<PathBuf>,
  pub(crate) output: OutputCoreArgs,
  pub(crate) output_selection: Vec<TimetreeOutputSelection>,
  pub(crate) topology_order: TopologyOrderArgs,
  pub(crate) seed_args: SeedArgs,
  #[expect(
    dead_code,
    reason = "parsed but not implemented, see kb/issues/M-cli-flags-parsed-but-ignored.md"
  )]
  pub(crate) aa: bool,
  #[expect(
    dead_code,
    reason = "parsed but not implemented, see kb/issues/M-cli-flags-parsed-but-ignored.md"
  )]
  pub(crate) custom_gtr: Option<PathBuf>,
  #[expect(
    dead_code,
    reason = "parsed but not implemented, see kb/issues/M-cli-flags-parsed-but-ignored.md"
  )]
  pub(crate) clock_filter_method: Option<String>,
  pub(crate) gen_per_year: f64,
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
      .ok_or_else(|| missing_required_args::<TreetimeTimetreeArgsRaw>(&["tree"]))?;
    Ok(Self {
      alignment: raw.alignment,
      tree,
      tree_dialect: raw.tree_dialect,
      vcf_reference: raw.vcf_reference,
      metadata: raw.metadata,
      metadata_id: raw.metadata_id,
      date_column: raw.date_column,
      sequence_length: raw.sequence_length,
      clock_rate: raw.clock_rate,
      clock_std_dev: raw.clock_std_dev,
      branch_length_mode: raw.branch_length_mode.into(),
      time_marginal: raw.time_marginal.into(),
      confidence: raw.confidence,
      max_grid_points: raw.max_grid_points.unwrap_or_default(),
      keep_polytomies: raw.keep_polytomies,
      resolve_polytomies: raw.resolve_polytomies,
      relax: raw.relax,
      max_iter: raw.max_iter,
      coalescent: raw.coalescent,
      coalescent_opt: raw.coalescent_opt,
      coalescent_skyline: raw.coalescent_skyline,
      skyline_n_points: raw.skyline_n_points,
      skyline_stiffness: raw.skyline_stiffness,
      coalescent_confidence: raw.coalescent_confidence,
      n_branches_posterior: raw.n_branches_posterior,
      tip_labels: raw.tip_labels,
      no_tip_labels: raw.no_tip_labels,
      clock_filter: raw.clock_filter,
      n_iqd: raw.n_iqd,
      reroot: raw.reroot,
      keep_root: raw.keep_root,
      allow_negative_rate: raw.allow_negative_rate,
      tip_slack: raw.tip_slack,
      covariation: raw.covariation,
      model_args: raw.model_args,
      method_anc: raw.method_anc.into(),
      alphabet_args: raw.alphabet_args,
      dense: raw.dense,
      gap_fill_args: raw.gap_fill_args,
      zero_based: raw.zero_based,
      include_leaves: raw.include_leaves || raw.reconstruct_tip_states,
      impute_missing_data: raw.impute_missing_data || raw.reconstruct_tip_states,
      report_ambiguous: raw.report_ambiguous,
      no_indels: raw.no_indels,
      divergence_units: raw.divergence_units,
      output_augur_node_data: raw.output_augur_node_data,
      output_gtr: raw.output_gtr,
      output_reconstructed_nuc_fasta: raw.output_reconstructed_nuc_fasta,
      output_clock_model: raw.output_clock_model,
      output_clock_csv: raw.output_clock_csv,
      output_confidence_tsv: raw.output_confidence_tsv,
      output_tracelog: raw.output_tracelog,
      output_coalescent_tsv: raw.output_coalescent_tsv,
      output_coalescent_csv: raw.output_coalescent_csv,
      output_coalescent_json: raw.output_coalescent_json,
      output: raw.output,
      output_selection: raw.output_selection,
      topology_order: raw.topology_order,
      seed_args: raw.seed_args,
      aa: raw.aa,
      custom_gtr: raw.custom_gtr,
      clock_filter_method: raw.clock_filter_method,
      gen_per_year: raw.gen_per_year,
    })
  }
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

  /// Tree in Newick format.
  #[cfg_attr(feature = "clap", clap(long, short = 't', help_heading = "Input data"))]
  #[cfg_attr(feature = "clap", clap(value_hint = ValueHint::FilePath))]
  #[schemars(extend("x-path" = "input"))]
  pub tree: Option<PathBuf>,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub tree_dialect: TreeDialectArgs,

  /// Only for vcf input: fasta file of the sequence the VCF was mapped to.
  #[cfg_attr(feature = "clap", clap(long, short = 'r', help_heading = "Input data"))]
  #[cfg_attr(feature = "clap", clap(value_hint = ValueHint::FilePath))]
  #[schemars(extend("x-path" = "input"))]
  pub vcf_reference: Option<PathBuf>,

  /// CSV/TSV file with metadata including sampling dates
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
  /// via quadrature sum. Requires `--covariation` or `--clock-std-dev`.
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

  /// maximal number of iterations the inference cycle is run. For polytomy resolution and
  /// coalescence models max_iter should be at least 2
  #[default = 2]
  #[cfg_attr(feature = "clap", clap(long, default_value_t = TreetimeTimetreeArgsRaw::default().max_iter, help_heading = "Dating"))]
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

  /// ignore tips that don't follow a loose clock, 'clock-filter=number of inter-quartile ranges from
  /// regression'. Default=3.0, set to 0 to switch off.
  #[cfg_attr(
    feature = "clap",
    clap(long, default_value = "3.0", help_heading = "Molecular clock")
  )]
  #[default = 3.0]
  pub clock_filter: f64,

  /// Number of IQD (interquartile distance) for clock filter outlier detection
  #[cfg_attr(feature = "clap", clap(long, help_heading = "Molecular clock"))]
  pub n_iqd: Option<f64>,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub reroot: RerootArgs,

  /// don't reroot the tree. Otherwise, reroot to minimize the residual of the regression of
  /// root-to-tip distance and sampling time
  #[cfg_attr(feature = "clap", clap(long, conflicts_with_all = ["reroot", "reroot_tips"], help_heading = "Rooting"))]
  pub keep_root: bool,

  /// By default, rates are forced to be positive. For trees with little temporal signal it is advisable to remove this restriction to achieve essentially mid-point rooting.
  #[cfg_attr(feature = "clap", clap(long, help_heading = "Molecular clock"))]
  pub allow_negative_rate: bool,

  /// excess variance associated with terminal nodes accounting for overdispersion of the molecular
  /// clock
  #[cfg_attr(feature = "clap", clap(long, help_heading = "Molecular clock"))]
  pub tip_slack: Option<f64>,

  /// Account for covariation when estimating rates or rerooting using root-to-tip regression
  #[cfg_attr(feature = "clap", clap(long, help_heading = "Molecular clock"))]
  pub covariation: bool,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub model_args: ModelArgs,

  /// Method used for reconstructing ancestral sequences
  #[cfg_attr(feature = "clap", clap(long, value_enum, default_value_t = MethodAncestralCli::default(), help_heading = "Ancestral reconstruction"))]
  pub method_anc: MethodAncestralCli,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub alphabet_args: AlphabetArgs,

  /// Use dense representation for sequences (store full probability distributions)
  #[cfg_attr(feature = "clap", clap(long, help_heading = "Ancestral reconstruction"))]
  pub dense: Option<bool>,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub gap_fill_args: GapFillArgs,

  /// Zero-based mutation indexing
  #[cfg_attr(feature = "clap", clap(long, help_heading = "Input data"))]
  pub zero_based: bool,

  /// Emit reconstructed leaf (tip) sequences in addition to internal nodes.
  #[cfg_attr(feature = "clap", clap(long, help_heading = "Ancestral reconstruction"))]
  pub include_leaves: bool,

  /// Resolve ambiguous and unknown tip states (`N` and IUPAC codes such as `R`) to the most likely
  /// inferred state. Gaps are left as deletions.
  #[cfg_attr(feature = "clap", clap(long, help_heading = "Ancestral reconstruction"))]
  pub impute_missing_data: bool,

  /// v0-compatible alias for `--include-leaves --impute-missing-data`.
  #[cfg_attr(feature = "clap", clap(long, help_heading = "Ancestral reconstruction"))]
  pub reconstruct_tip_states: bool,

  /// Include branch mutations from or to the fully ambiguous nucleotide state `N`.
  ///
  /// By default these mutations are omitted from the branch mutations of the Newick and Nexus
  /// annotations and the Auspice JSON. Other ambiguity codes, such as `K` or `R`, are always
  /// reported. The MAT outputs store `N` as missing data and never contain these mutations.
  #[cfg_attr(feature = "clap", clap(long, help_heading = "Ancestral reconstruction"))]
  pub report_ambiguous: bool,

  /// Disable indel (insertion/deletion) contributions to branch-length
  /// optimization and branch-length distributions.
  ///
  /// When set, branch-length optimization uses substitution-only likelihood
  /// and timetree branch distributions exclude the Poisson indel term.
  /// Matches standard phylogenetic tools (RAxML, IQ-TREE, PhyML, BEAST)
  /// and enables v0 parity testing. Default: indels enabled.
  #[cfg_attr(feature = "clap", clap(long, help_heading = "Branch lengths"))]
  pub no_indels: bool,

  /// Units for divergence values in augur node data JSON and auspice output.
  ///
  /// `mutations-per-site` (default): branch divergence as substitutions per site.
  /// `mutations`: count of reconstructed substitutions per branch that change the state:
  /// gaps, `N`, and ambiguity codes compatible with the parent state do not count.
  /// Requires ancestral reconstruction (incompatible with `--branch-length-mode=input`).
  #[cfg_attr(feature = "clap", clap(long, value_enum, default_value_t = DivergenceUnits::default(), help_heading = "Output"))]
  pub divergence_units: DivergenceUnits,

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

  /// Path to output GTR model JSON.
  ///
  /// Takes precedence over paths configured with `--output-all` and `--output-selection`.
  #[cfg_attr(feature = "clap", clap(long, value_hint = ValueHint::FilePath, help_heading = "Output"))]
  #[schemars(extend("x-path" = "output"))]
  pub output_gtr: Option<PathBuf>,

  /// Path to output reconstructed ancestral-sequence nucleotide FASTA.
  ///
  /// The v1 equivalent of TreeTime v0's `ancestral_sequences.fasta`: internal-node sequences
  /// reconstructed by the marginal pass, plus reconstructed tip sequences when `--include-leaves`
  /// (or `--reconstruct-tip-states`) is set. `--impute-missing-data` resolves ambiguous tip states.
  ///
  /// Takes precedence over paths configured with `--output-all` and `--output-selection`.
  #[cfg_attr(feature = "clap", clap(long, value_hint = ValueHint::FilePath, help_heading = "Output"))]
  #[schemars(extend("x-path" = "output"))]
  pub output_reconstructed_nuc_fasta: Option<PathBuf>,

  /// Path to output clock model JSON.
  ///
  /// Takes precedence over paths configured with `--output-all` and `--output-selection`.
  #[cfg_attr(feature = "clap", clap(long, value_hint = ValueHint::FilePath, help_heading = "Output"))]
  #[schemars(extend("x-path" = "output"))]
  pub output_clock_model: Option<PathBuf>,

  /// Path to output clock regression CSV.
  ///
  /// One row per sample as the final clock model saw it: the date the regression used, marked
  /// `input` or `inferred`, the root-to-tip distance it regressed on, the date the model predicts
  /// from that distance, and whether the clock filter excluded the sample. Not written by default.
  /// Takes precedence over paths configured with `--output-all` and `--output-selection`.
  #[cfg_attr(feature = "clap", clap(long, value_hint = ValueHint::FilePath, help_heading = "Output"))]
  #[schemars(extend("x-path" = "output"))]
  pub output_clock_csv: Option<PathBuf>,

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

  /// Comma-separated list of outputs to produce with `--output-all`.
  ///
  /// Restricts which outputs `--output-all` writes. Special value `all` expands to every output
  /// available for this command. Requires `--output-all`. Per-file flags are always honored
  /// regardless of this selection.
  ///
  /// A selected output that the run has no data for is skipped without a message, for example
  /// the substitution model of a run that fits none. A per-file flag for such an output fails.
  #[cfg_attr(
    feature = "clap",
    clap(long, value_delimiter = ',', requires = "output_all", help_heading = "Output")
  )]
  pub output_selection: Vec<TimetreeOutputSelection>,

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

impl From<TimeMarginalModeCli> for TimeMarginalMode {
  fn from(mode: TimeMarginalModeCli) -> Self {
    match mode {
      TimeMarginalModeCli::Never => TimeMarginalMode::Never,
      TimeMarginalModeCli::Always => TimeMarginalMode::Always,
      TimeMarginalModeCli::OnlyFinal => TimeMarginalMode::OnlyFinal,
    }
  }
}
