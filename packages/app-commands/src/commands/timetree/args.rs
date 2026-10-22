use crate::commands::shared::alignment::AlignmentArgs;
use crate::commands::shared::alphabet::AlphabetArgs;
use crate::commands::shared::ambiguity::TimetreeAmbiguityArgs;
use crate::commands::shared::argument_error::argument_error;
use crate::commands::shared::branch_lengths::BranchLengthArgs;
use crate::commands::shared::clock_filter::ClockFilterArgs;
use crate::commands::shared::config::ConfigArgs;
use crate::commands::shared::dependent_settings::{check_applies, setting_flag};
use crate::commands::shared::documented_default::documented_default;
use crate::commands::shared::gap_fill::GapFillArgs;
use crate::commands::shared::gtr_output::GtrOutputArgs;
use crate::commands::shared::jobs::JobsArgs;
use crate::commands::shared::metadata::{DateColumnArgs, MetadataIdArgs};
use crate::commands::shared::model::{ModelArgs, PartitionRatesArgs};
use crate::commands::shared::output_args::{ClockOutputArgs, OutputCoreArgs, TimetreeOutputSelectionArgs};
use crate::commands::shared::reconstructed_fasta::ReconstructedFastaOutputArgs;
use crate::commands::shared::representation::RepresentationArgs;
use crate::commands::shared::required::missing_required_args;
use crate::commands::shared::reroot::{RerootArgs, RootingArgs};
use crate::commands::shared::seed::SeedArgs;
use crate::commands::shared::sequence_inputs::MissingSequenceArgs;
use crate::commands::shared::tip_sequences::{ImputeMissingDataArgs, TipSequenceArgs};
use crate::commands::shared::topology_order_args::TopologyOrderArgs;
use crate::commands::shared::tree_input::{TreeArgs, TreeDialectArgs};
use crate::commands::shared::zero_based::ZeroBasedArgs;
use crate::commands::timetree::branch_length_mode::{BranchLengthModeCli, ResolvedBranchLengthMode};
use crate::commands::timetree::coalescent::CoalescentArg;
#[cfg(feature = "clap")]
use crate::commands::timetree::coalescent::parse_coalescent;
use app_output::ambiguity_outputs::AmbiguityOutputs;
#[cfg(feature = "clap")]
use clap::ValueHint;
use derive_more::{Deref, DerefMut};
use deser::{Deserialize, Serialize};
use eyre::Report;
use schemars::JsonSchema;
use smart_default::SmartDefault;
use std::path::PathBuf;
use treetime::gtr::get_gtr::GtrModel;
use treetime::optimize::indel::IndelPolicy;
use treetime_grid::MaxGridPoints;
use treetime_schema::{schema_defaults, skip_serializing_optionals};
use treetime_utils::o;

pub const N_SKYLINE_DEFAULT: usize = 20;

const N_SKYLINE_MIN: usize = 2;

pub const SKYLINE_STIFFNESS_DEFAULT: f64 = 2.0;

pub const COALESCENT_CONFIDENCE_DEFAULT: f64 = 2.0;

pub const GEN_PER_YEAR_DEFAULT: f64 = 50.0;

#[derive(Debug, Clone, Deref, DerefMut)]
pub struct TreetimeTimetreeArgs {
  #[deref]
  #[deref_mut]
  pub(crate) raw: TreetimeTimetreeArgsRaw,
  pub(crate) tree: PathBuf,
  pub(crate) indel_policy: IndelPolicy,
  pub(crate) max_grid_points: MaxGridPoints,
  pub(crate) include_leaves: bool,
  pub(crate) impute_missing_data: bool,
  pub(crate) ambiguity_outputs: AmbiguityOutputs,
  pub(crate) model: GtrModel,
}

impl TryFrom<TreetimeTimetreeArgsRaw> for TreetimeTimetreeArgs {
  type Error = Report;

  fn try_from(raw: TreetimeTimetreeArgsRaw) -> Result<Self, Report> {
    let tree = raw
      .tree
      .tree
      .clone()
      .ok_or_else(|| missing_required_args::<TreetimeTimetreeArgsRaw>(&["tree"]))?;
    check_sequence_source(&raw)?;
    check_dependent_settings(&raw)?;
    Ok(Self {
      tree,
      indel_policy: raw.branch_lengths.indel_policy(),
      max_grid_points: raw.max_grid_points.unwrap_or_default(),
      include_leaves: raw.tip_sequences.include_leaves || raw.tip_sequences.reconstruct_tip_states,
      impute_missing_data: raw.impute_missing_data.impute_missing_data || raw.tip_sequences.reconstruct_tip_states,
      ambiguity_outputs: raw.ambiguity.resolve()?,
      model: raw.model_args.model()?,
      raw,
    })
  }
}

pub(crate) fn check_branch_length_settings(
  raw: &TreetimeTimetreeArgsRaw,
  resolved: &ResolvedBranchLengthMode,
) -> Result<(), Report> {
  for (setting, given) in [
    ("representation", raw.representation.is_given()),
    ("partition_rates", raw.partition_rates_args.is_given()),
    ("ambiguity_threshold", raw.ambiguity.ambiguity_threshold.is_some()),
  ] {
    if given {
      resolved.require_marginal(&[setting], &setting_flag(setting))?;
    }
  }
  Ok(())
}

fn check_sequence_source(raw: &TreetimeTimetreeArgsRaw) -> Result<(), Report> {
  if !raw.alignment.alignment.is_empty() {
    return check_applies(
      "sequence_length",
      raw.sequence_length.is_some(),
      "alignment",
      false,
      || {
        o!(
          "applies only without --alignment: the run takes the sequence length from the alignment, and an \
            alignment of variable sites with a separate genome length is not supported"
        )
      },
    );
  }
  if raw.branch_length_mode == BranchLengthModeCli::Marginal {
    return Err(argument_error(
      &["alignment", "branch_length_mode"],
      "--branch-length-mode=marginal requires an alignment: add --alignment or use --branch-length-mode=input",
    ));
  }
  if raw.sequence_length.is_none_or(|length| length == 0) {
    return Err(argument_error(
      &["alignment", "sequence_length"],
      "--sequence-length is required without an alignment: the time inference needs the sequence length",
    ));
  }
  Ok(())
}

fn check_dependent_settings(raw: &TreetimeTimetreeArgsRaw) -> Result<(), Report> {
  check_coalescent_settings(raw)?;
  check_applies(
    "clock_std_dev",
    raw.clock_std_dev.is_some(),
    "confidence",
    raw.confidence,
    || o!("applies only with --confidence: only the date intervals use the clock rate uncertainty"),
  )?;
  raw
    .partition_rates_args
    .check_partitions(raw.alignment.alignment.len())?;
  if raw.branch_length_mode != BranchLengthModeCli::Auto {
    let has_alignment = !raw.alignment.alignment.is_empty();
    check_branch_length_settings(
      raw,
      &ResolvedBranchLengthMode::new(raw.branch_length_mode, has_alignment, None),
    )?;
  }
  Ok(())
}

fn check_coalescent_settings(raw: &TreetimeTimetreeArgsRaw) -> Result<(), Report> {
  let coalescent = raw.coalescent;
  let actual = || coalescent.map_or_else(|| o!("not set"), |coalescent| coalescent.to_string());
  let skyline = matches!(coalescent, Some(CoalescentArg::Skyline));
  for (setting, given) in [
    ("n_skyline", raw.n_skyline.is_some()),
    ("skyline_stiffness", raw.skyline_stiffness.is_some()),
  ] {
    check_applies(setting, given, "coalescent", skyline, || {
      format!("applies only to --coalescent skyline, but --coalescent is {}", actual())
    })?;
  }
  if let Some(n_skyline) = raw.n_skyline.filter(|&n_skyline| n_skyline < N_SKYLINE_MIN) {
    return Err(argument_error(
      &["n_skyline"],
      format!("--n-skyline is {n_skyline}, but the skyline needs at least {N_SKYLINE_MIN} grid points"),
    ));
  }
  let fitted = matches!(coalescent, Some(CoalescentArg::Opt | CoalescentArg::Skyline));
  check_applies(
    "coalescent_confidence",
    raw.coalescent_confidence.is_some(),
    "coalescent",
    fitted,
    || {
      format!(
        "applies only to --coalescent opt and --coalescent skyline, but --coalescent is {}",
        actual()
      )
    },
  )?;
  check_applies(
    "gen_per_year",
    raw.gen_per_year.is_some(),
    "coalescent",
    coalescent.is_some(),
    || o!("applies only with --coalescent, which is not set"),
  )
}

impl TreetimeTimetreeArgsRaw {
  pub fn n_skyline(&self) -> usize {
    self.n_skyline.unwrap_or(N_SKYLINE_DEFAULT)
  }

  pub fn skyline_stiffness(&self) -> f64 {
    self.skyline_stiffness.unwrap_or(SKYLINE_STIFFNESS_DEFAULT)
  }

  pub fn coalescent_confidence(&self) -> f64 {
    self.coalescent_confidence.unwrap_or(COALESCENT_CONFIDENCE_DEFAULT)
  }

  pub fn gen_per_year(&self) -> f64 {
    self.gen_per_year.unwrap_or(GEN_PER_YEAR_DEFAULT)
  }
}

#[derive(Debug, Clone, SmartDefault, JsonSchema, Serialize, Deserialize)]
#[deser(skip_serializing_optionals)]
#[schemars(transform = skip_serializing_optionals)]
#[schemars(default, deny_unknown_fields)]
#[schemars(transform = schema_defaults::<Self>)]
#[deser(default, deny_unknown_fields)]
#[cfg_attr(feature = "clap", derive(clap::Parser))]
#[cfg_attr(
  feature = "clap",
  command(
    mut_arg("n_skyline", |arg| documented_default(&N_SKYLINE_DEFAULT).help(arg)),
    mut_arg("skyline_stiffness", |arg| documented_default(&SKYLINE_STIFFNESS_DEFAULT).help(arg)),
    mut_arg("coalescent_confidence", |arg| documented_default(&COALESCENT_CONFIDENCE_DEFAULT).help(arg)),
    mut_arg("gen_per_year", |arg| documented_default(&GEN_PER_YEAR_DEFAULT).help(arg)),
  )
)]
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
  pub newick_dialect: TreeDialectArgs,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub missing_sequences: MissingSequenceArgs,

  /// Metadata table (CSV or TSV) with one row per sample
  ///
  /// `--metadata-id-columns` names the column that links a row to a tree tip. The other columns hold the data of the
  /// sample, such as its sampling date or its traits.
  #[cfg_attr(feature = "clap", clap(long, short = 'm', help_heading = "Input data"))]
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

  /// Length of the sequences, used to calculate the expected variation of the branch lengths
  ///
  /// Required without `--alignment`, and applies only then: with an alignment, the run takes the length from the
  /// alignment.
  #[cfg_attr(feature = "clap", clap(long, help_heading = "Input data"))]
  pub sequence_length: Option<usize>,

  /// If specified, the rate of the molecular clock won't be optimized.
  #[schemars(example = 0.001)]
  #[cfg_attr(feature = "clap", clap(long, help_heading = "Molecular clock"))]
  pub clock_rate: Option<f64>,

  /// Standard deviation of the provided clock rate estimate
  ///
  /// Applies only with `--confidence`, whose date intervals include the clock rate uncertainty.
  #[cfg_attr(feature = "clap", clap(long, help_heading = "Molecular clock"))]
  pub clock_std_dev: Option<f64>,

  /// How the branch lengths of the time inference are obtained
  ///
  /// `auto` keeps the input branch lengths without an alignment, and also when the longest input branch is longer
  /// than 0.1 substitutions per site, because the re-estimated branch lengths are accurate only at short evolutionary
  /// distances; otherwise it re-estimates them as `marginal` does. When no input branch has a length, `auto` chooses
  /// `marginal`. The run logs the choice.
  #[cfg_attr(feature = "clap", clap(long, value_enum, default_value_t = BranchLengthModeCli::default(), help_heading = "Branch lengths"))]
  pub branch_length_mode: BranchLengthModeCli,

  /// Estimate the 90% intervals of the node dates
  ///
  /// Runs a final marginal round of the time inference after the iterations and writes the 90% interval of every
  /// dated node. The intervals include the clock rate uncertainty when one is known: `--clock-std-dev` for a given
  /// rate, or the Gaussian estimate of an estimated rate. The rate uncertainty is added by inferring the dates again
  /// at the rate plus and minus one standard deviation and combining both parts in quadrature.
  #[cfg_attr(feature = "clap", clap(long, help_heading = "Dating"))]
  pub confidence: bool,

  /// Largest number of points of one probability grid during time inference.
  ///
  /// A run stops with an error when a grid would need more points. The value bounds single grids, not the total
  /// memory of a run. When unset, the run takes the limit of the app or server that runs it, otherwise 1000000.
  #[schemars(example = 1000000)]
  #[cfg_attr(feature = "clap", clap(long, help_heading = "Dating"))]
  pub max_grid_points: Option<MaxGridPoints>,

  /// Keep the polytomies of the input tree
  ///
  /// Without this flag, each polytomy is resolved into binary mergers by the node dates. The resolution is
  /// stochastic: the run logs the seed it draws, and `--seed` reproduces the run.
  #[cfg_attr(feature = "clap", clap(long, help_heading = "Polytomies"))]
  pub keep_polytomies: bool,

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

  /// Coalescent prior: a fixed time scale in substitutions per site (as in TreeTime 0.x and augur, for example 0.01)
  /// or as a time (10y, 26w, 90d); opt fits a constant time scale; skyline fits a piecewise-constant one
  ///
  /// A time is a positive number directly followed by the unit y (years), w (weeks) or d (days). A time scale in
  /// substitutions per site is converted to years with the clock rate of each iteration. Without this setting the node
  /// times have no coalescent prior.
  #[schemars(example = 0.01)]
  #[cfg_attr(
    feature = "clap",
    clap(long, value_parser = parse_coalescent, value_name = "TC|opt|skyline", help_heading = "Coalescent prior")
  )]
  pub coalescent: Option<CoalescentArg>,

  /// Number of grid points in skyline coalescent model.
  ///
  /// Applies only to --coalescent skyline. Defines how many piecewise linear segments
  /// are used to model Tc(t) over time. Must be at least 2.
  #[cfg_attr(feature = "clap", clap(long, help_heading = "Coalescent prior"))]
  #[schemars(range(min = N_SKYLINE_MIN), transform = documented_default(&N_SKYLINE_DEFAULT))]
  pub n_skyline: Option<usize>,

  /// Smoothing stiffness for the skyline coalescent.
  ///
  /// Penalizes log-fold-changes of the coalescent time scale Tc between adjacent
  /// skyline segments: with z = ln Tc, the objective adds
  /// `(stiffness/2) * Σ (ln(Tc_{i+1}/Tc_i))^2`. Because it acts on log Tc, the
  /// stiffness is dimensionless and scale-independent. Larger values enforce a
  /// smoother Tc(t). Applies only to --coalescent skyline.
  #[cfg_attr(feature = "clap", clap(long, help_heading = "Coalescent prior"))]
  #[schemars(transform = documented_default(&SKYLINE_STIFFNESS_DEFAULT))]
  pub skyline_stiffness: Option<f64>,

  /// Confidence level for coalescent time scale (Tc) bands, in standard deviations.
  ///
  /// Applies only to the inferred coalescent priors (--coalescent opt and --coalescent skyline).
  /// The band spans `Tc * exp(±confidence * σ)`, where `σ` is
  /// the standard deviation of `ln Tc` from the coalescent likelihood curvature. A fixed
  /// time scale is not inferred and therefore has no band.
  #[cfg_attr(feature = "clap", clap(long, help_heading = "Coalescent prior"))]
  #[schemars(transform = documented_default(&COALESCENT_CONFIDENCE_DEFAULT))]
  pub coalescent_confidence: Option<f64>,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub clock_filter: ClockFilterArgs,

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
  #[cfg_attr(feature = "clap", clap(long, value_hint = ValueHint::FilePath, help_heading = "Output"))]
  #[schemars(extend("x-path" = "output"))]
  pub output_tracelog: Option<PathBuf>,

  /// Path to output the coalescent time scale as a flat TSV (one row per skyline segment).
  ///
  /// Written when a coalescent prior is set with `--coalescent`. A fixed time scale writes one
  /// band-less segment over the tree span. Takes precedence over paths configured with `--output-all` and `--output-selection`.
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

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(skip)]
  #[deser(skip)]
  pub jobs_args: JobsArgs,

  /// Generations per year for converting the coalescent time scale Tc into an effective
  /// population size.
  ///
  /// The coalescent output reports an effective population size `N_e = Tc * gen_per_year`. Tc is
  /// already in calendar years, so this factor rescales it into generation units, the standard
  /// axis of a skyline plot. Only affects the reported `N_e`; it does not enter the inference. Applies only with
  /// `--coalescent`.
  #[cfg_attr(feature = "clap", clap(long, help_heading = "Coalescent prior"))]
  #[schemars(transform = documented_default(&GEN_PER_YEAR_DEFAULT))]
  pub gen_per_year: Option<f64>,
}
