use crate::commands::ancestral::args::SampleModeCli;
use crate::commands::shared::alignment::AlignmentArgs;
use crate::commands::shared::alphabet::AlphabetArgs;
use crate::commands::shared::config::ConfigArgs;
use crate::commands::shared::gap_fill::GapFillArgs;
use crate::commands::shared::method_anc::MethodAncestralCli;
use crate::commands::shared::model::ModelArgs;
use crate::commands::shared::output_args::{HomoplasyOutputSelection, OutputCoreArgs};
use crate::commands::shared::required::missing_required_args;
use crate::commands::shared::seed::SeedArgs;
use crate::commands::shared::topology_order_args::TopologyOrderArgs;
#[cfg(feature = "clap")]
use clap::ValueHint;
use eyre::Report;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use serde_with::skip_serializing_none;
use smart_default::SmartDefault;
use std::path::PathBuf;
use treetime::ancestral::params::{AncestralParams, MethodAncestral};
use treetime::make_error;
use treetime::partition::marginal::sample::SampleMode;

pub fn homoplasy_ancestral_params(args: &TreetimeHomoplasyArgs, seed: u64) -> AncestralParams {
  AncestralParams {
    method: args.method_anc,
    model: args.model_args.model_name(),
    dense: args.dense,
    include_leaves: false,
    report_ambiguous: true,
    impute_missing_data: args.impute_missing_data,
    gtr_iterations: args.gtr_iterations,
    site_specific_gtr: args.site_specific_gtr,
    seed,
    sample_from_profile: args.sample_from_profile,
  }
}

#[derive(Debug, Clone)]
pub struct TreetimeHomoplasyArgs {
  pub(crate) alignment: AlignmentArgs,
  pub(crate) tree: PathBuf,
  pub(crate) alphabet_args: AlphabetArgs,
  pub(crate) model_args: ModelArgs,
  pub(crate) method_anc: MethodAncestral,
  pub(crate) dense: Option<bool>,
  pub(crate) gap_fill_args: GapFillArgs,
  pub(crate) zero_based: bool,
  pub(crate) impute_missing_data: bool,
  pub(crate) ignore_missing_alns: bool,
  pub(crate) gtr_iterations: usize,
  pub(crate) site_specific_gtr: bool,
  pub(crate) sample_from_profile: SampleMode,
  pub(crate) seed_args: SeedArgs,
  pub(crate) constant_sites: usize,
  pub(crate) rescale: f64,
  pub(crate) detailed: bool,
  pub(crate) drms: Option<PathBuf>,
  pub(crate) num_mut: usize,
  pub(crate) output: OutputCoreArgs,
  pub(crate) output_homoplasy_stats: Option<PathBuf>,
  pub(crate) output_homoplasy_report: Option<PathBuf>,
  pub(crate) output_selection: Vec<HomoplasyOutputSelection>,
  pub(crate) topology_order: TopologyOrderArgs,
}

impl TryFrom<TreetimeHomoplasyArgsRaw> for TreetimeHomoplasyArgs {
  type Error = Report;

  fn try_from(raw: TreetimeHomoplasyArgsRaw) -> Result<Self, Report> {
    let tree = raw
      .tree
      .ok_or_else(|| missing_required_args::<TreetimeHomoplasyArgsRaw>(&["tree"]))?;
    if !raw.rescale.is_finite() || raw.rescale <= 0.0 {
      return make_error!("--rescale must be a positive finite number, but got {}", raw.rescale);
    }
    Ok(Self {
      alignment: raw.alignment,
      tree,
      alphabet_args: raw.alphabet_args,
      model_args: raw.model_args,
      method_anc: raw.method_anc.into(),
      dense: raw.dense,
      gap_fill_args: raw.gap_fill_args,
      zero_based: raw.zero_based,
      impute_missing_data: raw.impute_missing_data,
      ignore_missing_alns: raw.ignore_missing_alns,
      gtr_iterations: raw.gtr_iterations,
      site_specific_gtr: raw.site_specific_gtr,
      sample_from_profile: raw.sample_from_profile.into(),
      seed_args: raw.seed_args,
      constant_sites: raw.constant_sites,
      rescale: raw.rescale,
      detailed: raw.detailed,
      drms: raw.drms,
      num_mut: raw.num_mut,
      output: raw.output,
      output_homoplasy_stats: raw.output_homoplasy_stats,
      output_homoplasy_report: raw.output_homoplasy_report,
      output_selection: raw.output_selection,
      topology_order: raw.topology_order,
    })
  }
}

#[skip_serializing_none]
#[derive(Debug, Clone, SmartDefault, Serialize, Deserialize, JsonSchema, deser::Serialize, deser::Deserialize)]
#[deser(skip_serializing_optionals)]
#[serde(default, deny_unknown_fields)]
#[deser(default, deny_unknown_fields)]
#[cfg_attr(feature = "clap", derive(clap::Parser))]
#[schemars(rename = "HomoplasyConfig")]
pub struct TreetimeHomoplasyArgsRaw {
  #[cfg_attr(feature = "clap", clap(flatten))]
  #[serde(skip)]
  #[deser(skip)]
  pub config_args: ConfigArgs,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[serde(flatten)]
  #[deser(flatten)]
  #[schemars(extend("x-path" = "input"))]
  pub alignment: AlignmentArgs,

  /// Tree in Newick format.
  #[cfg_attr(feature = "clap", clap(long, short = 't', help_heading = "Input data"))]
  #[cfg_attr(feature = "clap", clap(value_hint = ValueHint::FilePath))]
  #[schemars(extend("x-path" = "input"))]
  pub tree: Option<PathBuf>,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[serde(flatten)]
  #[deser(flatten)]
  pub alphabet_args: AlphabetArgs,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[serde(flatten)]
  #[deser(flatten)]
  pub model_args: ModelArgs,

  /// Method used for reconstructing ancestral sequences, which places the mutations on the branches
  #[cfg_attr(feature = "clap", clap(long, value_enum, default_value_t = MethodAncestralCli::default(), help_heading = "Ancestral reconstruction"))]
  pub method_anc: MethodAncestralCli,

  /// Use dense representation (stores full probability vectors at each position)
  ///
  /// When combined with `--model infer`, marginal reconstruction runs twice: once to populate
  /// profiles for GTR inference, and again with the inferred GTR.
  #[cfg_attr(feature = "clap", clap(long, help_heading = "Ancestral reconstruction"))]
  pub dense: Option<bool>,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[serde(flatten)]
  #[deser(flatten)]
  pub gap_fill_args: GapFillArgs,

  /// Report sequence positions counted from 0 instead of 1.
  ///
  /// Applies to the positions of the report, the statistics JSON, and the `GENOMIC_POSITION` column of
  /// `--drms`, which is always counted from 1.
  #[cfg_attr(feature = "clap", clap(long, help_heading = "Input data"))]
  pub zero_based: bool,

  /// Resolve ambiguous and unknown tip states (`N` and IUPAC codes such as `R`) to the most likely
  /// inferred state.
  ///
  /// Changes involving ambiguous characters on terminal branches then disappear from the report.
  /// Only defined for marginal reconstruction; a no-op with a warning under `--method-anc=parsimony`.
  #[cfg_attr(feature = "clap", clap(long, help_heading = "Ancestral reconstruction"))]
  pub impute_missing_data: bool,

  /// Treat tree tips that have no sequence in the alignment as fully ambiguous (missing data)
  /// instead of aborting.
  ///
  /// Without this flag the run aborts when more than one third of the tips lack a sequence, matching
  /// TreeTime v0.
  #[cfg_attr(feature = "clap", clap(long, help_heading = "Input data"))]
  pub ignore_missing_alns: bool,

  /// Number of outer GTR refinement iterations.
  ///
  /// Re-estimates the rate matrix from marginal posterior profiles after each
  /// reconstruction pass. Only effective with `--model infer`. Default 0 preserves
  /// the current single-pass behavior.
  #[cfg_attr(
    feature = "clap",
    clap(long, default_value_t = 0, help_heading = "Substitution model")
  )]
  pub gtr_iterations: usize,

  /// Use site-specific GTR model with per-site equilibrium frequencies.
  ///
  /// Requires `--model infer` and `--dense true`. Incompatible with sequence compression
  /// (sparse representation). When enabled, each alignment position gets its own
  /// eigendecomposition based on position-specific base composition.
  #[cfg_attr(feature = "clap", clap(long, hide = true, help_heading = "Substitution model"))]
  pub site_specific_gtr: bool,

  /// How to pick ancestral states from the marginal posterior profile.
  ///
  /// 'argmax': most likely state at every node (deterministic, default).
  /// 'root': sample from the posterior at the root only, argmax elsewhere.
  /// 'all': sample from the posterior at every node.
  ///
  /// Only affects marginal reconstruction (`--method-anc=marginal`). Use `--seed` for reproducible
  /// draws.
  #[cfg_attr(feature = "clap", clap(long, value_enum, default_value_t = SampleModeCli::default(), help_heading = "Ancestral reconstruction"))]
  pub sample_from_profile: SampleModeCli,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[serde(flatten)]
  #[deser(flatten)]
  pub seed_args: SeedArgs,

  /// Number of constant sites that the alignment leaves out.
  ///
  /// Added to the alignment length to give the number of sites of the genome, which sets the
  /// number of sites without mutations and the rate of the Poisson comparison.
  #[cfg_attr(
    feature = "clap",
    clap(long = "const", default_value_t = 0, help_heading = "Homoplasy")
  )]
  pub constant_sites: usize,

  /// Factor that multiplies every branch length of the input tree.
  ///
  /// Use it when the tree is not in substitutions per site, for example `--rescale=0.001` for a
  /// tree in substitutions per thousand sites. Scales the reconstruction, the tree outputs, and the
  /// reported tree lengths.
  #[cfg_attr(feature = "clap", clap(long, default_value_t = 1.0, help_heading = "Homoplasy"))]
  #[default = 1.0]
  pub rescale: f64,

  /// Add the mutations on terminal branches and the taxa that carry recurrent mutations to the
  /// report.
  ///
  /// The statistics JSON always contains them.
  #[cfg_attr(feature = "clap", clap(long, help_heading = "Homoplasy"))]
  pub detailed: bool,

  /// TSV file of drug resistance mutations (DRM) that annotates the report.
  ///
  /// Columns: `GENOMIC_POSITION` (counted from 1), `ALT_BASE`, `DRUG`, `GENE`, `SUBSTITUTION`.
  /// A mutation at a listed position gets the gene and the drug, and the substitution when its
  /// derived base is a listed `ALT_BASE`.
  #[cfg_attr(feature = "clap", clap(long, help_heading = "Homoplasy"))]
  #[cfg_attr(feature = "clap", clap(value_hint = ValueHint::FilePath))]
  #[schemars(extend("x-path" = "input"))]
  pub drms: Option<PathBuf>,

  /// Number of rows of each list in the report.
  ///
  /// The statistics JSON always contains the complete lists.
  #[cfg_attr(
    feature = "clap",
    clap(long, short = 'n', default_value_t = 10, help_heading = "Homoplasy")
  )]
  #[default = 10]
  pub num_mut: usize,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[serde(flatten)]
  #[deser(flatten)]
  pub output: OutputCoreArgs,

  /// Path to output homoplasy statistics JSON.
  ///
  /// Contains the counts, histograms, Poisson comparison, and complete ranked lists of the report,
  /// for substitutions, changes involving ambiguous characters, and insertions and deletions.
  ///
  /// Takes precedence over paths configured with `--output-all` and `--output-selection`.
  #[cfg_attr(feature = "clap", clap(long, value_hint = ValueHint::FilePath, help_heading = "Output"))]
  #[schemars(extend("x-path" = "output"))]
  pub output_homoplasy_stats: Option<PathBuf>,

  /// Path to output homoplasy report text.
  ///
  /// The run also logs the report at info level (`-v`).
  ///
  /// Takes precedence over paths configured with `--output-all` and `--output-selection`.
  #[cfg_attr(feature = "clap", clap(long, value_hint = ValueHint::FilePath, help_heading = "Output"))]
  #[schemars(extend("x-path" = "output"))]
  pub output_homoplasy_report: Option<PathBuf>,

  /// Comma-separated list of outputs to produce with `--output-all`.
  ///
  /// Restricts which outputs `--output-all` writes. Special value `all` expands to every output
  /// available for this command. Requires `--output-all`. Per-file flags are always honored
  /// regardless of this selection.
  #[cfg_attr(
    feature = "clap",
    clap(long, value_delimiter = ',', requires = "output_all", help_heading = "Output")
  )]
  pub output_selection: Vec<HomoplasyOutputSelection>,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[serde(flatten)]
  #[deser(flatten)]
  pub topology_order: TopologyOrderArgs,
}
