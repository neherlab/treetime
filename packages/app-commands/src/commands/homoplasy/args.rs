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
use crate::commands::shared::tree_input::TreeDialectArgs;
#[cfg(feature = "clap")]
use clap::ValueHint;
use derive_more::{Deref, DerefMut};
use deser::{Deserialize, Serialize};
use eyre::Report;
use schemars::JsonSchema;
use smart_default::SmartDefault;
use std::path::PathBuf;
use treetime::ancestral::params::{AncestralParams, MethodAncestral};
use treetime::make_error;
use treetime::partition::marginal::sample::SampleMode;
use treetime_schema::{schema_defaults, skip_serializing_optionals};

pub fn homoplasy_ancestral_params(args: &TreetimeHomoplasyArgs, seed: u64) -> AncestralParams {
  AncestralParams {
    method: args.method_anc,
    model: args.model_args.model_name(),
    dense: args.dense,
    include_leaves: false,
    report_ambiguous: true,
    impute_missing_data: args.impute_missing_data,
    site_specific_gtr: args.site_specific_gtr,
    seed,
    sample_from_profile: args.sample_from_profile,
  }
}

#[derive(Debug, Clone, Deref, DerefMut)]
pub struct TreetimeHomoplasyArgs {
  #[deref]
  #[deref_mut]
  pub(crate) raw: TreetimeHomoplasyArgsRaw,
  pub(crate) tree: PathBuf,
  pub(crate) method_anc: MethodAncestral,
  pub(crate) sample_from_profile: SampleMode,
}

impl TryFrom<TreetimeHomoplasyArgsRaw> for TreetimeHomoplasyArgs {
  type Error = Report;

  fn try_from(raw: TreetimeHomoplasyArgsRaw) -> Result<Self, Report> {
    let tree = raw
      .tree
      .clone()
      .ok_or_else(|| missing_required_args::<TreetimeHomoplasyArgsRaw>(&["tree"]))?;
    if !raw.rescale.is_finite() || raw.rescale <= 0.0 {
      return make_error!("--rescale must be a positive finite number, but got {}", raw.rescale);
    }
    let n_alignments = raw.alignment.alignment.len();
    if n_alignments > 1 {
      return make_error!(
        "homoplasy reads exactly one --alignment file, but {n_alignments} were given: \
         the reported positions refer to the coordinates of one alignment"
      );
    }
    Ok(Self {
      tree,
      method_anc: raw.method_anc.into(),
      sample_from_profile: raw.sample_from_profile.into(),
      raw,
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
#[schemars(rename = "HomoplasyConfig")]
pub struct TreetimeHomoplasyArgsRaw {
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

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub alphabet_args: AlphabetArgs,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub model_args: ModelArgs,

  /// Method used for reconstructing ancestral sequences, which places the mutations on the branches
  #[cfg_attr(feature = "clap", clap(long, value_enum, default_value_t = MethodAncestralCli::default(), help_heading = "Ancestral reconstruction"))]
  pub method_anc: MethodAncestralCli,

  /// Use dense representation (stores full probability vectors at each position)
  #[cfg_attr(feature = "clap", clap(long, help_heading = "Ancestral reconstruction"))]
  pub dense: Option<bool>,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
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
  #[schemars(flatten)]
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
  #[schemars(flatten)]
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
  #[schemars(flatten)]
  #[deser(flatten)]
  pub topology_order: TopologyOrderArgs,
}
