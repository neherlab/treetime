use crate::ancestral::params::MethodAncestral;
use crate::ancestral::sample::SampleMode;
use crate::commands::ancestral::aa_model::AaModelName;
use crate::commands::shared::alignment::AlignmentArgs;
use crate::commands::shared::alphabet::AlphabetArgs;
use crate::commands::shared::config::ConfigArgs;
use crate::commands::shared::gap_fill::GapFillArgs;
use crate::commands::shared::model::ModelArgs;
use crate::commands::shared::output::{AncestralOutputSelection, OutputCoreArgs, TopologyOrderArgs};
use crate::commands::shared::required::missing_required_args;
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
pub struct TreetimeAncestralArgsRaw {
  #[cfg_attr(feature = "clap", clap(flatten))]
  #[serde(skip)]
  pub config_args: ConfigArgs,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[serde(flatten)]
  pub alignment: AlignmentArgs,

  /// FASTA file of the sequence the VCF was mapped to (only for vcf input)
  #[cfg_attr(feature = "clap", clap(long, short = 'r'))]
  #[cfg_attr(feature = "clap", clap(value_hint = ValueHint::FilePath))]
  pub vcf_reference: Option<PathBuf>,

  /// Name of file containing the tree in newick, nexus, or phylip format.
  ///
  /// If none is provided, treetime will attempt to build a tree from the alignment using fasttree, iqtree, or raxml (assuming they are installed)
  #[cfg_attr(feature = "clap", clap(long, short = 't'))]
  #[cfg_attr(feature = "clap", clap(value_hint = ValueHint::FilePath))]
  pub tree: Option<PathBuf>,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[serde(flatten)]
  pub alphabet_args: AlphabetArgs,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[serde(flatten)]
  pub model_args: ModelArgs,

  /// Method used for reconstructing ancestral sequences
  #[cfg_attr(feature = "clap", clap(long, value_enum, default_value_t = MethodAncestral::default()))]
  pub method_anc: MethodAncestral,

  /// Use dense representation (stores full probability vectors at each position)
  ///
  /// When combined with `--model infer`, marginal reconstruction runs twice: once to populate
  /// profiles for GTR inference, and again with the inferred GTR.
  #[cfg_attr(feature = "clap", clap(long))]
  pub dense: Option<bool>,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[serde(flatten)]
  pub gap_fill_args: GapFillArgs,

  /// Zero-based mutation indexing
  #[cfg_attr(feature = "clap", clap(long))]
  pub zero_based: bool,

  /// Emit reconstructed leaf (tip) sequences in addition to internal nodes.
  #[cfg_attr(feature = "clap", clap(long))]
  pub include_leaves: bool,

  /// Resolve ambiguous and unknown tip states (`N` and IUPAC codes such as `R`) to the most likely
  /// inferred state.
  ///
  /// Gaps are left as deletions (inferred structure, not missing data). Only defined for marginal
  /// reconstruction; a no-op with a warning under `--method-anc=parsimony`.
  #[cfg_attr(feature = "clap", clap(long))]
  pub impute_missing_data: bool,

  /// v0-compatible alias for `--include-leaves --impute-missing-data`.
  ///
  /// Emits tip sequences and resolves ambiguous/unknown tip states to the most likely inferred state.
  #[cfg_attr(feature = "clap", clap(long))]
  pub reconstruct_tip_states: bool,

  /// Include transitions involving ambiguous states
  #[cfg_attr(feature = "clap", clap(long))]
  pub report_ambiguous: bool,

  /// Treat tree tips that have no sequence in the alignment as fully ambiguous (missing data)
  /// instead of aborting.
  ///
  /// Without this flag the run aborts when more than one third of the tips lack a sequence, matching
  /// TreeTime v0. Useful when consuming per-CDS translations where some samples have no peptide for a
  /// given CDS.
  #[cfg_attr(feature = "clap", clap(long))]
  pub ignore_missing_alns: bool,

  /// Path to output augur-compatible node data JSON.
  ///
  /// Contains per-node nucleotide mutations, reconstructed sequences, the alignment
  /// mask, genome annotations, and the reference (root) sequence. The output is
  /// compatible with augur export v2 --node-data for Nextstrain pipeline integration.
  ///
  /// Takes precedence over paths configured with `--output-all` and `--output-selection`.
  #[cfg_attr(feature = "clap", clap(long, value_hint = ValueHint::FilePath, help_heading = "Output"))]
  pub output_augur_node_data: Option<PathBuf>,

  /// Path to output GTR model JSON.
  ///
  /// Takes precedence over paths configured with `--output-all` and `--output-selection`.
  #[cfg_attr(feature = "clap", clap(long, value_hint = ValueHint::FilePath, help_heading = "Output"))]
  pub output_gtr: Option<PathBuf>,

  /// Path to output reconstructed nucleotide FASTA.
  ///
  /// Takes precedence over paths configured with `--output-all` and `--output-selection`.
  #[cfg_attr(feature = "clap", clap(long, value_hint = ValueHint::FilePath, help_heading = "Output"))]
  pub output_reconstructed_nuc_fasta: Option<PathBuf>,

  /// Path template for per-CDS amino-acid FASTA alignments.
  ///
  /// The template must contain a CDS placeholder, replaced with each value from `--cdses` (or each
  /// CDS in `--annotation` when `--cdses` is omitted). Both `{cds}` (Nextclade
  /// `--output-translations`) and `%GENE` (augur) placeholders are accepted.
  #[cfg_attr(feature = "clap", clap(long))]
  #[cfg_attr(feature = "clap", clap(value_hint = ValueHint::FilePath))]
  pub translations: Option<String>,

  /// Comma-separated CDS names to reconstruct from `--translations`.
  ///
  /// When omitted, the CDS set is derived from `--annotation`.
  #[cfg_attr(
    feature = "clap",
    clap(long = "cdses", visible_alias = "genes", value_name = "CDS", value_delimiter = ',')
  )]
  pub cdses: Vec<String>,

  /// GFF3 file with CDS coordinates for Augur node data annotations.
  ///
  /// Also supplies the CDS set when `--cdses` is omitted.
  #[cfg_attr(feature = "clap", clap(long, alias = "annotation-gff"))]
  #[cfg_attr(feature = "clap", clap(value_hint = ValueHint::FilePath))]
  pub annotation: Option<PathBuf>,

  /// FASTA file with one amino-acid root/reference sequence per CDS.
  #[cfg_attr(feature = "clap", clap(long))]
  #[cfg_attr(feature = "clap", clap(value_hint = ValueHint::FilePath))]
  pub aa_root_sequence: Option<PathBuf>,

  /// Amino-acid substitution model. Mirrors the nucleotide `--model`; default `infer` matches augur.
  #[cfg_attr(feature = "clap", clap(long, value_enum, default_value_t = AaModelName::default()))]
  pub aa_model: AaModelName,

  /// Path template for per-CDS reconstructed amino-acid FASTA output (including internal nodes).
  ///
  /// Off by default. When set, the reconstructed sequence of every node is written per CDS. Accepts
  /// the same `{cds}`/`%GENE` placeholders as `--translations`.
  ///
  /// Takes precedence over paths configured with `--output-all` and `--output-selection`.
  #[cfg_attr(feature = "clap", clap(long, value_hint = ValueHint::FilePath, help_heading = "Output"))]
  pub output_reconstructed_aa_fasta: Option<String>,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[serde(flatten)]
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
  pub output_selection: Vec<AncestralOutputSelection>,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[serde(flatten)]
  pub topology_order: TopologyOrderArgs,

  /// Number of outer GTR refinement iterations.
  ///
  /// Re-estimates the rate matrix from marginal posterior profiles after each
  /// reconstruction pass. Only effective with `--model infer`. Default 0 preserves
  /// the current single-pass behavior. Mugration uses 5 by default.
  #[cfg_attr(feature = "clap", clap(long, default_value_t = 0))]
  pub gtr_iterations: usize,

  /// Use site-specific GTR model with per-site equilibrium frequencies.
  ///
  /// Requires `--model infer` and `--dense true`. Incompatible with sequence compression
  /// (sparse representation). When enabled, each alignment position gets its own
  /// eigendecomposition based on position-specific base composition.
  #[cfg_attr(feature = "clap", clap(long, hide = true))]
  pub site_specific_gtr: bool,

  /// Random seed
  #[cfg_attr(feature = "clap", clap(long, visible_alias = "rng-seed"))]
  pub seed: Option<u64>,

  /// Use amino-acid alphabet (v0 compat, equivalent to `--alphabet=aa`)
  #[cfg_attr(feature = "clap", clap(long, hide = true))]
  pub aa: bool,

  /// Shortcut for `--method-anc=marginal` (v0 compat)
  #[cfg_attr(feature = "clap", clap(long, hide = true))]
  pub marginal: bool,

  /// Load a custom GTR model from file (not yet implemented)
  #[cfg_attr(feature = "clap", clap(long, hide = true))]
  #[cfg_attr(feature = "clap", clap(value_hint = ValueHint::FilePath))]
  pub custom_gtr: Option<PathBuf>,

  /// How to pick ancestral states from the marginal posterior profile.
  ///
  /// 'argmax': most likely state at every node (deterministic, default).
  /// 'root': sample from the posterior at the root only, argmax elsewhere (matches augur's
  /// `sample_from_profile='root'`). Use `--seed` for reproducible draws.
  /// 'all': sample from the posterior at every node.
  ///
  /// Only affects marginal reconstruction (`--method-anc=marginal`).
  #[cfg_attr(feature = "clap", clap(long, value_enum, default_value_t = SampleMode::default()))]
  pub sample_from_profile: SampleMode,
}

/// Ancestral reconstruction arguments with required inputs proven present.
///
/// Produced from [`TreetimeAncestralArgsRaw`] by [`TryFrom`] once the `--config` overlay has run, so
/// the run code reads `tree` without an `Option`.
#[derive(Debug, Clone)]
pub struct TreetimeAncestralArgs {
  pub alignment: AlignmentArgs,
  pub vcf_reference: Option<PathBuf>,
  pub tree: PathBuf,
  pub alphabet_args: AlphabetArgs,
  pub model_args: ModelArgs,
  pub method_anc: MethodAncestral,
  pub dense: Option<bool>,
  pub gap_fill_args: GapFillArgs,
  pub zero_based: bool,
  pub include_leaves: bool,
  pub impute_missing_data: bool,
  pub reconstruct_tip_states: bool,
  pub report_ambiguous: bool,
  pub ignore_missing_alns: bool,
  pub output_augur_node_data: Option<PathBuf>,
  pub output_gtr: Option<PathBuf>,
  pub output_reconstructed_nuc_fasta: Option<PathBuf>,
  pub translations: Option<String>,
  pub cdses: Vec<String>,
  pub annotation: Option<PathBuf>,
  pub aa_root_sequence: Option<PathBuf>,
  pub aa_model: AaModelName,
  pub output_reconstructed_aa_fasta: Option<String>,
  pub output: OutputCoreArgs,
  pub output_selection: Vec<AncestralOutputSelection>,
  pub topology_order: TopologyOrderArgs,
  pub gtr_iterations: usize,
  pub site_specific_gtr: bool,
  pub seed: Option<u64>,
  pub aa: bool,
  pub marginal: bool,
  pub custom_gtr: Option<PathBuf>,
  pub sample_from_profile: SampleMode,
}

impl TreetimeAncestralArgs {
  /// Input tree path.
  pub fn tree(&self) -> &Path {
    &self.tree
  }
}

impl TryFrom<TreetimeAncestralArgsRaw> for TreetimeAncestralArgs {
  type Error = Report;

  fn try_from(raw: TreetimeAncestralArgsRaw) -> Result<Self, Report> {
    let tree = raw
      .tree
      .ok_or_else(|| missing_required_args::<TreetimeAncestralArgsRaw>(&["tree"]))?;
    Ok(Self {
      alignment: raw.alignment,
      vcf_reference: raw.vcf_reference,
      tree,
      alphabet_args: raw.alphabet_args,
      model_args: raw.model_args,
      method_anc: raw.method_anc,
      dense: raw.dense,
      gap_fill_args: raw.gap_fill_args,
      zero_based: raw.zero_based,
      include_leaves: raw.include_leaves,
      impute_missing_data: raw.impute_missing_data,
      reconstruct_tip_states: raw.reconstruct_tip_states,
      report_ambiguous: raw.report_ambiguous,
      ignore_missing_alns: raw.ignore_missing_alns,
      output_augur_node_data: raw.output_augur_node_data,
      output_gtr: raw.output_gtr,
      output_reconstructed_nuc_fasta: raw.output_reconstructed_nuc_fasta,
      translations: raw.translations,
      cdses: raw.cdses,
      annotation: raw.annotation,
      aa_root_sequence: raw.aa_root_sequence,
      aa_model: raw.aa_model,
      output_reconstructed_aa_fasta: raw.output_reconstructed_aa_fasta,
      output: raw.output,
      output_selection: raw.output_selection,
      topology_order: raw.topology_order,
      gtr_iterations: raw.gtr_iterations,
      site_specific_gtr: raw.site_specific_gtr,
      seed: raw.seed,
      aa: raw.aa,
      marginal: raw.marginal,
      custom_gtr: raw.custom_gtr,
      sample_from_profile: raw.sample_from_profile,
    })
  }
}
