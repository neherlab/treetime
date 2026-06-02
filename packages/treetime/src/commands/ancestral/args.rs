use crate::ancestral::params::MethodAncestral;
use crate::ancestral::sample::SampleMode;
use crate::commands::shared::alignment::AlignmentArgs;
use crate::commands::shared::alphabet::AlphabetArgs;
use crate::commands::shared::gap_fill::GapFillArgs;
use crate::commands::shared::model::ModelArgs;
use crate::commands::shared::output::OutputArgs;
#[cfg(feature = "clap")]
use clap::ValueHint;
use serde::{Deserialize, Serialize};
use smart_default::SmartDefault;
use std::fmt::Debug;
use std::path::PathBuf;

#[derive(Debug, SmartDefault, Serialize, Deserialize)]
#[serde(default)]
#[cfg_attr(feature = "clap", derive(clap::Parser))]
pub struct TreetimeAncestralArgs {
  #[cfg_attr(feature = "clap", clap(flatten))]
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
  pub tree: PathBuf,

  #[cfg_attr(feature = "clap", clap(flatten))]
  pub alphabet_args: AlphabetArgs,

  #[cfg_attr(feature = "clap", clap(flatten))]
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

  /// Write augur-compatible node data JSON to this path.
  ///
  /// Contains per-node nucleotide mutations, reconstructed sequences, the alignment
  /// mask, genome annotations, and the reference (root) sequence. The output is
  /// compatible with augur export v2 --node-data for Nextstrain pipeline integration.
  #[cfg_attr(feature = "clap", clap(long))]
  #[cfg_attr(feature = "clap", clap(value_hint = ValueHint::FilePath))]
  pub output_augur_node_data: Option<PathBuf>,

  /// Path template for per-gene amino-acid FASTA alignments.
  ///
  /// The template must contain `{cds}`, which is replaced with each value from
  /// `--genes`. This matches Nextclade's `--output-translations` template.
  #[cfg_attr(feature = "clap", clap(long))]
  #[cfg_attr(feature = "clap", clap(value_hint = ValueHint::FilePath))]
  pub translations: Option<String>,

  /// Gene/CDS names to reconstruct from `--translations`.
  #[cfg_attr(feature = "clap", clap(long = "genes", value_name = "GENE"))]
  pub genes: Vec<String>,

  /// GFF3 file with CDS coordinates for Augur node data annotations.
  #[cfg_attr(feature = "clap", clap(long))]
  #[cfg_attr(feature = "clap", clap(value_hint = ValueHint::FilePath))]
  pub annotation_gff: Option<PathBuf>,

  /// FASTA file with one amino-acid root/reference sequence per gene.
  #[cfg_attr(feature = "clap", clap(long))]
  #[cfg_attr(feature = "clap", clap(value_hint = ValueHint::FilePath))]
  pub aa_root_sequence: Option<PathBuf>,

  #[cfg_attr(feature = "clap", clap(flatten))]
  pub output: OutputArgs,

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
  #[cfg_attr(feature = "clap", clap(long))]
  pub seed: Option<u64>,

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
