use crate::alphabet::alphabet::AlphabetName;
use crate::ancestral::params::MethodAncestral;
use crate::gtr::get_gtr::GtrModelName;
use crate::seq::gap_fill::GapFill;
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
  /// Path to one or multiple FASTA files with aligned input sequences
  ///
  /// Accepts plain or compressed FASTA files. If a compressed fasta file is provided, it will be transparently
  /// decompressed. Supported compression formats: `gz`, `bz2`, `xz`, `zstd`. Decompressor is chosen based on file
  /// extension. If there's multiple input files, then different files can have different compression formats.
  ///
  /// If no input files provided, the plain fasta input is read from standard input (stdin).
  ///
  /// See: https://en.wikipedia.org/wiki/FASTA_format
  #[cfg_attr(feature = "clap", clap(value_hint = ValueHint::FilePath))]
  #[cfg_attr(feature = "clap", clap(display_order = 1))]
  pub input_fastas: Vec<PathBuf>,

  /// REMOVED. Use positional arguments instead.
  ///
  /// Example: treetime ancestral seq1.fasta seq2.fasta
  #[cfg_attr(feature = "clap", clap(long, visible_alias("aln")))]
  #[cfg_attr(feature = "clap", clap(value_hint = ValueHint::FilePath))]
  #[cfg_attr(feature = "clap", clap(hide_long_help = true, hide_short_help = true))]
  pub aln: Option<PathBuf>,

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

  /// Alphabet
  ///
  #[cfg_attr(feature = "clap", arg(long, short = 'a', value_enum))]
  pub alphabet: Option<AlphabetName>,

  /// GTR model to use
  ///
  /// '--gtr infer' will infer a model from the data. Alternatively, specify the model type. If the specified model requires additional options, use '--gtr-params' to specify those.
  #[cfg_attr(feature = "clap", clap(long = "model", short = 'g', value_enum, default_value_t = GtrModelName::Infer))]
  pub model_name: GtrModelName,

  /// GTR parameters for the model specified by the --gtr argument. The parameters should be feed as 'key=value' list of parameters.
  ///
  /// Example: '--gtr K80 --gtr-params kappa=0.2 pis=0.25,0.25,0.25,0.25'.
  ///
  /// See the exact definitions of the parameters in the GTR creation methods in treetime/nuc_models.py or treetime/aa_models.py
  #[cfg_attr(feature = "clap", clap(long))]
  pub gtr_params: Vec<String>,

  /// Method used for reconstructing ancestral sequences
  #[cfg_attr(feature = "clap", clap(long, value_enum, default_value_t = MethodAncestral::default()))]
  pub method_anc: MethodAncestral,

  /// Use dense representation (stores full probability vectors at each position)
  ///
  /// When combined with `--model infer`, marginal reconstruction runs twice: once to populate
  /// profiles for GTR inference, and again with the inferred GTR.
  #[cfg_attr(feature = "clap", clap(long))]
  pub dense: Option<bool>,

  /// Use aminoacid alphabet
  #[cfg_attr(feature = "clap", clap(long))]
  pub aa: bool,

  /// How to handle gap characters in input sequences
  ///
  /// 'only-terminal': replace leading and trailing gap characters with the ambiguous character (default, matches v0).
  /// 'all': replace all gap characters with the ambiguous character.
  /// 'none': leave all gap characters unchanged.
  #[cfg_attr(feature = "clap", clap(long, value_enum, default_value_t = GapFill::default(), conflicts_with = "keep_overhangs"))]
  pub gap_fill: GapFill,

  /// Do not fill terminal gaps (deprecated: use --gap-fill=none)
  #[cfg_attr(feature = "clap", clap(long, hide = true))]
  pub keep_overhangs: bool,

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

  /// Directory to write the output to
  #[cfg_attr(feature = "clap", clap(long, short = 'O'))]
  pub outdir: PathBuf,

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
}

impl TreetimeAncestralArgs {
  pub fn effective_gap_fill(&self) -> GapFill {
    if self.keep_overhangs {
      GapFill::None
    } else {
      self.gap_fill
    }
  }
}
