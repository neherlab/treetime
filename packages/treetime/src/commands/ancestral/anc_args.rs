use crate::alphabet::alphabet::AlphabetName;
use crate::gtr::get_gtr::GtrModelName;
use clap::{ArgEnum, Parser, ValueHint};
use std::fmt::Debug;
use std::path::PathBuf;

#[derive(Copy, Debug, Clone, PartialEq, Eq, PartialOrd, Ord, ArgEnum)]
#[clap(rename = "kebab-case")]
pub enum MethodAncestral {
  MaximumLikelihoodJoint,
  MaximumLikelihoodMarginal,
  Parsimony,
}

impl Default for MethodAncestral {
  fn default() -> Self {
    Self::MaximumLikelihoodJoint
  }
}

#[derive(Parser, Debug)]
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
  #[clap(value_hint = ValueHint::FilePath)]
  #[clap(display_order = 1)]
  pub input_fastas: Vec<PathBuf>,

  /// REMOVED. Use positional arguments instead.
  ///
  /// Example: treetime ancestral seq1.fasta seq2.fasta
  #[clap(long, visible_alias("aln"))]
  #[clap(value_hint = ValueHint::FilePath)]
  #[clap(hide_long_help = true, hide_short_help = true)]
  pub aln: Option<PathBuf>,

  /// FASTA file of the sequence the VCF was mapped to (only for vcf input)
  #[clap(long, short = 'r')]
  #[clap(value_hint = ValueHint::FilePath)]
  pub vcf_reference: Option<PathBuf>,

  /// Name of file containing the tree in newick, nexus, or phylip format.
  ///
  /// If none is provided, treetime will attempt to build a tree from the alignment using fasttree, iqtree, or raxml (assuming they are installed)
  #[clap(long, short = 't')]
  #[clap(value_hint = ValueHint::FilePath)]
  pub tree: PathBuf,

  /// Alphabet
  ///
  #[clap(long, short = 'a', arg_enum)]
  pub alphabet: Option<AlphabetName>,

  /// GTR model to use
  ///
  /// '--gtr infer' will infer a model from the data. Alternatively, specify the model type. If the specified model requires additional options, use '--gtr-params' to specify those.
  #[clap(long = "model", short = 'g', arg_enum, default_value_t = GtrModelName::Infer)]
  pub model_name: GtrModelName,

  /// GTR parameters for the model specified by the --gtr argument. The parameters should be feed as 'key=value' list of parameters.
  ///
  /// Example: '--gtr K80 --gtr-params kappa=0.2 pis=0.25,0.25,0.25,0.25'.
  ///
  /// See the exact definitions of the parameters in the GTR creation methods in treetime/nuc_models.py or treetime/aa_models.py
  #[clap(long)]
  pub gtr_params: Vec<String>,

  /// Method used for reconstructing ancestral sequences
  #[clap(long, arg_enum, default_value_t = MethodAncestral::default())]
  pub method_anc: MethodAncestral,

  /// Use dense representation
  ///
  /// TODO: explain this better
  #[clap(long)]
  pub dense: Option<bool>,

  /// Use aminoacid alphabet
  #[clap(long)]
  pub aa: bool,

  /// Do not fill terminal gaps
  #[clap(long)]
  pub keep_overhangs: bool,

  /// Zero-based mutation indexing
  #[clap(long)]
  pub zero_based: bool,

  /// Overwrite ambiguous states on tips with the most likely inferred state
  #[clap(long)]
  pub reconstruct_tip_states: bool,

  /// Include transitions involving ambiguous states
  #[clap(long)]
  pub report_ambiguous: bool,

  /// Directory to write the output to
  #[clap(long, short = 'O')]
  pub outdir: PathBuf,

  /// Random seed
  #[clap(long)]
  pub seed: Option<u64>,
}
