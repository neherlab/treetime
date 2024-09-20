use crate::alphabet::alphabet::AlphabetName;
use crate::gtr::get_gtr::GtrModelName;
use clap::{Parser, ValueHint};
use std::fmt::Debug;
use std::path::PathBuf;

#[derive(Parser, Debug)]
pub struct TreetimeOptimizeArgs {
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

  /// Use dense representation
  ///
  /// TODO: explain this better
  #[clap(long)]
  pub dense: Option<bool>,

  /// Directory to write the output to
  #[clap(long, short = 'O')]
  pub outdir: PathBuf,

  /// Maximum number of iterations
  #[clap(long, default_value_t = 20)]
  pub max_iter: usize,

  /// Small allowable difference in the likelihood between iterations to determine if the loop should terminate
  #[clap(long, default_value_t = 1e-3)]
  pub dp: f64,
}
