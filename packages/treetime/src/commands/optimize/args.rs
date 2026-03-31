use crate::alphabet::alphabet::AlphabetName;
use crate::gtr::get_gtr::GtrModelName;
use clap::{Parser, ValueEnum, ValueHint};
use serde::{Deserialize, Serialize};
use std::fmt::Debug;
use std::path::PathBuf;

/// Controls when the discrete-count initial branch length estimate runs
/// before Newton optimization.
///
/// The initial guess computes `#substitutions / effective_length` per edge,
/// which is a crude starting point. When input trees already carry
/// well-calibrated branch lengths (e.g. from RAxML, IQ-TREE, or a previous
/// TreeTime run), skipping the initial guess preserves those values and
/// lets Newton converge from a better starting position.
#[derive(Copy, Debug, Clone, PartialEq, Eq, PartialOrd, Ord, ValueEnum, Serialize, Deserialize)]
#[value(rename_all = "kebab-case")]
#[derive(Default)]
pub enum InitialGuessMode {
  /// Run only when any edge lacks a branch length
  #[default]
  Auto,
  /// Run unconditionally (overwrites input branch lengths)
  Always,
  /// Skip entirely (use input branch lengths as-is)
  Never,
}

#[derive(Parser, Debug, Serialize)]
pub struct TreetimeOptimizeArgs {
  /// Path to one or multiple FASTA files with aligned input sequences
  ///
  /// Accepts plain or compressed FASTA files. If a compressed fasta file is provided, it will be transparently
  /// decompressed. Supported compression formats: `gz`, `bz2`, `xz`, `zstd`. Decompressor is chosen based on file
  /// extension. If there's multiple input files, then different files can have different compression formats.
  ///
  /// Use '-' to read uncompressed FASTA from standard input (stdin).
  ///
  /// See: https://en.wikipedia.org/wiki/FASTA_format
  #[clap(long = "aln", value_hint = ValueHint::FilePath, value_name = "FILEPATH")]
  pub input_fastas: Vec<PathBuf>,

  /// Name of file containing the tree in newick, nexus, or phylip format.
  ///
  /// If none is provided, treetime will attempt to build a tree from the alignment using fasttree, iqtree, or raxml (assuming they are installed)
  #[clap(long, short = 't')]
  #[clap(value_hint = ValueHint::FilePath)]
  pub tree: PathBuf,

  /// Alphabet
  ///
  #[clap(long, short = 'a', value_enum)]
  pub alphabet: Option<AlphabetName>,

  /// GTR model to use
  ///
  /// '--model infer' will infer a model from the data.
  ///
  /// TODO: Alternatively, specify the model type. If the specified model requires additional options, use '--gtr-params' to specify those.
  #[clap(long = "model", short = 'g', value_enum, default_value_t = GtrModelName::Infer)]
  pub model_name: GtrModelName,

  /// Use dense representation of sequences on the tree, useful if branches are long.
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
  #[clap(long, default_value_t = 1e-2)]
  pub dp: f64,

  /// Damping factor for outer-loop branch length updates.
  ///
  /// Controls how aggressively new branch lengths replace old ones during
  /// iterative optimization. At each iteration i, the update is:
  ///   bl = bl_new * (1 - damping^(i+1)) + bl_old * damping^(i+1)
  ///
  /// Higher values are more conservative (slower convergence, less oscillation).
  /// Set to 0.0 to disable damping (full Newton update each iteration).
  /// Must be in [0.0, 1.0).
  #[clap(long, default_value_t = 0.75)]
  pub damping: f64,

  /// When to compute the discrete-count initial branch length estimate
  /// before Newton optimization.
  ///
  /// - auto:   run only when any edge has no branch length (default)
  /// - always: run unconditionally, overwriting input values
  /// - never:  skip, use input branch lengths as-is
  #[clap(long, value_enum, default_value_t = InitialGuessMode::Auto)]
  pub initial_guess: InitialGuessMode,
}
