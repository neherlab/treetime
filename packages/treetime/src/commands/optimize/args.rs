use crate::alphabet::alphabet::AlphabetName;
use crate::gtr::get_gtr::GtrModelName;
use crate::optimize::params::{BranchOptMethod, InitialGuessMode};
use crate::seq::gap_fill::GapFill;
use serde::{Deserialize, Serialize};
use smart_default::SmartDefault;
use std::fmt::Debug;
use std::path::PathBuf;
#[cfg(feature = "clap")]
use clap::ValueHint;

#[derive(Debug, SmartDefault, Serialize, Deserialize)]
#[serde(default)]
#[cfg_attr(feature = "clap", derive(clap::Parser))]
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
  #[cfg_attr(feature = "clap", clap(long = "aln", value_hint = ValueHint::FilePath, value_name = "FILEPATH"))]
  pub input_fastas: Vec<PathBuf>,

  /// Name of file containing the tree in newick, nexus, or phylip format.
  ///
  /// If none is provided, treetime will attempt to build a tree from the alignment using fasttree, iqtree, or raxml (assuming they are installed)
  #[cfg_attr(feature = "clap", clap(long, short = 't'))]
  #[cfg_attr(feature = "clap", clap(value_hint = ValueHint::FilePath))]
  pub tree: PathBuf,

  /// Alphabet
  ///
  #[cfg_attr(feature = "clap", clap(long, short = 'a', value_enum))]
  pub alphabet: Option<AlphabetName>,

  /// GTR model to use
  ///
  /// '--model infer' will infer a model from the data. Alternatively, specify the model
  /// type. If the specified model requires additional options, use '--gtr-params' to
  /// specify those.
  #[cfg_attr(feature = "clap", clap(long = "model", short = 'g', value_enum, default_value_t = GtrModelName::Infer))]
  pub model_name: GtrModelName,

  /// Use dense representation of sequences on the tree
  ///
  /// Dense mode stores full probability vectors at every alignment position for each
  /// node. Sparse mode stores only variable positions. Dense is more accurate when
  /// branches are long and many sites change, but uses more memory.
  #[cfg_attr(feature = "clap", clap(long))]
  pub dense: Option<bool>,

  /// Directory to write the output to
  #[cfg_attr(feature = "clap", clap(long, short = 'O'))]
  pub outdir: PathBuf,

  /// Maximum number of iterations
  #[cfg_attr(feature = "clap", clap(long, default_value_t = 10))]
  #[default = 10]
  pub max_iter: usize,

  /// Likelihood convergence threshold. The loop stops when successive
  /// likelihoods differ by less than this value, or when a 2-cycle with
  /// amplitude below this value is detected.
  #[cfg_attr(feature = "clap", clap(long, default_value_t = 0.1))]
  #[default = 0.1]
  pub dp: f64,

  /// Damping factor for outer-loop branch length updates.
  ///
  /// Controls how aggressively new branch lengths replace old ones during
  /// iterative optimization. At each iteration i, the update is:
  ///   bl = bl_new * (1 - d) + bl_old * d
  /// where d = max(damping^(i+1), 0.01). The 1% floor prevents fully
  /// undamped late iterations on non-monotone objectives.
  ///
  /// Higher values are more conservative (slower convergence, less oscillation).
  /// Set to 0.0 to disable damping (full update each iteration, no floor).
  /// Must be in [0.0, 1.0).
  #[cfg_attr(feature = "clap", clap(long, default_value_t = 0.75))]
  #[default = 0.75]
  pub damping: f64,

  /// Initial branch length estimate before Newton optimization.
  ///
  /// - auto: estimate only edges with missing or invalid branch lengths,
  ///   preserve valid input values (default)
  /// - always: estimate all edges, overwriting input branch lengths
  /// - never: use input branch lengths as-is; fails if any are missing
  #[cfg_attr(feature = "clap", clap(long = "branch-length-initial-guess", value_enum, default_value_t = InitialGuessMode::Auto))]
  pub branch_length_initial_guess: InitialGuessMode,

  /// Per-edge branch length optimization method.
  ///
  /// Algorithm x parameterization:
  /// - brent: Brent's method in t space (derivative-free)
  /// - brent-sqrt: Brent's method in sqrt(t) space (default, matches v0)
  /// - brent-log: Brent's method in ln(t) space
  /// - newton: Newton-Raphson in t space
  /// - newton-sqrt: Newton-Raphson in sqrt(t) space
  /// - newton-log: Newton-Raphson in ln(t) space
  #[cfg_attr(feature = "clap", clap(long = "opt-method", value_enum, default_value_t = BranchOptMethod::default()))]
  pub opt_method: BranchOptMethod,

  /// Disable indel (insertion/deletion) contributions to branch-length
  /// optimization.
  ///
  /// When set, the optimizer uses substitution-only likelihood, matching
  /// standard phylogenetic tools (RAxML, IQ-TREE, PhyML, BEAST) and
  /// enabling v0 parity testing. Default: indels enabled.
  #[cfg_attr(feature = "clap", clap(long))]
  pub no_indels: bool,

  /// Gap fill strategy for terminal and internal gaps
  #[cfg_attr(feature = "clap", clap(long, value_enum, default_value_t = GapFill::default(), conflicts_with = "keep_overhangs"))]
  pub gap_fill: GapFill,

  /// Do not fill terminal gaps (deprecated: use --gap-fill=none)
  #[cfg_attr(feature = "clap", clap(long, hide = true))]
  pub keep_overhangs: bool,
}

impl TreetimeOptimizeArgs {
  pub fn effective_gap_fill(&self) -> GapFill {
    if self.keep_overhangs {
      GapFill::None
    } else {
      self.gap_fill
    }
  }
}
