use crate::alphabet::alphabet::AlphabetName;
use crate::gtr::get_gtr::GtrModelName;
use clap::{Parser, ValueEnum, ValueHint};
use serde::{Deserialize, Serialize};
use smart_default::SmartDefault;
use std::fmt::Debug;
use std::path::PathBuf;

/// Per-edge branch length optimization method.
///
/// Controls how `run_optimize_mixed()` finds the maximum-likelihood branch
/// length for each edge. Two orthogonal axes: algorithm (Newton-Raphson
/// vs Brent's method) and parameterization ($t$, $\sqrt{t}$, $\ln(t)$).
#[derive(Copy, Clone, Debug, PartialEq, Eq, ValueEnum, SmartDefault, Serialize, Deserialize)]
#[value(rename_all = "kebab-case")]
pub enum BranchOptMethod {
  /// Brent's method in $t$ space (derivative-free, bracket-based).
  ///
  /// Finds the maximum within a bracket derived from the grid search bounds.
  /// Convergence is independent of Hessian conditioning. Uses `argmin::BrentOpt`.
  /// Included for completeness; `brent-sqrt` dominates for convergence speed.
  Brent,

  /// Brent's method in $\sqrt{t}$ space.
  ///
  /// Matches v0 exactly (same algorithm, same parameterization). The $\sqrt{t}$
  /// reparameterization smooths the objective, giving parabolic interpolation
  /// a better fit. Default method for golden master comparison against v0.
  #[default]
  BrentSqrt,

  /// Brent's method in $\ln(t)$ space.
  ///
  /// Smoothest objective of all parameterizations, giving the best parabolic
  /// interpolation. Requires a finite lower bound in log-space.
  BrentLog,

  /// Newton-Raphson in $t$ space.
  ///
  /// Baseline Newton method matching RAxML-NG/IQ-TREE. The Poisson indel
  /// Hessian ($-k/t^2$) can dominate the substitution Hessian on short
  /// branches, causing the step-size convergence criterion to fire before
  /// the combined gradient reaches zero.
  Newton,

  /// Newton-Raphson in $\sqrt{t}$ space.
  ///
  /// Reparameterizes the optimization variable as $s = \sqrt{t}$ and applies
  /// the chain rule to transform derivatives. Reduces the indel Hessian
  /// singularity from $O(1/t^2)$ to $O(1/t)$. Residual dominance on extreme
  /// cases ($t < 0.001$, $k > 10$).
  NewtonSqrt,

  /// Newton-Raphson in $\ln(t)$ space.
  ///
  /// Eliminates the indel singularity entirely ($\ell''_{\text{indel}} = -\mu t$,
  /// bounded). Natural relative tolerance. Best conditioning of all Newton
  /// variants.
  NewtonLog,
}

/// Controls the initial branch length estimate that runs before Newton
/// optimization.
///
/// The estimate computes `#substitutions / effective_alignment_length` per
/// edge from the marginal reconstruction. When input trees already carry
/// well-calibrated branch lengths (e.g. from RAxML, IQ-TREE, or a previous
/// TreeTime run), preserving those values lets Newton converge from a
/// better starting position.
#[derive(Copy, Debug, Clone, PartialEq, Eq, PartialOrd, Ord, ValueEnum, Default, Serialize, Deserialize)]
#[value(rename_all = "kebab-case")]
pub enum InitialGuessMode {
  /// Estimate only edges with missing or invalid branch lengths, preserve
  /// valid input values. No-op when all edges have finite branch lengths.
  #[default]
  Auto,
  /// Estimate all edges, overwriting input branch lengths.
  Always,
  /// Use input branch lengths as-is. Fails if any edge has a missing or
  /// invalid branch length.
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

  /// Initial branch length estimate before Newton optimization.
  ///
  /// - auto: estimate only edges with missing or invalid branch lengths,
  ///   preserve valid input values (default)
  /// - always: estimate all edges, overwriting input branch lengths
  /// - never: use input branch lengths as-is; fails if any are missing
  #[clap(long = "branch-length-initial-guess", value_enum, default_value_t = InitialGuessMode::Auto)]
  pub branch_length_initial_guess: InitialGuessMode,

  /// Per-edge branch length optimization method.
  ///
  /// Algorithm x parameterization:
  /// - brent: Brent's method in t space (derivative-free)
  /// - brent-sqrt: Brent's method in sqrt(t) space (default, matches v0)
  /// - brent-log: Brent's method in ln(t) space
  /// - newton: Newton-Raphson in t space
  /// - newton-sqrt: Newton-Raphson in sqrt(t) space
  /// - newton-log: Newton-Raphson in ln(t) space (not yet implemented)
  #[clap(long = "opt-method", value_enum, default_value_t = BranchOptMethod::default())]
  pub opt_method: BranchOptMethod,
}
