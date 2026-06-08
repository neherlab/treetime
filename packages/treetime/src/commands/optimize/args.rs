use crate::commands::shared::alignment::AlignmentArgs;
use crate::commands::shared::alphabet::AlphabetArgs;
use crate::commands::shared::gap_fill::GapFillArgs;
use crate::commands::shared::model::ModelArgs;
use crate::commands::shared::output::{DivergenceUnits, OutputArgs};
use crate::optimize::params::{BranchOptMethod, InitialGuessMode};
#[cfg(feature = "clap")]
use clap::ValueHint;
use serde::{Deserialize, Serialize};
use smart_default::SmartDefault;
use std::fmt::Debug;
use std::path::PathBuf;

#[derive(Debug, SmartDefault, Serialize, Deserialize)]
#[serde(default)]
#[cfg_attr(feature = "clap", derive(clap::Parser))]
pub struct TreetimeOptimizeArgs {
  #[cfg_attr(feature = "clap", clap(flatten))]
  pub alignment: AlignmentArgs,

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

  /// Use dense representation of sequences on the tree
  ///
  /// Dense mode stores full probability vectors at every alignment position for each
  /// node. Sparse mode stores only variable positions. Dense is more accurate when
  /// branches are long and many sites change, but uses more memory.
  #[cfg_attr(feature = "clap", clap(long))]
  pub dense: Option<bool>,

  #[cfg_attr(feature = "clap", clap(flatten))]
  pub output: OutputArgs,

  /// Units for divergence values in augur node data JSON output.
  ///
  /// `mutations-per-site` (default): branch divergence as substitutions per site.
  /// `mutations`: absolute count of reconstructed substitutions per branch,
  /// excluding ambiguous and gap positions.
  #[cfg_attr(feature = "clap", clap(long, value_enum, default_value_t = DivergenceUnits::default()))]
  pub divergence_units: DivergenceUnits,

  /// Write augur-compatible node data JSON to this path
  ///
  /// Contains per-node optimized branch lengths (divergence, substitutions per
  /// site) and the input alignment and tree paths. The output is compatible with
  /// augur export v2 --node-data, equivalent to `augur refine` run without
  /// `--timetree` (no clock or date fields). Defaults to
  /// `<outdir>/optimize.augur-node-data.json`.
  #[cfg_attr(feature = "clap", clap(long))]
  #[cfg_attr(feature = "clap", clap(value_hint = ValueHint::FilePath))]
  pub output_augur_node_data: Option<PathBuf>,

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

  #[cfg_attr(feature = "clap", clap(flatten))]
  pub gap_fill_args: GapFillArgs,
}
