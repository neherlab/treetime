use crate::commands::shared::alignment::AlignmentArgs;
use crate::commands::shared::alphabet::AlphabetArgs;
use crate::commands::shared::config::ConfigArgs;
use crate::commands::shared::gap_fill::GapFillArgs;
use crate::commands::shared::model::ModelArgs;
use crate::commands::shared::output_args::{DivergenceUnits, OptimizeOutputSelection, OutputCoreArgs};
use crate::commands::shared::required::missing_required_args;
use crate::commands::shared::topology_order_args::TopologyOrderArgs;
#[cfg(feature = "clap")]
use clap::ValueHint;
use eyre::Report;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use smart_default::SmartDefault;
use std::fmt::Debug;
use std::path::{Path, PathBuf};
use treetime::clock::find_best_root::params::{RerootMethod, RerootSpec};
use treetime::optimize::params::{BranchOptMethod, InitialGuessMode, TopologyOps};

/// Reroot methods available in the optimize command.
///
/// Only date-free methods are valid here because optimize has no sampling dates.
/// Date-dependent methods (least-squares, oldest, clock-filter) are available
/// in the timetree and clock commands.
#[derive(Copy, Debug, Clone, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
#[cfg_attr(feature = "clap", derive(clap::ValueEnum))]
#[serde(rename_all = "kebab-case")]
pub enum OptimizeRerootMethod {
  MinDev,
}

impl From<OptimizeRerootMethod> for RerootMethod {
  fn from(m: OptimizeRerootMethod) -> Self {
    match m {
      OptimizeRerootMethod::MinDev => RerootMethod::MinDev,
    }
  }
}

/// Per-edge branch length optimization method.
///
/// Controls how `run_optimize_mixed()` finds the maximum-likelihood branch
/// length for each edge. Two orthogonal axes: algorithm (Newton-Raphson
/// vs Brent's method) and parameterization ($t$, $\sqrt{t}$, $\ln(t)$).
#[derive(Copy, Clone, Debug, PartialEq, Eq, SmartDefault, Serialize, Deserialize, JsonSchema)]
#[cfg_attr(feature = "clap", derive(clap::ValueEnum))]
#[serde(rename_all = "kebab-case")]
#[schemars(rename = "BranchOptMethod")]
pub enum BranchOptMethodCli {
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

impl From<BranchOptMethodCli> for BranchOptMethod {
  fn from(method: BranchOptMethodCli) -> Self {
    match method {
      BranchOptMethodCli::Brent => BranchOptMethod::Brent,
      BranchOptMethodCli::BrentSqrt => BranchOptMethod::BrentSqrt,
      BranchOptMethodCli::BrentLog => BranchOptMethod::BrentLog,
      BranchOptMethodCli::Newton => BranchOptMethod::Newton,
      BranchOptMethodCli::NewtonSqrt => BranchOptMethod::NewtonSqrt,
      BranchOptMethodCli::NewtonLog => BranchOptMethod::NewtonLog,
    }
  }
}

/// Controls whether marginal reconstruction estimates initial branch lengths
/// from substitutions divided by effective alignment length. Preserving valid
/// input lengths can provide a better Newton starting point.
#[derive(Copy, Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Default, Serialize, Deserialize, JsonSchema)]
#[cfg_attr(feature = "clap", derive(clap::ValueEnum))]
#[serde(rename_all = "kebab-case")]
#[schemars(rename = "InitialGuessMode")]
pub enum InitialGuessModeCli {
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

impl From<InitialGuessModeCli> for InitialGuessMode {
  fn from(mode: InitialGuessModeCli) -> Self {
    match mode {
      InitialGuessModeCli::Auto => InitialGuessMode::Auto,
      InitialGuessModeCli::Always => InitialGuessMode::Always,
      InitialGuessModeCli::Never => InitialGuessMode::Never,
    }
  }
}

#[derive(Debug, Clone, SmartDefault, Serialize, Deserialize, JsonSchema)]
#[serde(default, deny_unknown_fields)]
#[cfg_attr(feature = "clap", derive(clap::Parser))]
pub struct TreetimeOptimizeArgsRaw {
  #[cfg_attr(feature = "clap", clap(flatten))]
  #[serde(skip)]
  pub config_args: ConfigArgs,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[serde(flatten)]
  pub alignment: AlignmentArgs,

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

  /// Use dense representation of sequences on the tree
  ///
  /// Dense mode stores full probability vectors at every alignment position for each
  /// node. Sparse mode stores only variable positions. Dense is more accurate when
  /// branches are long and many sites change, but uses more memory.
  #[cfg_attr(feature = "clap", clap(long))]
  pub dense: Option<bool>,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[serde(flatten)]
  pub output: OutputCoreArgs,

  /// Units for divergence values in augur node data JSON output.
  ///
  /// `mutations-per-site` (default): branch divergence as substitutions per site.
  /// `mutations`: absolute count of reconstructed substitutions per branch,
  /// excluding ambiguous and gap positions.
  #[cfg_attr(feature = "clap", clap(long, value_enum, default_value_t = DivergenceUnits::default()))]
  pub divergence_units: DivergenceUnits,

  /// Path to output augur-compatible node data JSON.
  ///
  /// Contains per-node optimized branch lengths (divergence, substitutions per
  /// site) and the input alignment and tree paths. The output is compatible with
  /// augur export v2 --node-data.
  ///
  /// Takes precedence over paths configured with `--output-all` and `--output-selection`.
  #[cfg_attr(feature = "clap", clap(long, value_hint = ValueHint::FilePath, help_heading = "Output"))]
  pub output_augur_node_data: Option<PathBuf>,

  /// Path to output GTR model JSON.
  ///
  /// Takes precedence over paths configured with `--output-all` and `--output-selection`.
  #[cfg_attr(feature = "clap", clap(long, value_hint = ValueHint::FilePath, help_heading = "Output"))]
  pub output_gtr: Option<PathBuf>,

  /// Comma-separated list of outputs to produce with `--output-all`.
  ///
  /// Restricts which outputs `--output-all` writes. Special value `all` expands to every output
  /// available for this command. Requires `--output-all`. Per-file flags are always honored
  /// regardless of this selection.
  #[cfg_attr(
    feature = "clap",
    clap(long, value_delimiter = ',', requires = "output_all", help_heading = "Output")
  )]
  pub output_selection: Vec<OptimizeOutputSelection>,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[serde(flatten)]
  pub topology_order: TopologyOrderArgs,

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

  /// Damping factor $d$ for outer-loop updates: $b=b_{new}(1-d)+b_{old}d$,
  /// where $d=max(damping^{i+1},0.01)$. Higher values reduce oscillation;
  /// zero disables damping. Must be in $[0,1)$.
  #[cfg_attr(feature = "clap", clap(long, default_value_t = 0.75))]
  #[default = 0.75]
  pub damping: f64,

  /// Initial branch length estimate before Newton optimization.
  ///
  /// - auto: estimate only edges with missing or invalid branch lengths,
  ///   preserve valid input values (default)
  /// - always: estimate all edges, overwriting input branch lengths
  /// - never: use input branch lengths as-is; fails if any are missing
  #[cfg_attr(feature = "clap", clap(long = "branch-length-initial-guess", value_enum, default_value_t = InitialGuessModeCli::Auto))]
  pub branch_length_initial_guess: InitialGuessModeCli,

  /// Per-edge optimizer and parameterization. Values are `brent`, `brent-sqrt`
  /// (default and v0-compatible), `brent-log`, `newton`, `newton-sqrt`, and
  /// `newton-log`; suffixes select $\sqrt{t}$ or $\ln(t)$.
  #[cfg_attr(feature = "clap", clap(long = "opt-method", value_enum, default_value_t = BranchOptMethodCli::default()))]
  pub opt_method: BranchOptMethodCli,

  /// Disable indel (insertion/deletion) contributions to branch-length
  /// optimization.
  ///
  /// When set, the optimizer uses substitution-only likelihood, matching
  /// standard phylogenetic tools (RAxML, IQ-TREE, PhyML, BEAST) and
  /// enabling v0 parity testing. Default: indels enabled.
  #[cfg_attr(feature = "clap", clap(long))]
  pub no_indels: bool,

  /// Reroot the tree by minimizing root-to-tip divergence variance.
  ///
  /// By default, optimize keeps the input root. Pass --reroot or --reroot=min-dev
  /// to enable divergence-based rerooting. Date-dependent methods (least-squares,
  /// oldest, clock-filter) are available in the timetree and clock commands.
  #[cfg_attr(feature = "clap", clap(
    long,
    value_enum,
    num_args = 0..=1,
    default_missing_value = "min-dev",
    conflicts_with = "reroot_tips",
  ))]
  pub reroot: Option<OptimizeRerootMethod>,

  /// Reroot on the branch leading to a tip or the MRCA of a comma-separated tip list.
  #[cfg_attr(feature = "clap", clap(long, value_delimiter = ',', conflicts_with = "reroot"))]
  pub reroot_tips: Vec<String>,

  /// Keep the input tree root instead of rerooting.
  ///
  /// Optimize keeps the input root by default; this flag is the explicit form and
  /// is mutually exclusive with the reroot options.
  #[cfg_attr(feature = "clap", clap(long, conflicts_with_all = ["reroot", "reroot_tips"]))]
  pub keep_root: bool,

  /// Disable collapsing of internal branches whose optimized length is zero.
  ///
  /// By default the optimize loop contracts internal edges the per-edge optimizer drove to
  /// exactly zero that carry no substitutions or indels, turning the resulting binary nodes into
  /// polytomies. When set, such edges are kept in the output tree with length zero.
  #[cfg_attr(feature = "clap", clap(long))]
  pub no_collapse_short_branches: bool,

  /// Disable merging of polytomy siblings that share substitutions.
  ///
  /// By default, sibling branches in a polytomy that carry identical substitutions are grouped
  /// under a new internal node. Requires the sparse sequence representation, so this step has no
  /// effect under `--dense` and the flag is then a no-op.
  #[cfg_attr(feature = "clap", clap(long))]
  pub no_merge_siblings: bool,

  /// Disable the flip-parent-child step (reversion hoist) in polytomies.
  ///
  /// By default, when a child branch carries the exact reversion of a substitution on the node's
  /// parent branch, a new node is inserted that groups the node with that child, removing one
  /// mutation per reverted position. When set, those reversions are left in place. Requires the
  /// sparse sequence representation, so this step has no effect under `--dense` and the flag is
  /// then a no-op.
  #[cfg_attr(feature = "clap", clap(long))]
  pub no_flip_parent_child: bool,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[serde(flatten)]
  pub gap_fill_args: GapFillArgs,
}

#[derive(Debug, Clone)]
pub struct TreetimeOptimizeArgs {
  pub alignment: AlignmentArgs,
  pub tree: PathBuf,
  pub alphabet_args: AlphabetArgs,
  pub model_args: ModelArgs,
  pub dense: Option<bool>,
  pub output: OutputCoreArgs,
  pub divergence_units: DivergenceUnits,
  pub output_augur_node_data: Option<PathBuf>,
  pub output_gtr: Option<PathBuf>,
  pub output_selection: Vec<OptimizeOutputSelection>,
  pub topology_order: TopologyOrderArgs,
  pub max_iter: usize,
  pub dp: f64,
  pub damping: f64,
  pub branch_length_initial_guess: InitialGuessMode,
  pub opt_method: BranchOptMethod,
  pub no_indels: bool,
  pub reroot: Option<OptimizeRerootMethod>,
  pub reroot_tips: Vec<String>,
  pub keep_root: bool,
  pub topology_ops: TopologyOps,
  pub gap_fill_args: GapFillArgs,
}

impl TreetimeOptimizeArgs {
  pub fn tree(&self) -> &Path {
    &self.tree
  }

  pub fn reroot_spec(&self) -> Option<RerootSpec> {
    if self.keep_root {
      return None;
    }

    if let Some(method) = self.reroot {
      return Some(RerootSpec::Method(RerootMethod::from(method)));
    }

    if !self.reroot_tips.is_empty() {
      return Some(RerootSpec::Tips(self.reroot_tips.clone()));
    }

    None
  }
}

impl TryFrom<TreetimeOptimizeArgsRaw> for TreetimeOptimizeArgs {
  type Error = Report;

  fn try_from(raw: TreetimeOptimizeArgsRaw) -> Result<Self, Report> {
    let tree = raw
      .tree
      .ok_or_else(|| missing_required_args::<TreetimeOptimizeArgsRaw>(&["tree"]))?;
    Ok(Self {
      alignment: raw.alignment,
      tree,
      alphabet_args: raw.alphabet_args,
      model_args: raw.model_args,
      dense: raw.dense,
      output: raw.output,
      divergence_units: raw.divergence_units,
      output_augur_node_data: raw.output_augur_node_data,
      output_gtr: raw.output_gtr,
      output_selection: raw.output_selection,
      topology_order: raw.topology_order,
      max_iter: raw.max_iter,
      dp: raw.dp,
      damping: raw.damping,
      branch_length_initial_guess: raw.branch_length_initial_guess.into(),
      opt_method: raw.opt_method.into(),
      no_indels: raw.no_indels,
      reroot: raw.reroot,
      reroot_tips: raw.reroot_tips,
      keep_root: raw.keep_root,
      topology_ops: TopologyOps {
        collapse_short_branches: !raw.no_collapse_short_branches,
        merge_siblings: !raw.no_merge_siblings,
        flip_parent_child: !raw.no_flip_parent_child,
      },
      gap_fill_args: raw.gap_fill_args,
    })
  }
}
