use crate::commands::shared::alignment::AlignmentArgs;
use crate::commands::shared::alphabet::AlphabetArgs;
use crate::commands::shared::config::ConfigArgs;
use crate::commands::shared::gap_fill::GapFillArgs;
use crate::commands::shared::mirror_enum::impl_from_mirror;
use crate::commands::shared::model::{ModelArgs, PartitionRatesArgs};
use crate::commands::shared::output_args::{DivergenceUnits, OptimizeOutputSelection, OutputCoreArgs};
use crate::commands::shared::required::missing_required_args;
use crate::commands::shared::reroot::resolve_reroot_tips;
use crate::commands::shared::topology_order_args::TopologyOrderArgs;
use crate::commands::shared::tree_input::TreeDialectArgs;
#[cfg(feature = "clap")]
use clap::ValueHint;
use derive_more::{Deref, DerefMut};
use deser::{Deserialize, Serialize};
use eyre::Report;
use schemars::JsonSchema;
use smart_default::SmartDefault;
use std::path::{Path, PathBuf};
use treetime::clock::find_best_root::params::{RerootMethod, RerootSpec};
use treetime::optimize::indel::IndelPolicy;
use treetime::optimize::params::{BranchOptMethod, InitialGuessMode, TopologyOps};
use treetime_graph::graph::Graph;
use treetime_graph::node_names::NodeNames;
use treetime_schema::{schema_defaults, skip_serializing_optionals};

#[derive(Debug, Clone, Deref, DerefMut)]
pub struct TreetimeOptimizeArgs {
  #[deref]
  #[deref_mut]
  pub(crate) raw: TreetimeOptimizeArgsRaw,
  pub(crate) tree: PathBuf,
  pub(crate) branch_length_initial_guess: InitialGuessMode,
  pub(crate) opt_method: BranchOptMethod,
  pub(crate) indel_policy: IndelPolicy,
  pub(crate) topology_ops: TopologyOps,
}

impl TreetimeOptimizeArgs {
  pub fn tree(&self) -> &Path {
    &self.tree
  }

  pub fn reroot_spec(&self, graph: &Graph, names: &NodeNames) -> Result<Option<RerootSpec>, Report> {
    if self.keep_root {
      return Ok(None);
    }

    if let Some(method) = self.reroot {
      return Ok(Some(RerootSpec::Method(RerootMethod::from(method))));
    }

    if !self.reroot_tips.is_empty() {
      return Ok(Some(RerootSpec::Tips(resolve_reroot_tips(
        &self.reroot_tips,
        graph,
        names,
      )?)));
    }

    Ok(None)
  }
}

impl TryFrom<TreetimeOptimizeArgsRaw> for TreetimeOptimizeArgs {
  type Error = Report;

  fn try_from(raw: TreetimeOptimizeArgsRaw) -> Result<Self, Report> {
    let tree = raw
      .tree
      .clone()
      .ok_or_else(|| missing_required_args::<TreetimeOptimizeArgsRaw>(&["tree"]))?;
    Ok(Self {
      tree,
      branch_length_initial_guess: raw.branch_length_initial_guess.into(),
      opt_method: raw.opt_method.into(),
      indel_policy: if raw.no_indels {
        IndelPolicy::Exclude
      } else {
        IndelPolicy::Include
      },
      topology_ops: TopologyOps {
        collapse_short_branches: !raw.no_collapse_short_branches,
        merge_siblings: !raw.no_merge_siblings,
        flip_parent_child: !raw.no_flip_parent_child,
      },
      raw,
    })
  }
}

#[derive(Debug, Clone, SmartDefault, JsonSchema, Serialize, Deserialize)]
#[deser(skip_serializing_optionals)]
#[schemars(transform = skip_serializing_optionals)]
#[schemars(default, deny_unknown_fields)]
#[schemars(transform = schema_defaults::<Self>)]
#[deser(default, deny_unknown_fields)]
#[cfg_attr(feature = "clap", derive(clap::Parser))]
#[schemars(rename = "OptimizeConfig")]
pub struct TreetimeOptimizeArgsRaw {
  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(skip)]
  #[deser(skip)]
  pub config_args: ConfigArgs,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  #[schemars(extend("x-path" = "input"))]
  pub alignment: AlignmentArgs,

  /// Tree in Newick format.
  #[cfg_attr(feature = "clap", clap(long, short = 't', help_heading = "Input data"))]
  #[cfg_attr(feature = "clap", clap(value_hint = ValueHint::FilePath))]
  #[schemars(extend("x-path" = "input"))]
  pub tree: Option<PathBuf>,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub tree_dialect: TreeDialectArgs,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub alphabet_args: AlphabetArgs,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub model_args: ModelArgs,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub partition_rates_args: PartitionRatesArgs,

  /// Use dense representation of sequences on the tree
  ///
  /// Dense mode stores full probability vectors at every alignment position for each
  /// node. Sparse mode stores only variable positions. Dense is more accurate when
  /// branches are long and many sites change, but uses more memory.
  #[cfg_attr(feature = "clap", clap(long, help_heading = "Ancestral reconstruction"))]
  pub dense: Option<bool>,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub output: OutputCoreArgs,

  /// Units for divergence values in augur node data JSON output.
  ///
  /// `mutations-per-site` (default): branch divergence as substitutions per site.
  /// `mutations`: count of reconstructed substitutions per branch that change the state:
  /// gaps, `N`, and ambiguity codes compatible with the parent state do not count.
  #[cfg_attr(feature = "clap", clap(long, value_enum, default_value_t = DivergenceUnits::default(), help_heading = "Output"))]
  pub divergence_units: DivergenceUnits,

  /// Path to output augur-compatible node data JSON.
  ///
  /// Contains per-node optimized branch lengths (divergence, substitutions per
  /// site) and the input alignment and tree paths. The output is compatible with
  /// augur export v2 --node-data.
  ///
  /// Takes precedence over paths configured with `--output-all` and `--output-selection`.
  #[cfg_attr(feature = "clap", clap(long, value_hint = ValueHint::FilePath, help_heading = "Output"))]
  #[schemars(extend("x-path" = "output"))]
  pub output_augur_node_data: Option<PathBuf>,

  /// Path to output GTR model JSON.
  ///
  /// Takes precedence over paths configured with `--output-all` and `--output-selection`.
  #[cfg_attr(feature = "clap", clap(long, value_hint = ValueHint::FilePath, help_heading = "Output"))]
  #[schemars(extend("x-path" = "output"))]
  pub output_gtr: Option<PathBuf>,

  /// Comma-separated list of outputs to produce with `--output-all`.
  ///
  /// Restricts which outputs `--output-all` writes. Special value `all` expands to every output
  /// available for this command. Requires `--output-all`. Per-file flags are always honored
  /// regardless of this selection.
  ///
  /// A selected output that the run has no data for is skipped without a message, for example
  /// the substitution model of a run that fits none. A per-file flag for such an output fails.
  #[cfg_attr(
    feature = "clap",
    clap(long, value_delimiter = ',', requires = "output_all", help_heading = "Output")
  )]
  pub output_selection: Vec<OptimizeOutputSelection>,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub topology_order: TopologyOrderArgs,

  /// Maximum number of iterations
  #[cfg_attr(feature = "clap", clap(long, default_value_t = 10, help_heading = "Branch lengths"))]
  #[default = 10]
  pub max_iter: usize,

  /// Likelihood convergence threshold. The loop stops when successive
  /// likelihoods differ by less than this value, or when a 2-cycle with
  /// amplitude below this value is detected.
  #[cfg_attr(feature = "clap", clap(long, default_value_t = 0.1, help_heading = "Branch lengths"))]
  #[default = 0.1]
  pub dp: f64,

  /// Damping factor $d$ for outer-loop updates: $b=b_{new}(1-d)+b_{old}d$,
  /// where $d=max(damping^{i+1},0.01)$. Higher values reduce oscillation;
  /// zero disables damping. Must be in $[0,1)$.
  #[cfg_attr(feature = "clap", clap(long, default_value_t = 0.0, help_heading = "Branch lengths"))]
  #[default = 0.0]
  pub damping: f64,

  /// Initial branch length estimate before the per-edge optimization.
  ///
  /// - auto: estimate only edges with missing or invalid branch lengths,
  ///   preserve valid input values (default)
  /// - always: estimate all edges, overwriting input branch lengths
  /// - never: use input branch lengths as-is; fails if any are missing
  #[cfg_attr(feature = "clap", clap(long = "branch-length-initial-guess", value_enum, default_value_t = InitialGuessModeCli::Auto, help_heading = "Branch lengths"))]
  pub branch_length_initial_guess: InitialGuessModeCli,

  /// Per-edge optimizer. Both methods optimize in $\sqrt{t}$ space. Values are
  /// `brent` (default) and `newton`.
  #[cfg_attr(feature = "clap", clap(long = "opt-method", value_enum, default_value_t = BranchOptMethodCli::default(), help_heading = "Branch lengths"))]
  pub opt_method: BranchOptMethodCli,

  /// Disable indel (insertion/deletion) contributions to branch-length
  /// optimization.
  ///
  /// When set, the optimizer uses substitution-only likelihood, matching
  /// standard phylogenetic tools (RAxML, IQ-TREE, PhyML, BEAST) and
  /// enabling v0 parity testing. Default: indels enabled.
  #[cfg_attr(feature = "clap", clap(long, help_heading = "Branch lengths"))]
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
    help_heading = "Rooting",
  ))]
  pub reroot: Option<OptimizeRerootMethod>,

  /// Reroot on the branch leading to a tip or the MRCA of a comma-separated tip list.
  #[cfg_attr(
    feature = "clap",
    clap(long, value_delimiter = ',', conflicts_with = "reroot", help_heading = "Rooting")
  )]
  pub reroot_tips: Vec<String>,

  /// Keep the input tree root instead of rerooting.
  ///
  /// Optimize keeps the input root by default; this flag is the explicit form and
  /// is mutually exclusive with the reroot options.
  #[cfg_attr(feature = "clap", clap(long, conflicts_with_all = ["reroot", "reroot_tips"], help_heading = "Rooting"))]
  pub keep_root: bool,

  /// Disable collapsing of internal branches whose optimized length is zero.
  ///
  /// By default the optimize loop contracts internal edges the per-edge optimizer drove to
  /// exactly zero that carry no substitutions or indels, turning the resulting binary nodes into
  /// polytomies. When set, such edges are kept in the output tree with length zero.
  #[cfg_attr(feature = "clap", clap(long, help_heading = "Branch lengths"))]
  pub no_collapse_short_branches: bool,

  /// Disable merging of polytomy siblings that share substitutions.
  ///
  /// By default, sibling branches in a polytomy that carry identical substitutions are grouped
  /// under a new internal node. Requires the sparse sequence representation, so this step has no
  /// effect under `--dense` and the flag is then a no-op.
  #[cfg_attr(feature = "clap", clap(long, help_heading = "Branch lengths"))]
  pub no_merge_siblings: bool,

  /// Disable the flip-parent-child step (reversion hoist) in polytomies.
  ///
  /// By default, when a child branch carries the exact reversion of a substitution on the node's
  /// parent branch, a new node is inserted that groups the node with that child, removing one
  /// mutation per reverted position. When set, those reversions are left in place. Requires the
  /// sparse sequence representation, so this step has no effect under `--dense` and the flag is
  /// then a no-op.
  #[cfg_attr(feature = "clap", clap(long, help_heading = "Branch lengths"))]
  pub no_flip_parent_child: bool,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub gap_fill_args: GapFillArgs,
}

/// Reroot methods available in the optimize command.
///
/// Only date-free methods are valid here because optimize has no sampling dates.
/// Date-dependent methods (least-squares, oldest, clock-filter) are available
/// in the timetree and clock commands.
#[derive(Copy, Debug, Clone, PartialEq, Eq, JsonSchema, Serialize, Deserialize)]
#[cfg_attr(feature = "clap", derive(clap::ValueEnum))]
#[schemars(rename_all = "kebab-case")]
#[deser(rename_all = "kebab-case")]
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
/// Controls how the maximum-likelihood branch length of each edge is found:
/// Brent's method or Newton-Raphson, both in $\sqrt{t}$ space.
#[derive(Copy, Clone, Debug, PartialEq, Eq, SmartDefault, JsonSchema, Serialize, Deserialize)]
#[cfg_attr(feature = "clap", derive(clap::ValueEnum))]
#[schemars(rename_all = "kebab-case")]
#[deser(rename_all = "kebab-case")]
#[schemars(rename = "BranchOptMethod")]
pub enum BranchOptMethodCli {
  /// Brent's method in $\sqrt{t}$ space.
  ///
  /// Derivative-free line search that combines parabolic interpolation with
  /// golden-section steps. Optimizing $s=\sqrt{t}$ instead of $t$ smooths the
  /// likelihood near short branches, so the interpolation fits better. v0
  /// uses the same method and parameterization, with its own search bounds
  /// and stopping rule, so branch lengths can differ from v0. This is the
  /// default.
  #[default]
  Brent,

  /// Newton-Raphson in $\sqrt{t}$ space.
  ///
  /// Uses the first and second derivatives of the likelihood, including the
  /// indel term, and converges in few iterations on well-behaved edges.
  /// Choose it to compare against the default; branch lengths can differ
  /// from the default and from v0.
  Newton,
}

impl_from_mirror!(BranchOptMethodCli => BranchOptMethod { Brent, Newton });

/// Controls whether marginal reconstruction estimates initial branch lengths
/// from substitutions divided by effective alignment length. Preserving valid
/// input lengths can provide a better Newton starting point.
#[derive(Copy, Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Default, JsonSchema, Serialize, Deserialize)]
#[cfg_attr(feature = "clap", derive(clap::ValueEnum))]
#[schemars(rename_all = "kebab-case")]
#[deser(rename_all = "kebab-case")]
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

impl_from_mirror!(InitialGuessModeCli => InitialGuessMode { Auto, Always, Never });
