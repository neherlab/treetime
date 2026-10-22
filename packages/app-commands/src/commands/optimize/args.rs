use crate::commands::shared::alignment::AlignmentArgs;
use crate::commands::shared::alphabet::AlphabetArgs;
use crate::commands::shared::branch_lengths::BranchLengthArgs;
use crate::commands::shared::config::ConfigArgs;
use crate::commands::shared::dependent_settings::check_applies;
use crate::commands::shared::gap_fill::GapFillArgs;
use crate::commands::shared::gtr_output::GtrOutputArgs;
use crate::commands::shared::jobs::JobsArgs;
use crate::commands::shared::mirror_enum::impl_from_mirror;
use crate::commands::shared::model::{ModelArgs, PartitionRatesArgs};
use crate::commands::shared::output_args::{OptimizeOutputSelectionArgs, OutputCoreArgs};
use crate::commands::shared::representation::RepresentationArgs;
use crate::commands::shared::required::missing_required_args;
use crate::commands::shared::reroot::{RootingArgs, resolve_outgroup};
use crate::commands::shared::sequence_inputs::MissingSequenceArgs;
use crate::commands::shared::topology_order_args::TopologyOrderArgs;
use crate::commands::shared::tree_input::{TreeArgs, TreeDialectArgs};
use crate::commands::shared::zero_based::ZeroBasedArgs;
#[cfg(feature = "clap")]
use clap::ValueHint;
use derive_more::{Deref, DerefMut};
use deser::{Deserialize, Serialize};
use eyre::Report;
use schemars::JsonSchema;
use smart_default::SmartDefault;
use std::path::{Path, PathBuf};
use treetime::clock::find_best_root::params::{RerootMethod, RerootSpec};
use treetime::gtr::get_gtr::GtrModel;
use treetime::optimize::indel::IndelPolicy;
use treetime::optimize::params::{BranchOptMethod, InitialGuessMode, TopologyOps};
use treetime_graph::graph::Graph;
use treetime_graph::node_names::NodeNames;
use treetime_schema::{schema_defaults, skip_serializing_optionals};
use treetime_utils::o;

#[derive(Debug, Clone, Deref, DerefMut)]
pub struct TreetimeOptimizeArgs {
  #[deref]
  #[deref_mut]
  pub(crate) raw: TreetimeOptimizeArgsRaw,
  pub(crate) tree: PathBuf,
  pub(crate) initial_branch_lengths: InitialGuessMode,
  pub(crate) optimizer: BranchOptMethod,
  pub(crate) indel_policy: IndelPolicy,
  pub(crate) topology_ops: TopologyOps,
  pub(crate) model: GtrModel,
}

impl TreetimeOptimizeArgs {
  pub fn tree(&self) -> &Path {
    &self.tree
  }

  pub fn reroot_spec(&self, graph: &Graph, names: &NodeNames) -> Result<Option<RerootSpec>, Report> {
    if self.rooting.keep_root {
      return Ok(None);
    }

    if let Some(method) = self.reroot {
      return Ok(Some(RerootSpec::Method(RerootMethod::from(method))));
    }

    if !self.rooting.outgroup.is_empty() {
      return Ok(Some(RerootSpec::Tips(resolve_outgroup(
        &self.rooting.outgroup,
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
      .tree
      .clone()
      .ok_or_else(|| missing_required_args::<TreetimeOptimizeArgsRaw>(&["tree"]))?;
    check_dependent_settings(&raw)?;
    Ok(Self {
      tree,
      initial_branch_lengths: raw.initial_branch_lengths.into(),
      optimizer: raw.optimizer.into(),
      indel_policy: raw.branch_lengths.indel_policy(),
      topology_ops: TopologyOps {
        collapse_short_branches: !raw.no_collapse_zero_branches,
        merge_siblings: !raw.no_merge_shared_mutations,
        flip_parent_child: !raw.no_hoist_reversions,
      },
      model: raw.model_args.model()?,
      raw,
    })
  }
}

fn check_dependent_settings(raw: &TreetimeOptimizeArgsRaw) -> Result<(), Report> {
  raw
    .partition_rates_args
    .check_partitions(raw.alignment.alignment.len())?;
  let dense = raw.representation.dense() == Some(true);
  for (setting, given) in [
    ("no_merge_shared_mutations", raw.no_merge_shared_mutations),
    ("no_hoist_reversions", raw.no_hoist_reversions),
  ] {
    check_applies(setting, given, "representation", !dense, || {
      o!("applies only to the sparse representation, but --representation is dense")
    })?;
  }
  Ok(())
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

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub tree: TreeArgs,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub newick_dialect: TreeDialectArgs,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub missing_sequences: MissingSequenceArgs,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub alphabet_args: AlphabetArgs,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub zero_based: ZeroBasedArgs,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub model_args: ModelArgs,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub partition_rates_args: PartitionRatesArgs,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub representation: RepresentationArgs,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub output: OutputCoreArgs,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub branch_lengths: BranchLengthArgs,

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

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub gtr_output: GtrOutputArgs,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub output_selection: OptimizeOutputSelectionArgs,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub topology_order: TopologyOrderArgs,

  /// Maximum number of iterations of the inference loop
  ///
  /// Polytomy resolution and coalescent models of the time inference need at least 2 iterations.
  #[cfg_attr(feature = "clap", clap(long, default_value_t = 10, help_heading = "Branch lengths"))]
  #[default = 10]
  pub max_iter: usize,

  /// Likelihood convergence threshold. The loop stops when successive
  /// likelihoods differ by less than this value, or when a 2-cycle with
  /// amplitude below this value is detected.
  #[cfg_attr(feature = "clap", clap(long, default_value_t = 0.1, help_heading = "Branch lengths"))]
  #[default = 0.1]
  pub tolerance: f64,

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
  #[cfg_attr(feature = "clap", clap(long, value_enum, default_value_t = InitialGuessModeCli::Auto, help_heading = "Branch lengths"))]
  pub initial_branch_lengths: InitialGuessModeCli,

  /// Per-edge optimizer. Both methods optimize in $\sqrt{t}$ space. Values are
  /// `brent` (default) and `newton`.
  #[cfg_attr(feature = "clap", clap(long, value_enum, default_value_t = BranchOptMethodCli::default(), help_heading = "Branch lengths"))]
  pub optimizer: BranchOptMethodCli,

  /// Reroot the tree by minimizing root-to-tip divergence variance.
  ///
  /// By default, optimize keeps the input root. Pass --reroot=min-dev to enable divergence-based rerooting.
  /// Date-dependent methods (least-squares, oldest) are available in the timetree and clock commands.
  #[cfg_attr(
    feature = "clap",
    clap(long, value_enum, conflicts_with = "outgroup", help_heading = "Rooting")
  )]
  pub reroot: Option<OptimizeRerootMethod>,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub rooting: RootingArgs,

  /// Disable collapsing of internal branches whose optimized length is zero.
  ///
  /// By default the optimize loop contracts internal edges the per-edge optimizer drove to
  /// exactly zero that carry no substitutions or indels, turning the resulting binary nodes into
  /// polytomies. When set, such edges are kept in the output tree with length zero.
  #[cfg_attr(feature = "clap", clap(long, help_heading = "Branch lengths"))]
  pub no_collapse_zero_branches: bool,

  /// Disable merging of polytomy siblings that share substitutions.
  ///
  /// By default, sibling branches in a polytomy that carry identical substitutions are grouped
  /// under a new internal node. The step and the flag apply only to the sparse representation.
  #[cfg_attr(feature = "clap", clap(long, help_heading = "Branch lengths"))]
  pub no_merge_shared_mutations: bool,

  /// Disable the reversion hoist in polytomies.
  ///
  /// By default, when a child branch carries the exact reversion of a substitution on the node's
  /// parent branch, a new node is inserted that groups the node with that child, removing one
  /// mutation per reverted position. When set, those reversions are left in place. The step and the flag apply
  /// only to the sparse representation.
  #[cfg_attr(feature = "clap", clap(long, help_heading = "Branch lengths"))]
  pub no_hoist_reversions: bool,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub gap_fill_args: GapFillArgs,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(skip)]
  #[deser(skip)]
  pub jobs_args: JobsArgs,
}

/// Reroot methods available in the optimize command.
///
/// Only date-free methods are valid here because optimize has no sampling dates.
/// Date-dependent methods (least-squares, oldest) are available
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
