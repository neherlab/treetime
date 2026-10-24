use crate::commands::shared::alignment::{AlignmentArgs, check_required_alignment_paths};
use crate::commands::shared::alphabet::AlphabetArgs;
use crate::commands::shared::branch_lengths::{BranchLengthArgs, indel_policy};
use crate::commands::shared::config::ConfigArgs;
use crate::commands::shared::dependent_settings::check_applies;
use crate::commands::shared::gap_fill::GapFillArgs;
use crate::commands::shared::gtr_output::GtrOutputArgs;
use crate::commands::shared::jobs::JobsArgs;
use crate::commands::shared::mirror_enum::impl_from_mirror;
use crate::commands::shared::model::{
  ModelArgs, PartitionRatesArgs, check_partition_rates, resolve_model, resolve_partition_rates,
};
use crate::commands::shared::output_args::{
  DivergenceUnits, OptimizeOutputSelectionArgs, OutputCoreArgs, OutputRequest, TreeOutputPaths, selected_outputs,
};
use crate::commands::shared::partition::{PartitionArgs, partition_entries};
use crate::commands::shared::representation::{RepresentationArgs, representation_dense};
use crate::commands::shared::required::missing_required_args;
use crate::commands::shared::reroot::{RerootTarget, RootingArgs};
use crate::commands::shared::sequence_inputs::{MissingSequenceArgs, SequenceInputs, missing_leaf_limit};
use crate::commands::shared::topology_order_args::{TopologyOrder, TopologyOrderArgs};
use crate::commands::shared::tree_input::{TreeArgs, TreeDialectArgs};
#[cfg(feature = "clap")]
use crate::commands::shared::value_names::value_name_from_hint;
use crate::commands::shared::zero_based::ZeroBasedArgs;
use app_output::command_outputs::CommandKind;
use app_output::output_conditions::OutputSettings;
use app_output::output_plan::OutputSelection;
#[cfg(feature = "clap")]
use clap::ValueHint;
use deser::{Deserialize, Serialize};
use eyre::Report;
use schemars::JsonSchema;
use smart_default::SmartDefault;
use std::path::PathBuf;
use treetime::gtr::get_gtr::GtrModel;
use treetime::optimize::params::{BranchOptMethod, InitialGuessMode, OptimizeLoopSettings, TopologyOps};
use treetime::reroot::params::RerootMethod;
use treetime::seq::mutation::PositionBase;
use treetime_io::nwk::NewickDialect;
use treetime_schema::{schema_defaults, skip_serializing_optionals};
use treetime_utils::o;

pub const OPTIMIZE_ABOUT: &str = "Optimize the branch lengths of a tree for the likelihood of an alignment";

#[derive(Debug, Clone)]
pub struct TreetimeOptimizeArgs {
  pub(crate) tree: PathBuf,
  pub(crate) newick_dialect: NewickDialect,
  pub(crate) sequences: SequenceInputs,
  pub(crate) positions: PositionBase,
  pub(crate) model: GtrModel,
  pub(crate) dense: Option<bool>,
  pub(crate) initial_branch_lengths: InitialGuessMode,
  pub(crate) reroot: Option<RerootTarget>,
  pub(crate) loop_settings: OptimizeLoopSettings,
  pub(crate) mutation_units: bool,
  pub(crate) topology_order: TopologyOrder,
  pub(crate) outputs: OutputRequest,
}

impl TreetimeOptimizeArgs {
  pub fn alignments(&self) -> &[PathBuf] {
    &self.sequences.alignments
  }
}

impl TryFrom<TreetimeOptimizeArgsRaw> for TreetimeOptimizeArgs {
  type Error = Report;

  #[deny(unused_variables)]
  #[expect(
    clippy::unneeded_field_pattern,
    reason = "the CLI reads its own settings before the conversion, and `..` would hide every setting that no run reads"
  )]
  fn try_from(raw: TreetimeOptimizeArgsRaw) -> Result<Self, Report> {
    let TreetimeOptimizeArgsRaw {
      config_args: _,
      alignment: AlignmentArgs { alignment },
      partition: PartitionArgs { partition },
      tree: TreeArgs { tree },
      newick_dialect: TreeDialectArgs { newick_dialect },
      missing_sequences: MissingSequenceArgs {
        allow_missing_sequences,
      },
      alphabet_args: AlphabetArgs { alphabet },
      zero_based: ZeroBasedArgs { zero_based },
      model_args: ModelArgs { model, model_params },
      partition_rates_args: PartitionRatesArgs { partition_rates },
      representation: RepresentationArgs { representation },
      output,
      branch_lengths: BranchLengthArgs {
        no_indels,
        divergence_units,
      },
      output_node_data,
      gtr_output,
      output_selection,
      topology_order,
      max_iter,
      tolerance,
      damping,
      initial_branch_lengths,
      optimizer,
      reroot,
      rooting: RootingArgs { outgroup, keep_root },
      no_collapse_zero_branches,
      no_merge_shared_mutations,
      no_hoist_reversions,
      gap_fill_args: GapFillArgs { gap_fill },
      jobs_args: _,
    } = raw;
    let tree = tree.ok_or_else(|| missing_required_args::<TreetimeOptimizeArgsRaw>(&["tree"]))?;
    check_required_alignment_paths(&alignment)?;
    check_partition_rates(partition_rates, alignment.len(), partition.len())?;
    let dense = representation_dense(representation);
    check_sparse_only_settings(dense, no_merge_shared_mutations, no_hoist_reversions)?;
    Ok(Self {
      tree,
      newick_dialect: newick_dialect.dialect(),
      outputs: optimize_outputs(output, output_selection, alignment.len(), output_node_data, gtr_output),
      sequences: SequenceInputs {
        alignments: alignment,
        alphabet: alphabet.map(Into::into),
        gap_fill: gap_fill.into(),
        missing_leaves: missing_leaf_limit(allow_missing_sequences),
        partition_entries: partition_entries(&partition)?,
      },
      positions: PositionBase::from_zero_based(zero_based),
      model: resolve_model(model, &model_params)?,
      dense,
      initial_branch_lengths: initial_branch_lengths.into(),
      reroot: optimize_reroot(keep_root, reroot, outgroup),
      loop_settings: OptimizeLoopSettings {
        max_iter,
        dp: tolerance,
        damping,
        opt_method: optimizer.into(),
        indel_policy: indel_policy(no_indels),
        topology_ops: TopologyOps {
          collapse_short_branches: !no_collapse_zero_branches,
          merge_siblings: !no_merge_shared_mutations,
          flip_parent_child: !no_hoist_reversions,
        },
        partition_rates: resolve_partition_rates(partition_rates),
      },
      mutation_units: divergence_units == DivergenceUnits::Mutations,
      topology_order: optimize_topology_order(topology_order)?,
    })
  }
}

#[deny(unused_variables)]
fn optimize_outputs(
  OutputCoreArgs {
    output_all,
    newick_annotations,
    output_tree_nwk,
    output_tree_nexus,
    output_tree_auspice,
    output_tree_mat_pb,
    output_tree_mat_json,
    output_tree_graph_json,
    output_tree_dot,
  }: OutputCoreArgs,
  OptimizeOutputSelectionArgs { output_selection }: OptimizeOutputSelectionArgs,
  alignment_count: usize,
  output_node_data: Option<PathBuf>,
  GtrOutputArgs { output_model }: GtrOutputArgs,
) -> OutputRequest {
  OutputRequest {
    command: CommandKind::Optimize,
    alignment_count,
    settings: OutputSettings::default(),
    output_all,
    newick_annotations,
    selection: selected_outputs(&output_selection),
    tree_paths: TreeOutputPaths {
      nwk: output_tree_nwk,
      nexus: output_tree_nexus,
      auspice: output_tree_auspice,
      mat_pb: output_tree_mat_pb,
      mat_json: output_tree_mat_json,
      graph_json: output_tree_graph_json,
      dot: output_tree_dot,
    },
    file_paths: vec![
      (OutputSelection::NodeData, output_node_data),
      (OutputSelection::Model, output_model),
    ],
  }
}

#[deny(unused_variables)]
fn optimize_topology_order(
  TopologyOrderArgs {
    topology_order,
    topology_order_tree,
    topology_order_list,
    topology_order_aggregate,
  }: TopologyOrderArgs,
) -> Result<TopologyOrder, Report> {
  TopologyOrder::new(
    topology_order,
    topology_order_tree,
    topology_order_list,
    topology_order_aggregate,
  )
}

fn check_sparse_only_settings(
  dense: Option<bool>,
  no_merge_shared_mutations: bool,
  no_hoist_reversions: bool,
) -> Result<(), Report> {
  let sparse = dense != Some(true);
  for (setting, given) in [
    ("no_merge_shared_mutations", no_merge_shared_mutations),
    ("no_hoist_reversions", no_hoist_reversions),
  ] {
    check_applies(setting, given, "representation", sparse, || {
      o!("applies only to the sparse representation, but --representation is dense")
    })?;
  }
  Ok(())
}

fn optimize_reroot(
  keep_root: bool,
  reroot: Option<OptimizeRerootMethod>,
  outgroup: Vec<String>,
) -> Option<RerootTarget> {
  match (keep_root, reroot) {
    (true, _) => None,
    (false, Some(method)) => Some(RerootTarget::Method(method.into())),
    (false, None) if outgroup.is_empty() => None,
    (false, None) => Some(RerootTarget::Outgroup(outgroup)),
  }
}

#[derive(Debug, Clone, SmartDefault, JsonSchema, Serialize, Deserialize)]
#[deser(skip_serializing_optionals)]
#[schemars(transform = skip_serializing_optionals)]
#[schemars(default, deny_unknown_fields)]
#[schemars(transform = schema_defaults::<Self>)]
#[deser(default, deny_unknown_fields)]
#[cfg_attr(feature = "clap", derive(clap::Parser))]
#[cfg_attr(feature = "clap", command(mut_args(value_name_from_hint)))]
#[schemars(rename = "OptimizeConfig", description = OPTIMIZE_ABOUT)]
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
  pub partition: PartitionArgs,

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

  /// Method that chooses the root of the tree from the branch lengths
  ///
  /// The tree is rerooted once, on the branch lengths after the initial guess, and every tip counts. Without this
  /// setting, the input root stays. The rerooting methods that need dates (`least-squares`, `oldest`) are available
  /// where the inputs have dates.
  #[cfg_attr(
    feature = "clap",
    clap(long, value_enum, conflicts_with = "outgroup", help_heading = "Rooting")
  )]
  pub reroot: Option<OptimizeRerootMethod>,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub rooting: RootingArgs,

  /// Maximum number of iterations of the inference loop
  ///
  /// Polytomy resolution and coalescent models of the time inference need at least 2 iterations.
  #[cfg_attr(
    feature = "clap",
    clap(long, value_name = "N", default_value_t = 10, help_heading = "Branch lengths")
  )]
  #[default = 10]
  pub max_iter: usize,

  /// Convergence threshold of the log-likelihood
  ///
  /// The loop stops when successive log-likelihoods differ by less than this value, or when it detects a 2-cycle with
  /// an amplitude below this value.
  #[cfg_attr(
    feature = "clap",
    clap(long, hide_short_help = true, default_value_t = 0.1, help_heading = "Branch lengths")
  )]
  #[default = 0.1]
  pub tolerance: f64,

  /// Damping of the branch length updates of the loop, from 0 (inclusive) to 1 (exclusive)
  ///
  /// Iteration i updates each branch length b to `b = b_new * (1 - d) + b_old * d`, with
  /// `d = max(damping^(i+1), 0.01)`. Higher values reduce oscillation; 0 switches damping off.
  #[cfg_attr(
    feature = "clap",
    clap(long, hide_short_help = true, default_value_t = 0.0, help_heading = "Branch lengths")
  )]
  #[default = 0.0]
  pub damping: f64,

  /// Initial branch lengths before the optimization of each branch
  #[cfg_attr(feature = "clap", clap(long, hide_short_help = true, value_enum, default_value_t = InitialGuessModeCli::Auto, help_heading = "Branch lengths"))]
  pub initial_branch_lengths: InitialGuessModeCli,

  /// Optimizer of each branch length; both methods optimize sqrt(t) in place of the branch length t
  #[cfg_attr(feature = "clap", clap(long, hide_short_help = true, value_enum, default_value_t = BranchOptMethodCli::default(), help_heading = "Branch lengths"))]
  pub optimizer: BranchOptMethodCli,

  /// Disable collapsing of internal branches whose optimized length is zero.
  ///
  /// By default the optimize loop contracts internal edges the per-edge optimizer drove to
  /// exactly zero that carry no substitutions or indels, turning the resulting binary nodes into
  /// polytomies. When set, such edges are kept in the output tree with length zero.
  #[cfg_attr(
    feature = "clap",
    clap(long, hide_short_help = true, help_heading = "Branch lengths")
  )]
  pub no_collapse_zero_branches: bool,

  /// Disable merging of polytomy siblings that share substitutions.
  ///
  /// By default, sibling branches in a polytomy that carry identical substitutions are grouped
  /// under a new internal node. The step and the flag apply only to the sparse representation.
  #[cfg_attr(
    feature = "clap",
    clap(long, hide_short_help = true, help_heading = "Branch lengths")
  )]
  pub no_merge_shared_mutations: bool,

  /// Disable the reversion hoist in polytomies.
  ///
  /// By default, when a child branch carries the exact reversion of a substitution on the node's
  /// parent branch, a new node is inserted that groups the node with that child, removing one
  /// mutation per reverted position. When set, those reversions are left in place. The step and the flag apply
  /// only to the sparse representation.
  #[cfg_attr(
    feature = "clap",
    clap(long, hide_short_help = true, help_heading = "Branch lengths")
  )]
  pub no_hoist_reversions: bool,

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
  pub gap_fill_args: GapFillArgs,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub output: OutputCoreArgs,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub output_selection: OptimizeOutputSelectionArgs,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub branch_lengths: BranchLengthArgs,

  /// Path to output augur-compatible node data JSON.
  ///
  /// Contains per-node optimized branch lengths (divergence, substitutions per site) and the input alignment and tree
  /// paths. The output is compatible with augur export v2 --node-data.
  #[cfg_attr(feature = "clap", clap(long, hide_short_help = true, value_hint = ValueHint::FilePath, help_heading = "Output"))]
  #[schemars(extend("x-path" = "output"))]
  pub output_node_data: Option<PathBuf>,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub gtr_output: GtrOutputArgs,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub topology_order: TopologyOrderArgs,

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
  /// Minimizes the variance of root-to-tip distances over all tips; more robust than midpoint
  MinDev,
  /// Roots at the middle of the longest tip-to-tip path; uses only the two farthest tips
  Midpoint,
}

impl From<OptimizeRerootMethod> for RerootMethod {
  fn from(m: OptimizeRerootMethod) -> Self {
    match m {
      OptimizeRerootMethod::MinDev => RerootMethod::MinDev,
      OptimizeRerootMethod::Midpoint => RerootMethod::Midpoint,
    }
  }
}

/// Per-edge branch length optimization method.
///
/// Controls how the maximum-likelihood branch length of each edge is found:
/// Optimizer of each branch length: Brent's method or Newton-Raphson, both on sqrt(t).
#[derive(Copy, Clone, Debug, PartialEq, Eq, SmartDefault, JsonSchema, Serialize, Deserialize)]
#[cfg_attr(feature = "clap", derive(clap::ValueEnum))]
#[schemars(rename_all = "kebab-case")]
#[deser(rename_all = "kebab-case")]
#[schemars(rename = "BranchOptMethod")]
pub enum BranchOptMethodCli {
  /// Brent's method on sqrt(t)
  ///
  /// Derivative-free line search that combines parabolic interpolation with golden-section steps. Optimizing
  /// s = sqrt(t) in place of t smooths the likelihood near short branches, so the interpolation fits better. TreeTime
  /// 0.x uses the same method and parameterization, with its own search bounds and stopping rule, so branch lengths can
  /// differ from TreeTime 0.x. This is the default.
  #[default]
  Brent,

  /// Newton-Raphson on sqrt(t)
  ///
  /// Uses the first and second derivatives of the likelihood, including the indel term, and converges in few
  /// iterations on well-behaved branches. Choose it to compare against the default; branch lengths can differ from the
  /// default and from TreeTime 0.x.
  Newton,
}

impl_from_mirror!(BranchOptMethodCli => BranchOptMethod { Brent, Newton });

/// Initial branch lengths: input lengths, or estimates from the number of substitutions divided by the effective
/// alignment length. Valid input lengths can give Newton-Raphson a better starting point.
#[derive(Copy, Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Default, JsonSchema, Serialize, Deserialize)]
#[cfg_attr(feature = "clap", derive(clap::ValueEnum))]
#[schemars(rename_all = "kebab-case")]
#[deser(rename_all = "kebab-case")]
#[schemars(rename = "InitialGuessMode")]
pub enum InitialGuessModeCli {
  /// Estimate only the branches with a missing or invalid length, and keep the valid input lengths
  #[default]
  Auto,
  /// Estimate every branch length, replacing the input lengths
  Always,
  /// Keep the input branch lengths; the run fails when a branch has a missing or invalid length
  Never,
}

impl_from_mirror!(InitialGuessModeCli => InitialGuessMode { Auto, Always, Never });
