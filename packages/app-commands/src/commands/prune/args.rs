use crate::commands::shared::alignment::AlignmentArgs;
use crate::commands::shared::alphabet::AlphabetArgs;
use crate::commands::shared::config::ConfigArgs;
use crate::commands::shared::gap_fill::GapFillArgs;
use crate::commands::shared::jobs::JobsArgs;
use crate::commands::shared::output_args::{OutputCoreArgs, PruneOutputSelectionArgs};
use crate::commands::shared::required::missing_required_args;
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
use treetime_schema::{schema_defaults, skip_serializing_optionals};

#[derive(Debug, Clone, Deref, DerefMut)]
pub struct TreetimePruneArgs {
  #[deref]
  #[deref_mut]
  pub(crate) raw: TreetimePruneArgsRaw,
  pub(crate) tree: PathBuf,
}

impl TreetimePruneArgs {
  pub fn tree(&self) -> &Path {
    &self.tree
  }
}

impl TryFrom<TreetimePruneArgsRaw> for TreetimePruneArgs {
  type Error = Report;

  fn try_from(raw: TreetimePruneArgsRaw) -> Result<Self, Report> {
    let tree = raw
      .tree
      .tree
      .clone()
      .ok_or_else(|| missing_required_args::<TreetimePruneArgsRaw>(&["tree"]))?;
    Ok(Self { raw, tree })
  }
}

#[derive(Debug, Clone, SmartDefault, JsonSchema, Serialize, Deserialize)]
#[deser(skip_serializing_optionals)]
#[schemars(transform = skip_serializing_optionals)]
#[schemars(default, deny_unknown_fields)]
#[schemars(transform = schema_defaults::<Self>)]
#[deser(default, deny_unknown_fields)]
#[cfg_attr(feature = "clap", derive(clap::Parser))]
#[schemars(rename = "PruneConfig")]
pub struct TreetimePruneArgsRaw {
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
  pub output: OutputCoreArgs,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub output_selection: PruneOutputSelectionArgs,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub topology_order: TopologyOrderArgs,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub gap_fill_args: GapFillArgs,

  /// Collapse branches shorter than this length
  ///
  /// Removes each internal branch with a length below LEN and attaches the children of its node to the parent node.
  #[cfg_attr(feature = "clap", clap(long, value_name = "LEN", help_heading = "Pruning"))]
  pub collapse_short: Option<f64>,

  /// Collapse branches without mutations
  ///
  /// Removes each internal branch that has no mutation or other state transition mapped to it.
  ///
  /// Requires --alignment
  #[cfg_attr(feature = "clap", clap(long, help_heading = "Pruning"))]
  pub collapse_empty: bool,

  /// Merge branches in polytomies that share identical mutations
  ///
  /// When sibling branches in a polytomy (node with >2 children) carry identical substitutions,
  /// they are grouped under a new internal node. The shared mutations move to the new branch
  /// (parent to new node), and only unique mutations remain on children's edges.
  /// Reduces tree builder artifacts from arbitrary binary resolution of polytomies.
  ///
  /// Requires --alignment
  #[cfg_attr(feature = "clap", clap(long, help_heading = "Pruning"))]
  pub merge_shared_mutations: bool,

  /// Names of nodes to remove from the tree
  ///
  /// Accepts several names, separated by commas or spaces.
  #[cfg_attr(
    feature = "clap",
    clap(long, num_args = 1.., value_delimiter = ',', value_name = "NAME", help_heading = "Pruning")
  )]
  pub remove: Vec<String>,

  /// File with names of nodes to remove from the tree, one name per line
  ///
  /// Surrounding whitespace is trimmed and empty lines are ignored. Use '-' to read from standard input (stdin).
  #[cfg_attr(feature = "clap", clap(long, value_hint = ValueHint::FilePath, value_name = "FILE", help_heading = "Pruning"))]
  #[schemars(extend("x-path" = "input"))]
  pub remove_file: Option<PathBuf>,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(skip)]
  #[deser(skip)]
  pub jobs_args: JobsArgs,
}
