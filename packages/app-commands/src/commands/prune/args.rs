use crate::commands::shared::alignment::AlignmentArgs;
use crate::commands::shared::alphabet::AlphabetArgs;
use crate::commands::shared::config::ConfigArgs;
use crate::commands::shared::gtr_output::GtrOutputArgs;
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
  pub gtr_output: GtrOutputArgs,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub output_selection: PruneOutputSelectionArgs,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub topology_order: TopologyOrderArgs,

  /// Threshold value for pruning of branches
  ///
  /// If set, prune branches with a length below this value
  #[cfg_attr(
    feature = "clap",
    clap(long, short = 's', value_name = "THRESHOLD", help_heading = "Pruning")
  )]
  pub prune_short: Option<f64>,

  /// Prune empty branches
  ///
  /// If set, prune any branch that does not have a mutation or other state transition mapped to it.
  ///
  /// Requires --alignment
  #[cfg_attr(feature = "clap", clap(long, short = 'e', help_heading = "Pruning"))]
  pub prune_empty: bool,

  /// Merge branches in polytomies that share identical mutations
  ///
  /// When sibling branches in a polytomy (node with >2 children) carry identical substitutions,
  /// they are grouped under a new internal node. The shared mutations move to the new branch
  /// (parent to new node), and only unique mutations remain on children's edges.
  /// Reduces tree builder artifacts from arbitrary binary resolution of polytomies.
  ///
  /// Requires --alignment
  #[cfg_attr(feature = "clap", clap(long, short = 'm', help_heading = "Pruning"))]
  pub merge_shared_mutations: bool,

  /// List of node names to prune
  ///
  /// List of node names to remove from the tree, comma-separated (,)
  ///
  /// Surrounding whitespace is trimmed and empty names are ignored. Use --prune-nodes-list-delimiter to specify a different delimiter.
  #[cfg_attr(
    feature = "clap",
    clap(long, short = 'n', value_name = "NODE_NAMES", help_heading = "Pruning")
  )]
  pub prune_nodes_list: Option<String>,

  /// Name separator for `--prune-nodes-list`
  ///
  /// String used to separate node names in the list given to (--prune-nodes-list). Make sure to correctly quote and escape the delimiter according to your shell.
  #[cfg_attr(
    feature = "clap",
    clap(long, default_value = ",", value_name = "DELIMITER", help_heading = "Pruning")
  )]
  #[default = ',']
  pub prune_nodes_list_delimiter: char,

  /// File containing list of node names to prune
  ///
  /// Path to a file containing node names to remove from the tree, newline-delimited (\n). Surrounding whitespace is trimmed and empty names are ignored.
  ///
  /// Use '-' to read from standard input (stdin).
  ///
  /// Use --prune-nodes-list-file-delimiter to specify a different delimiter.
  #[cfg_attr(feature = "clap", clap(long, short = 'N', value_hint = ValueHint::FilePath, value_name = "FILEPATH", help_heading = "Pruning"))]
  #[schemars(extend("x-path" = "input"))]
  pub prune_nodes_list_file: Option<PathBuf>,

  /// Separator for node names in the list file
  ///
  /// Character or string used to separate node names in the list file.
  #[cfg_attr(
    feature = "clap",
    clap(long, default_value = "\n", value_name = "DELIMITER", help_heading = "Pruning")
  )]
  #[default = '\n']
  pub prune_nodes_list_file_delimiter: char,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(skip)]
  #[deser(skip)]
  pub jobs_args: JobsArgs,
}
