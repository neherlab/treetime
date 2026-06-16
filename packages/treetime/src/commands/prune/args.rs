use crate::commands::shared::alignment::AlignmentArgs;
use crate::commands::shared::alphabet::AlphabetArgs;
use crate::commands::shared::output::{OutputCoreArgs, PruneOutputSelection, TopologyOrderArgs};
#[cfg(feature = "clap")]
use clap::ValueHint;
use serde::{Deserialize, Serialize};
use smart_default::SmartDefault;
use std::fmt::Debug;
use std::path::PathBuf;

#[derive(Debug, SmartDefault, Serialize, Deserialize)]
#[serde(default)]
#[cfg_attr(feature = "clap", derive(clap::Parser))]
pub struct TreetimePruneArgs {
  #[cfg_attr(feature = "clap", clap(flatten))]
  pub alignment: AlignmentArgs,

  /// Name of file containing the tree in newick, nexus, or phylip format.
  #[cfg_attr(feature = "clap", clap(long, short = 't'))]
  #[cfg_attr(feature = "clap", clap(value_hint = ValueHint::FilePath))]
  pub tree: PathBuf,

  #[cfg_attr(feature = "clap", clap(flatten))]
  pub alphabet_args: AlphabetArgs,

  #[cfg_attr(feature = "clap", clap(flatten))]
  pub output: OutputCoreArgs,

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
  pub output_selection: Vec<PruneOutputSelection>,

  #[cfg_attr(feature = "clap", clap(flatten))]
  pub topology_order: TopologyOrderArgs,

  /// Threshold value for pruning of branches
  ///
  /// If set, prune branches with a length below this value
  #[cfg_attr(feature = "clap", clap(long, short = 's', value_name = "THRESHOLD"))]
  pub prune_short: Option<f64>,

  /// Prune empty branches
  ///
  /// If set, prune any branch that does not have a mutation or other state transition mapped to it.
  ///
  /// Requires --alignment
  #[cfg_attr(feature = "clap", clap(long, short = 'e'))]
  pub prune_empty: bool,

  /// Merge branches in polytomies that share identical mutations
  ///
  /// When sibling branches in a polytomy (node with >2 children) carry identical substitutions,
  /// they are grouped under a new internal node. The shared mutations move to the new branch
  /// (parent to new node), and only unique mutations remain on children's edges.
  /// Reduces tree builder artifacts from arbitrary binary resolution of polytomies.
  ///
  /// Requires --alignment
  #[cfg_attr(feature = "clap", clap(long, short = 'm'))]
  pub merge_shared_mutations: bool,

  /// List of node names to prune
  ///
  /// List of node names to remove from the tree, comma-separated (,)
  ///
  /// Use --prune-nodes-list-delimiter to specify a different delimiter.
  #[cfg_attr(feature = "clap", clap(long, short = 'n', value_name = "NODE_NAMES"))]
  pub prune_nodes_list: Option<String>,

  /// Name separator for `--prune-nodes-list`
  ///
  /// String used to separate node names in the list given to (--prune-nodes-list). Make sure to correctly quote and escape the delimiter according to your shell.
  #[cfg_attr(feature = "clap", clap(long, default_value = ",", value_name = "DELIMITER"))]
  #[default = ',']
  pub prune_nodes_list_delimiter: char,

  /// File containing list of node names to prune
  ///
  /// Path to a file containing node names to remove from the tree, newline-delimited (\n).
  ///
  /// Use '-' to read from standard input (stdin).
  ///
  /// Use --prune-nodes-list-file-delimiter to specify a different delimiter.
  #[cfg_attr(feature = "clap", clap(long, short = 'N', value_hint = ValueHint::FilePath, value_name = "FILEPATH"))]
  pub prune_nodes_list_file: Option<PathBuf>,

  /// Separator for node names in the list file
  ///
  /// Character or string used to separate node names in the list file.
  #[cfg_attr(feature = "clap", clap(long, default_value = "\n", value_name = "DELIMITER"))]
  #[default = '\n']
  pub prune_nodes_list_file_delimiter: char,
}
