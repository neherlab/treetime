use crate::commands::shared::alignment::AlignmentArgs;
use crate::commands::shared::alphabet::AlphabetArgs;
use crate::commands::shared::argument_error::argument_error;
use crate::commands::shared::config::ConfigArgs;
use crate::commands::shared::gap_fill::GapFillArgs;
use crate::commands::shared::jobs::JobsArgs;
use crate::commands::shared::output_args::{
  OutputCoreArgs, OutputRequest, PruneOutputSelectionArgs, TreeOutputPaths, selected_outputs,
};
use crate::commands::shared::required::missing_required_args;
use crate::commands::shared::sequence_inputs::{MissingSequenceArgs, SequenceInputs, missing_leaf_limit};
use crate::commands::shared::topology_order_args::{TopologyOrder, TopologyOrderArgs};
use crate::commands::shared::tree_input::{TreeArgs, TreeDialectArgs};
#[cfg(feature = "clap")]
use crate::commands::shared::value_names::value_name_from_hint;
use crate::commands::shared::zero_based::ZeroBasedArgs;
use app_output::command_outputs::CommandKind;
use app_output::output_conditions::OutputSettings;
#[cfg(feature = "clap")]
use clap::ValueHint;
use deser::{Deserialize, Serialize};
use eyre::Report;
use schemars::JsonSchema;
use smart_default::SmartDefault;
use std::path::PathBuf;
use treetime::seq::mutation::PositionBase;
use treetime_io::nwk::NewickDialect;
use treetime_schema::{schema_defaults, skip_serializing_optionals};

pub const PRUNE_ABOUT: &str =
  "Collapse short branches and branches without mutations, and remove named nodes from a tree";

#[derive(Debug, Clone)]
pub struct TreetimePruneArgs {
  pub(crate) tree: PathBuf,
  pub(crate) newick_dialect: NewickDialect,
  pub(crate) sequences: SequenceInputs,
  pub(crate) positions: PositionBase,
  pub(crate) collapse_short: Option<f64>,
  pub(crate) collapse_empty: bool,
  pub(crate) merge_shared_mutations: bool,
  pub(crate) remove: Vec<String>,
  pub(crate) remove_file: Option<PathBuf>,
  pub(crate) topology_order: TopologyOrder,
  pub(crate) outputs: OutputRequest,
}

impl TreetimePruneArgs {
  pub fn alignments(&self) -> &[PathBuf] {
    &self.sequences.alignments
  }
}

impl TryFrom<TreetimePruneArgsRaw> for TreetimePruneArgs {
  type Error = Report;

  #[deny(unused_variables)]
  #[expect(
    clippy::unneeded_field_pattern,
    reason = "the CLI reads its own settings before the conversion, and `..` would hide every setting that no run reads"
  )]
  fn try_from(raw: TreetimePruneArgsRaw) -> Result<Self, Report> {
    let TreetimePruneArgsRaw {
      config_args: _,
      alignment: AlignmentArgs { alignment },
      tree: TreeArgs { tree },
      newick_dialect: TreeDialectArgs { newick_dialect },
      missing_sequences: MissingSequenceArgs {
        allow_missing_sequences,
      },
      alphabet_args: AlphabetArgs { alphabet },
      zero_based: ZeroBasedArgs { zero_based },
      output:
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
        },
      output_selection: PruneOutputSelectionArgs { output_selection },
      topology_order:
        TopologyOrderArgs {
          topology_order,
          topology_order_tree,
          topology_order_list,
          topology_order_aggregate,
        },
      gap_fill_args: GapFillArgs { gap_fill },
      collapse_short,
      collapse_empty,
      merge_shared_mutations,
      remove,
      remove_file,
      jobs_args: _,
    } = raw;
    let tree = tree.ok_or_else(|| missing_required_args::<TreetimePruneArgsRaw>(&["tree"]))?;
    check_sequence_settings(&alignment, collapse_empty, merge_shared_mutations)?;
    Ok(Self {
      tree,
      newick_dialect: newick_dialect.dialect(),
      outputs: OutputRequest {
        command: CommandKind::Prune,
        alignment_count: alignment.len(),
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
        file_paths: vec![],
      },
      sequences: SequenceInputs {
        alignments: alignment,
        alphabet: alphabet.map(Into::into),
        gap_fill: gap_fill.into(),
        missing_leaves: missing_leaf_limit(allow_missing_sequences),
        partition_entries: vec![],
      },
      positions: PositionBase::from_zero_based(zero_based),
      collapse_short,
      collapse_empty,
      merge_shared_mutations,
      remove,
      remove_file,
      topology_order: TopologyOrder::new(
        topology_order,
        topology_order_tree,
        topology_order_list,
        topology_order_aggregate,
      )?,
    })
  }
}

fn check_sequence_settings(
  alignments: &[PathBuf],
  collapse_empty: bool,
  merge_shared_mutations: bool,
) -> Result<(), Report> {
  if collapse_empty && alignments.is_empty() {
    return Err(argument_error(
      &["collapse_empty", "alignment"],
      "The --collapse-empty requires --alignment. Without sequence data, it's not possible to determine which branches lack mutations.",
    ));
  }

  if merge_shared_mutations && alignments.is_empty() {
    return Err(argument_error(
      &["merge_shared_mutations", "alignment"],
      "The --merge-shared-mutations requires --alignment. Without sequence data, it's not possible to determine which branches share mutations.",
    ));
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
#[cfg_attr(feature = "clap", command(mut_args(value_name_from_hint)))]
#[schemars(rename = "PruneConfig", description = PRUNE_ABOUT)]
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
  #[schemars(skip)]
  #[deser(skip)]
  pub jobs_args: JobsArgs,
}
