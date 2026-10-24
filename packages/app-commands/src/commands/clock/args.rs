use crate::commands::shared::clock_filter::ClockFilterArgs;
use crate::commands::shared::config::ConfigArgs;
use crate::commands::shared::dates_input::DatesInput;
use crate::commands::shared::jobs::JobsArgs;
use crate::commands::shared::metadata::{DateColumnArgs, MetadataArgs, MetadataIdArgs, MetadataInput};
use crate::commands::shared::output_args::{
  ClockOutputArgs, ClockOutputSelectionArgs, OutputCoreArgs, OutputRequest, TreeOutputPaths, selected_outputs,
};
use crate::commands::shared::required::missing_required_args;
use crate::commands::shared::reroot::{RerootArgs, RerootTarget, RootingArgs};
use crate::commands::shared::topology_order_args::{TopologyOrder, TopologyOrderArgs};
use crate::commands::shared::tree_input::{TreeArgs, TreeDialectArgs};
#[cfg(feature = "clap")]
use crate::commands::shared::value_names::value_name_from_hint;
use app_output::command_outputs::CommandKind;
use app_output::output_conditions::OutputSettings;
use app_output::output_plan::OutputSelection;
#[cfg(feature = "clap")]
use clap::ValueHint;
use deser::{Deserialize, Serialize};
use eyre::Report;
use itertools::Itertools;
use schemars::JsonSchema;
use smart_default::SmartDefault;
use std::path::PathBuf;
use treetime_io::nwk::NewickDialect;
use treetime_schema::{schema_defaults, skip_serializing_optionals};

pub const CLOCK_ABOUT: &str =
  "Fit a molecular clock by root-to-tip regression, and reroot the tree for the best clock signal";

#[derive(Debug, Clone)]
pub struct TreetimeClockArgs {
  pub(crate) tree: PathBuf,
  pub(crate) newick_dialect: NewickDialect,
  pub(crate) dates: DatesInput,
  pub(crate) clock_filter: f64,
  pub(crate) allow_negative_rate: bool,
  pub(crate) reroot: Option<RerootTarget>,
  pub(crate) topology_order: TopologyOrder,
  pub(crate) outputs: OutputRequest,
}

impl TryFrom<TreetimeClockArgsRaw> for TreetimeClockArgs {
  type Error = Report;

  #[deny(unused_variables)]
  #[expect(
    clippy::unneeded_field_pattern,
    reason = "the CLI reads its own settings before the conversion, and `..` would hide every setting that no run reads"
  )]
  fn try_from(raw: TreetimeClockArgsRaw) -> Result<Self, Report> {
    let TreetimeClockArgsRaw {
      config_args: _,
      tree: TreeArgs { tree },
      newick_dialect: TreeDialectArgs { newick_dialect },
      metadata: MetadataArgs { metadata },
      metadata_id: MetadataIdArgs {
        metadata_id_columns,
        metadata_delimiters,
      },
      date_column: DateColumnArgs {
        date_column,
        date_format,
      },
      clock_filter: ClockFilterArgs {
        clock_filter,
        allow_negative_rate,
      },
      reroot: RerootArgs { reroot },
      rooting: RootingArgs { outgroup, keep_root },
      output,
      clock_output,
      output_rtt_chart,
      output_selection,
      topology_order,
      jobs_args: _,
    } = raw;
    let (tree, metadata) = match (tree, metadata) {
      (Some(tree), Some(metadata)) => (tree, metadata),
      (tree, metadata) => {
        let missing = [("tree", tree.is_none()), ("metadata", metadata.is_none())]
          .into_iter()
          .filter_map(|(key, missing)| missing.then_some(key))
          .collect_vec();
        return Err(missing_required_args::<TreetimeClockArgsRaw>(&missing));
      },
    };
    Ok(Self {
      tree,
      newick_dialect: newick_dialect.dialect(),
      dates: DatesInput {
        metadata: MetadataInput {
          path: metadata,
          id_columns: metadata_id_columns,
          delimiters: metadata_delimiters,
        },
        date_column,
        date_format,
      },
      clock_filter,
      allow_negative_rate,
      reroot: RerootTarget::unless_keep_root(keep_root, reroot.map(Into::into), outgroup),
      topology_order: clock_topology_order(topology_order)?,
      outputs: clock_outputs(output, output_selection, clock_output, output_rtt_chart),
    })
  }
}

#[deny(unused_variables)]
fn clock_outputs(
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
  ClockOutputSelectionArgs { output_selection }: ClockOutputSelectionArgs,
  ClockOutputArgs {
    output_clock_model,
    output_rtt,
  }: ClockOutputArgs,
  output_rtt_chart: Option<PathBuf>,
) -> OutputRequest {
  OutputRequest {
    command: CommandKind::Clock,
    alignment_count: 0,
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
      (OutputSelection::ClockModel, output_clock_model),
      (OutputSelection::Rtt, output_rtt),
      (OutputSelection::RttChart, output_rtt_chart),
    ],
  }
}

#[deny(unused_variables)]
fn clock_topology_order(
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

#[derive(Debug, Clone, SmartDefault, JsonSchema, Serialize, Deserialize)]
#[deser(skip_serializing_optionals)]
#[schemars(transform = skip_serializing_optionals)]
#[schemars(default, deny_unknown_fields)]
#[schemars(transform = schema_defaults::<Self>)]
#[deser(default, deny_unknown_fields)]
#[cfg_attr(feature = "clap", derive(clap::Parser))]
#[cfg_attr(feature = "clap", command(mut_args(value_name_from_hint)))]
#[schemars(rename = "ClockConfig", description = CLOCK_ABOUT)]
pub struct TreetimeClockArgsRaw {
  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(skip)]
  #[deser(skip)]
  pub config_args: ConfigArgs,

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
  pub metadata: MetadataArgs,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub metadata_id: MetadataIdArgs,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub date_column: DateColumnArgs,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub clock_filter: ClockFilterArgs,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub reroot: RerootArgs,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub rooting: RootingArgs,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub output: OutputCoreArgs,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub output_selection: ClockOutputSelectionArgs,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub clock_output: ClockOutputArgs,

  /// Path to output the root-to-tip regression chart.
  ///
  /// The chart plots the root-to-tip divergence of each dated sample against its date, with the least-squares fit, its extrapolation to the root, and the samples that the clock filter excluded. It states the clock rate, the root date and R² of the fit. The extension chooses the image format: `.svg`, `.png` or `.pdf`.
  #[cfg_attr(feature = "clap", clap(long, hide_short_help = true, value_hint = ValueHint::FilePath, help_heading = "Output"))]
  #[schemars(extend("x-path" = "output"))]
  pub output_rtt_chart: Option<PathBuf>,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub topology_order: TopologyOrderArgs,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(skip)]
  #[deser(skip)]
  pub jobs_args: JobsArgs,
}
