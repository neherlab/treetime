use crate::commands::shared::clock_filter::ClockFilterArgs;
use crate::commands::shared::config::ConfigArgs;
use crate::commands::shared::jobs::JobsArgs;
use crate::commands::shared::metadata::{DateColumnArgs, MetadataIdArgs};
use crate::commands::shared::output_args::{ClockOutputArgs, ClockOutputSelectionArgs, OutputCoreArgs};
use crate::commands::shared::required::missing_required_args;
use crate::commands::shared::reroot::{RerootArgs, RootingArgs};
use crate::commands::shared::topology_order_args::TopologyOrderArgs;
use crate::commands::shared::tree_input::{TreeArgs, TreeDialectArgs};
#[cfg(feature = "clap")]
use clap::ValueHint;
use derive_more::{Deref, DerefMut};
use deser::{Deserialize, Serialize};
use eyre::Report;
use itertools::Itertools;
use schemars::JsonSchema;
use smart_default::SmartDefault;
use std::path::{Path, PathBuf};
use treetime_schema::{schema_defaults, skip_serializing_optionals};

#[derive(Debug, Clone, Deref, DerefMut)]
pub struct TreetimeClockArgs {
  #[deref]
  #[deref_mut]
  pub(crate) raw: TreetimeClockArgsRaw,
  pub(crate) tree: PathBuf,
  pub(crate) metadata: PathBuf,
}

impl TreetimeClockArgs {
  pub fn metadata(&self) -> &Path {
    &self.metadata
  }
}

impl TryFrom<TreetimeClockArgsRaw> for TreetimeClockArgs {
  type Error = Report;

  fn try_from(raw: TreetimeClockArgsRaw) -> Result<Self, Report> {
    let (tree, metadata) = match (raw.tree.tree.clone(), raw.metadata.clone()) {
      (Some(tree), Some(metadata)) => (tree, metadata),
      (tree, metadata) => {
        let missing = [("tree", tree.is_none()), ("metadata", metadata.is_none())]
          .into_iter()
          .filter_map(|(key, missing)| missing.then_some(key))
          .collect_vec();
        return Err(missing_required_args::<TreetimeClockArgsRaw>(&missing));
      },
    };
    Ok(Self { raw, tree, metadata })
  }
}

#[derive(Debug, Clone, SmartDefault, JsonSchema, Serialize, Deserialize)]
#[deser(skip_serializing_optionals)]
#[schemars(transform = skip_serializing_optionals)]
#[schemars(default, deny_unknown_fields)]
#[schemars(transform = schema_defaults::<Self>)]
#[deser(default, deny_unknown_fields)]
#[cfg_attr(feature = "clap", derive(clap::Parser))]
#[schemars(rename = "ClockConfig")]
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
  pub tree_dialect: TreeDialectArgs,

  /// Metadata table (CSV or TSV) with one row per sample
  ///
  /// `--metadata-id-columns` names the column that links a row to a tree tip. The other columns hold the data of the
  /// sample, such as its sampling date or its traits.
  #[cfg_attr(
    feature = "clap",
    clap(long = "metadata", visible_alias = "dates", short = 'd', help_heading = "Input data")
  )]
  #[cfg_attr(feature = "clap", clap(value_hint = ValueHint::FilePath))]
  #[schemars(extend("x-path" = "input"))]
  pub metadata: Option<PathBuf>,

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
  pub clock_output: ClockOutputArgs,

  /// Path to output root-to-tip regression chart in SVG format.
  ///
  /// Takes precedence over paths configured with `--output-all` and `--output-selection`.
  #[cfg_attr(feature = "clap", clap(long, value_hint = ValueHint::FilePath, help_heading = "Output"))]
  #[schemars(extend("x-path" = "output"))]
  pub output_clock_chart_svg: Option<PathBuf>,

  /// Path to output root-to-tip regression chart in PNG format.
  ///
  /// Takes precedence over paths configured with `--output-all` and `--output-selection`.
  #[cfg_attr(feature = "clap", clap(long, value_hint = ValueHint::FilePath, help_heading = "Output"))]
  #[schemars(extend("x-path" = "output"))]
  pub output_clock_chart_png: Option<PathBuf>,

  /// Path to output root-to-tip regression chart in PDF format.
  ///
  /// Takes precedence over paths configured with `--output-all` and `--output-selection`.
  #[cfg_attr(feature = "clap", clap(long, value_hint = ValueHint::FilePath, help_heading = "Output"))]
  #[schemars(extend("x-path" = "output"))]
  pub output_clock_chart_pdf: Option<PathBuf>,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub output_selection: ClockOutputSelectionArgs,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub topology_order: TopologyOrderArgs,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(skip)]
  #[deser(skip)]
  pub jobs_args: JobsArgs,
}
