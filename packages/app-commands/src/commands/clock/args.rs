use crate::commands::shared::alignment::AlignmentArgs;
use crate::commands::shared::branch_length_mode::BranchLengthModeCli;
use crate::commands::shared::config::ConfigArgs;
use crate::commands::shared::metadata::{DateColumnArgs, MetadataIdArgs};
use crate::commands::shared::method_anc::MethodAncestralCli;
use crate::commands::shared::model::ModelArgs;
use crate::commands::shared::output_args::{ClockOutputSelection, OutputCoreArgs};
use crate::commands::shared::required::missing_required_args;
use crate::commands::shared::reroot::RerootArgs;
use crate::commands::shared::topology_order_args::TopologyOrderArgs;
use crate::commands::shared::tree_input::TreeDialectArgs;
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
    let (tree, metadata) = match (raw.tree.clone(), raw.metadata.clone()) {
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

  /// Only for vcf input: fasta file of the sequence the VCF was mapped to.
  #[cfg_attr(feature = "clap", clap(long, short = 'r', help_heading = "Input data"))]
  #[cfg_attr(feature = "clap", clap(value_hint = ValueHint::FilePath))]
  #[schemars(extend("x-path" = "input"))]
  pub vcf_reference: Option<PathBuf>,

  /// CSV/TSV file with metadata including sampling dates
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
  pub model_args: ModelArgs,

  /// If set to 'input', the provided branch length will be used without modification. Note that branch lengths optimized by treetime are only accurate at short evolutionary distances.
  #[cfg_attr(feature = "clap", clap(long, value_enum, default_value_t = BranchLengthModeCli::default(), help_heading = "Branch lengths"))]
  pub branch_length_mode: BranchLengthModeCli,

  /// Method used for reconstructing ancestral sequences
  #[cfg_attr(feature = "clap", clap(long, value_enum, default_value_t = MethodAncestralCli::default(), help_heading = "Ancestral reconstruction"))]
  pub method_anc: MethodAncestralCli,

  /// ignore tips that don't follow a loose clock, 'clock-filter=number of interquartile ranges from regression'. Default=3.0, set to 0 to switch off.
  #[cfg_attr(
    feature = "clap",
    clap(long, default_value = "3.0", help_heading = "Molecular clock")
  )]
  #[default = 3.0]
  pub clock_filter: f64,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub reroot: RerootArgs,

  /// don't reroot the tree. Otherwise, reroot to minimize the residual of the regression of
  /// root-to-tip distance and sampling time
  #[cfg_attr(feature = "clap", clap(long, conflicts_with_all = ["reroot", "reroot_tips"], help_heading = "Rooting"))]
  pub keep_root: bool,

  #[cfg_attr(feature = "clap", clap(long, help_heading = "Branch lengths"))]
  pub prune_short: bool,

  /// By default, rates are forced to be positive. For trees with little temporal signal it is advisable to remove this restriction to achieve essentially mid-point rooting.
  #[cfg_attr(feature = "clap", clap(long, help_heading = "Molecular clock"))]
  pub allow_negative_rate: bool,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub output: OutputCoreArgs,

  /// Path to output clock model JSON.
  ///
  /// Takes precedence over paths configured with `--output-all` and `--output-selection`.
  #[cfg_attr(feature = "clap", clap(long, value_hint = ValueHint::FilePath, help_heading = "Output"))]
  #[schemars(extend("x-path" = "output"))]
  pub output_clock_model: Option<PathBuf>,

  /// Path to output clock regression CSV.
  ///
  /// Takes precedence over paths configured with `--output-all` and `--output-selection`.
  #[cfg_attr(feature = "clap", clap(long, value_hint = ValueHint::FilePath, help_heading = "Output"))]
  #[schemars(extend("x-path" = "output"))]
  pub output_clock_csv: Option<PathBuf>,

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
  pub output_selection: Vec<ClockOutputSelection>,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub topology_order: TopologyOrderArgs,

  /// Method for clock filter outlier detection (not yet implemented)
  #[cfg_attr(feature = "clap", clap(long, hide = true, help_heading = "Molecular clock"))]
  pub clock_filter_method: Option<String>,

  /// Filename to save root-to-tip regression plot (not yet implemented)
  #[cfg_attr(feature = "clap", clap(long, hide = true, help_heading = "Plots"))]
  #[cfg_attr(feature = "clap", clap(value_hint = ValueHint::FilePath))]
  #[schemars(extend("x-path" = "output"))]
  pub plot_rtt: Option<PathBuf>,

  /// Prune clock outlier tips from the tree (not yet implemented)
  #[cfg_attr(feature = "clap", clap(long, hide = true, help_heading = "Molecular clock"))]
  pub prune_outliers: bool,
}
