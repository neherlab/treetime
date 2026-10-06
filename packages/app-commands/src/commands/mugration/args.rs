use crate::commands::shared::config::ConfigArgs;
use crate::commands::shared::metadata::MetadataIdArgs;
use crate::commands::shared::output_args::{MugrationOutputSelection, OutputCoreArgs};
use crate::commands::shared::required::missing_required_args;
use crate::commands::shared::topology_order_args::TopologyOrderArgs;
#[cfg(feature = "clap")]
use clap::ValueHint;
use deser::{Deserialize, Serialize};
use eyre::Report;
use schemars::JsonSchema;
use smart_default::SmartDefault;
use std::path::{Path, PathBuf};
use treetime_schema::{schema_defaults, skip_serializing_optionals};

#[derive(Debug, Clone)]
pub struct TreetimeMugrationArgs {
  pub(crate) tree: PathBuf,
  pub(crate) attribute: String,
  pub(crate) metadata: PathBuf,
  pub(crate) weights: Option<PathBuf>,
  pub(crate) metadata_id: MetadataIdArgs,
  pub(crate) output_confidence_csv: Option<PathBuf>,
  pub(crate) pc: Option<f64>,
  pub(crate) missing_data: String,
  pub(crate) missing_weights_threshold: f64,
  pub(crate) iterations: usize,
  pub(crate) sampling_bias_correction: Option<f64>,
  pub(crate) smooth_initial_pi: bool,
  pub(crate) filter_uninformative_root: bool,
  pub(crate) output_augur_node_data: Option<PathBuf>,
  pub(crate) output_gtr: Option<PathBuf>,
  pub(crate) output_traits_csv: Option<PathBuf>,
  pub(crate) output: OutputCoreArgs,
  pub(crate) output_selection: Vec<MugrationOutputSelection>,
  pub(crate) topology_order: TopologyOrderArgs,
}

impl TreetimeMugrationArgs {
  pub fn metadata(&self) -> &Path {
    &self.metadata
  }

  pub fn attribute(&self) -> &str {
    &self.attribute
  }
}

impl TryFrom<TreetimeMugrationArgsRaw> for TreetimeMugrationArgs {
  type Error = Report;

  fn try_from(raw: TreetimeMugrationArgsRaw) -> Result<Self, Report> {
    match (raw.tree, raw.metadata, raw.attribute) {
      (Some(tree), Some(metadata), Some(attribute)) => Ok(Self {
        tree,
        attribute,
        metadata,
        weights: raw.weights,
        metadata_id: raw.metadata_id,
        output_confidence_csv: raw.output_confidence_csv,
        pc: raw.pc,
        missing_data: raw.missing_data,
        missing_weights_threshold: raw.missing_weights_threshold,
        iterations: raw.iterations,
        sampling_bias_correction: raw.sampling_bias_correction,
        smooth_initial_pi: raw.smooth_initial_pi,
        filter_uninformative_root: raw.filter_uninformative_root,
        output_augur_node_data: raw.output_augur_node_data,
        output_gtr: raw.output_gtr,
        output_traits_csv: raw.output_traits_csv,
        output: raw.output,
        output_selection: raw.output_selection,
        topology_order: raw.topology_order,
      }),
      (tree, metadata, attribute) => {
        let mut missing = Vec::new();
        if tree.is_none() {
          missing.push("tree");
        }
        if metadata.is_none() {
          missing.push("metadata");
        }
        if attribute.is_none() {
          missing.push("attribute");
        }
        Err(missing_required_args::<TreetimeMugrationArgsRaw>(&missing))
      },
    }
  }
}

#[derive(Debug, Clone, SmartDefault, JsonSchema, Serialize, Deserialize)]
#[deser(skip_serializing_optionals)]
#[schemars(transform = skip_serializing_optionals)]
#[schemars(default, deny_unknown_fields)]
#[schemars(transform = schema_defaults::<Self>)]
#[deser(default, deny_unknown_fields)]
#[cfg_attr(feature = "clap", derive(clap::Parser))]
#[schemars(rename = "MugrationConfig")]
pub struct TreetimeMugrationArgsRaw {
  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(skip)]
  #[deser(skip)]
  pub config_args: ConfigArgs,

  /// Tree in Newick format.
  #[cfg_attr(feature = "clap", clap(long, short = 't', help_heading = "Input data"))]
  #[cfg_attr(feature = "clap", clap(value_hint = ValueHint::FilePath))]
  #[schemars(extend("x-path" = "input"))]
  pub tree: Option<PathBuf>,

  /// Attribute to reconstruct, e.g. country
  #[cfg_attr(feature = "clap", clap(long, help_heading = "Input data"))]
  pub attribute: Option<String>,

  /// CSV or TSV file with discrete characters. #name,country,continent taxon1,micronesia,oceania ...
  #[cfg_attr(
    feature = "clap",
    clap(
      long = "metadata",
      visible_alias = "states",
      short = 's',
      help_heading = "Input data"
    )
  )]
  #[cfg_attr(feature = "clap", clap(value_hint = ValueHint::FilePath))]
  #[schemars(extend("x-path" = "input"))]
  pub metadata: Option<PathBuf>,

  /// CSV or TSV file with probabilities of that a randomly sampled sequence at equilibrium has a particular state. E.g. population of different continents or countries. E.g.: #country,weight micronesia,0.1 ...
  #[cfg_attr(feature = "clap", clap(long, short = 'w', help_heading = "Input data"))]
  #[cfg_attr(feature = "clap", clap(value_hint = ValueHint::FilePath))]
  #[schemars(extend("x-path" = "input"))]
  pub weights: Option<PathBuf>,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub metadata_id: MetadataIdArgs,

  /// Path to output state-probability-profile CSV.
  ///
  /// Takes precedence over paths configured with `--output-all` and `--output-selection`.
  #[cfg_attr(feature = "clap", clap(long, visible_alias = "confidence", value_hint = ValueHint::FilePath, help_heading = "Output"))]
  #[schemars(extend("x-path" = "output"))]
  pub output_confidence_csv: Option<PathBuf>,

  /// Pseudo-counts. Higher numbers results in 'flatter' models. Default: 1.0.
  #[cfg_attr(feature = "clap", clap(long, help_heading = "Substitution model"))]
  pub pc: Option<f64>,

  /// String indicating missing data
  #[cfg_attr(feature = "clap", clap(long, default_value = "?", help_heading = "Input data"))]
  #[default(_code = r#""?".to_owned()"#)]
  pub missing_data: String,

  /// Portion of attribute values that is allowed to not have weights in the weights file
  #[cfg_attr(
    feature = "clap",
    clap(long, default_value_t = 0.5, help_heading = "Substitution model")
  )]
  #[default = 0.5]
  pub missing_weights_threshold: f64,

  /// Number of iterations for GTR model refinement from data.
  #[cfg_attr(
    feature = "clap",
    clap(long, default_value_t = 5, help_heading = "Substitution model")
  )]
  #[default = 5]
  pub iterations: usize,

  /// Rough estimate of how many more events would have been observed if sequences represented an
  /// even sample.
  #[cfg_attr(feature = "clap", clap(long, help_heading = "Substitution model"))]
  pub sampling_bias_correction: Option<f64>,

  /// Smooth the initial equilibrium frequencies with the pseudo-count before the first
  /// reconstruction pass.
  ///
  /// Off by default (TreeTime v0 builds the initial model from raw frequencies and applies the
  /// pseudo-count only as infer_gtr regularization). Enabling this flattens the prior for the first
  /// pass; it only affects weighted models.
  #[cfg_attr(feature = "clap", clap(long, help_heading = "Substitution model"))]
  pub smooth_initial_pi: bool,

  /// Exclude near-uniform root positions from the equilibrium-frequency prior.
  ///
  /// Off by default (TreeTime v0 always folds the root's most-likely state into the prior). Enabling
  /// this drops root positions whose posterior carries no phylogenetic signal, removing a
  /// state-order-dependent bias at ambiguous roots.
  #[cfg_attr(feature = "clap", clap(long, help_heading = "Substitution model"))]
  pub filter_uninformative_root: bool,

  /// Path to output augur-compatible node data JSON.
  ///
  /// Contains per-node discrete trait assignments, confidence profiles, entropy,
  /// the inferred substitution model, and branch state-change labels. The output
  /// is compatible with augur export v2 --node-data for Nextstrain pipeline
  /// integration.
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

  /// Path to output traits CSV.
  ///
  /// Takes precedence over paths configured with `--output-all` and `--output-selection`.
  #[cfg_attr(feature = "clap", clap(long, value_hint = ValueHint::FilePath, help_heading = "Output"))]
  #[schemars(extend("x-path" = "output"))]
  pub output_traits_csv: Option<PathBuf>,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub output: OutputCoreArgs,

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
  pub output_selection: Vec<MugrationOutputSelection>,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub topology_order: TopologyOrderArgs,
}
