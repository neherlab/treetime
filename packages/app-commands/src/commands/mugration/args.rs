use crate::commands::shared::config::ConfigArgs;
use crate::commands::shared::metadata::MetadataIdArgs;
use crate::commands::shared::output_args::{MugrationOutputSelection, OutputCoreArgs};
use crate::commands::shared::required::missing_required_args;
use crate::commands::shared::topology_order_args::TopologyOrderArgs;
use crate::commands::shared::tree_input::TreeDialectArgs;
#[cfg(feature = "clap")]
use clap::ValueHint;
use derive_more::{Deref, DerefMut};
use deser::{Deserialize, Serialize};
use eyre::Report;
use schemars::JsonSchema;
use smart_default::SmartDefault;
use std::path::{Path, PathBuf};
use treetime::gtr::infer_gtr::common::DEFAULT_PSEUDO_COUNT;
use treetime_schema::{schema_defaults, skip_serializing_optionals};

#[derive(Debug, Clone, Deref, DerefMut)]
pub struct TreetimeMugrationArgs {
  #[deref]
  #[deref_mut]
  pub(crate) raw: TreetimeMugrationArgsRaw,
  pub(crate) tree: PathBuf,
  pub(crate) metadata: PathBuf,
  pub(crate) attribute: String,
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
    match (raw.tree.clone(), raw.metadata.clone(), raw.attribute.clone()) {
      (Some(tree), Some(metadata), Some(attribute)) => Ok(Self {
        raw,
        tree,
        metadata,
        attribute,
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

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub tree_dialect: TreeDialectArgs,

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

  /// Pseudo-count of the substitution model fit. Higher values give flatter models.
  #[cfg_attr(
    feature = "clap",
    clap(long, default_value_t = DEFAULT_PSEUDO_COUNT, help_heading = "Substitution model")
  )]
  #[default(DEFAULT_PSEUDO_COUNT)]
  pub pc: f64,

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
