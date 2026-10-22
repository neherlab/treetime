use crate::commands::shared::ambiguity::MugrationAmbiguityArgs;
use crate::commands::shared::config::ConfigArgs;
use crate::commands::shared::gtr_output::GtrOutputArgs;
use crate::commands::shared::jobs::JobsArgs;
use crate::commands::shared::metadata::MetadataIdArgs;
use crate::commands::shared::output_args::{MugrationOutputSelectionArgs, OutputCoreArgs};
use crate::commands::shared::required::missing_required_args;
use crate::commands::shared::topology_order_args::TopologyOrderArgs;
use crate::commands::shared::tree_input::{TreeArgs, TreeDialectArgs};
use app_output::ambiguity_outputs::AmbiguityOutputs;
use app_output::output_conditions::OutputSettings;
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
  pub(crate) ambiguity_outputs: AmbiguityOutputs,
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
    match (raw.tree.tree.clone(), raw.metadata.clone(), raw.attribute.clone()) {
      (Some(tree), Some(metadata), Some(attribute)) => Ok(Self {
        ambiguity_outputs: raw.ambiguity.resolve()?,
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

pub(crate) fn mugration_output_settings(raw: &TreetimeMugrationArgsRaw) -> OutputSettings {
  OutputSettings {
    trait_probabilities: raw.confidence,
    ..OutputSettings::default()
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

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub tree: TreeArgs,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub newick_dialect: TreeDialectArgs,

  /// Attribute to reconstruct, e.g. country
  #[cfg_attr(feature = "clap", clap(long, help_heading = "Input data"))]
  pub attribute: Option<String>,

  /// Metadata table (CSV or TSV) with one row per sample
  ///
  /// `--metadata-id-columns` names the column that links a row to a tree tip. The other columns hold the data of the
  /// sample, such as its sampling date or its traits.
  #[cfg_attr(feature = "clap", clap(long, short = 'm', help_heading = "Input data"))]
  #[cfg_attr(feature = "clap", clap(value_hint = ValueHint::FilePath))]
  #[schemars(extend("x-path" = "input"))]
  pub metadata: Option<PathBuf>,

  /// CSV or TSV file with probabilities of that a randomly sampled sequence at equilibrium has a particular state. E.g. population of different continents or countries. E.g.: #country,weight micronesia,0.1 ...
  #[cfg_attr(feature = "clap", clap(long, help_heading = "Input data"))]
  #[cfg_attr(feature = "clap", clap(value_hint = ValueHint::FilePath))]
  #[schemars(extend("x-path" = "input"))]
  pub weights: Option<PathBuf>,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub metadata_id: MetadataIdArgs,

  /// Write the probability of every state at every node
  ///
  /// Adds the `trait-probabilities` table to the outputs of `--output-all`, also when `--output-selection` lists
  /// other outputs. Without `--output-all`, give the file with `--output-trait-probabilities`.
  #[cfg_attr(feature = "clap", clap(long, help_heading = "Output"))]
  pub confidence: bool,

  /// Path to output the state probabilities of every node as a table.
  ///
  /// Takes precedence over paths configured with `--output-all` and `--output-selection`.
  #[cfg_attr(feature = "clap", clap(long, value_hint = ValueHint::FilePath, help_heading = "Output"))]
  #[schemars(extend("x-path" = "output"))]
  pub output_trait_probabilities: Option<PathBuf>,

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
  pub model_iterations: usize,

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
  pub output_node_data: Option<PathBuf>,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub gtr_output: GtrOutputArgs,

  /// Path to output the inferred trait state of every node as a table.
  ///
  /// Takes precedence over paths configured with `--output-all` and `--output-selection`.
  #[cfg_attr(feature = "clap", clap(long, value_hint = ValueHint::FilePath, help_heading = "Output"))]
  #[schemars(extend("x-path" = "output"))]
  pub output_traits: Option<PathBuf>,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub output: OutputCoreArgs,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub output_selection: MugrationOutputSelectionArgs,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub ambiguity: MugrationAmbiguityArgs,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub topology_order: TopologyOrderArgs,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(skip)]
  #[deser(skip)]
  pub jobs_args: JobsArgs,
}
