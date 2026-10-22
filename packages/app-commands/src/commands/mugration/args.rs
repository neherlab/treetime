use crate::commands::shared::ambiguity::{MugrationAmbiguityArgs, resolve_ambiguity_outputs};
use crate::commands::shared::config::ConfigArgs;
use crate::commands::shared::gtr_output::GtrOutputArgs;
use crate::commands::shared::jobs::JobsArgs;
use crate::commands::shared::metadata::{MetadataIdArgs, MetadataInput};
use crate::commands::shared::output_args::{
  MugrationOutputSelectionArgs, OutputCoreArgs, OutputRequest, TreeOutputPaths, selected_outputs,
};
use crate::commands::shared::required::missing_required_args;
use crate::commands::shared::topology_order_args::{TopologyOrder, TopologyOrderArgs};
use crate::commands::shared::tree_input::{TreeArgs, TreeDialectArgs};
use app_output::ambiguity_outputs::AmbiguityOutputs;
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
use treetime::gtr::infer_gtr::common::DEFAULT_PSEUDO_COUNT;
use treetime_io::nwk::NewickDialect;
use treetime_schema::{schema_defaults, skip_serializing_optionals};

#[derive(Debug, Clone)]
pub struct TreetimeMugrationArgs {
  pub(crate) tree: PathBuf,
  pub(crate) newick_dialect: NewickDialect,
  pub(crate) metadata: MetadataInput,
  pub(crate) attribute: String,
  pub(crate) weights: Option<PathBuf>,
  pub(crate) pc: f64,
  pub(crate) missing_data: String,
  pub(crate) missing_weights_threshold: f64,
  pub(crate) model_iterations: usize,
  pub(crate) sampling_bias_correction: Option<f64>,
  pub(crate) smooth_initial_pi: bool,
  pub(crate) filter_uninformative_root: bool,
  pub(crate) ambiguity_outputs: AmbiguityOutputs,
  pub(crate) topology_order: TopologyOrder,
  pub(crate) outputs: OutputRequest,
}

impl TryFrom<TreetimeMugrationArgsRaw> for TreetimeMugrationArgs {
  type Error = Report;

  #[deny(unused_variables)]
  #[expect(
    clippy::unneeded_field_pattern,
    reason = "the CLI reads its own settings before the conversion, and `..` would hide every setting that no run reads"
  )]
  fn try_from(raw: TreetimeMugrationArgsRaw) -> Result<Self, Report> {
    let TreetimeMugrationArgsRaw {
      config_args: _,
      tree: TreeArgs { tree },
      newick_dialect: TreeDialectArgs { newick_dialect },
      attribute,
      metadata,
      weights,
      metadata_id: MetadataIdArgs {
        metadata_id_columns,
        metadata_delimiters,
      },
      confidence,
      output_trait_probabilities,
      pc,
      missing_data,
      missing_weights_threshold,
      model_iterations,
      sampling_bias_correction,
      smooth_initial_pi,
      filter_uninformative_root,
      output_node_data,
      gtr_output,
      output_traits,
      output,
      output_selection,
      ambiguity: MugrationAmbiguityArgs {
        ambiguity_threshold,
        ambiguity_outputs,
      },
      topology_order,
      jobs_args: _,
    } = raw;
    let (tree, metadata, attribute) = match (tree, metadata, attribute) {
      (Some(tree), Some(metadata), Some(attribute)) => (tree, metadata, attribute),
      (tree, metadata, attribute) => {
        let missing = [
          ("tree", tree.is_none()),
          ("metadata", metadata.is_none()),
          ("attribute", attribute.is_none()),
        ]
        .into_iter()
        .filter_map(|(key, missing)| missing.then_some(key))
        .collect_vec();
        return Err(missing_required_args::<TreetimeMugrationArgsRaw>(&missing));
      },
    };
    Ok(Self {
      tree,
      newick_dialect: newick_dialect.dialect(),
      metadata: MetadataInput {
        path: metadata,
        id_columns: metadata_id_columns,
        delimiters: metadata_delimiters,
      },
      attribute,
      weights,
      pc,
      missing_data,
      missing_weights_threshold,
      model_iterations,
      sampling_bias_correction,
      smooth_initial_pi,
      filter_uninformative_root,
      ambiguity_outputs: resolve_ambiguity_outputs(ambiguity_threshold, &ambiguity_outputs)?,
      topology_order: mugration_topology_order(topology_order)?,
      outputs: mugration_outputs(
        output,
        output_selection,
        confidence,
        output_node_data,
        gtr_output,
        output_trait_probabilities,
        output_traits,
      ),
    })
  }
}

#[deny(unused_variables)]
fn mugration_outputs(
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
  MugrationOutputSelectionArgs { output_selection }: MugrationOutputSelectionArgs,
  confidence: bool,
  output_node_data: Option<PathBuf>,
  GtrOutputArgs { output_model }: GtrOutputArgs,
  output_trait_probabilities: Option<PathBuf>,
  output_traits: Option<PathBuf>,
) -> OutputRequest {
  OutputRequest {
    command: CommandKind::Mugration,
    alignment_count: 0,
    settings: mugration_output_settings(confidence),
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
      (OutputSelection::NodeData, output_node_data),
      (OutputSelection::Model, output_model),
      (OutputSelection::TraitProbabilities, output_trait_probabilities),
      (OutputSelection::Traits, output_traits),
    ],
  }
}

#[deny(unused_variables)]
fn mugration_topology_order(
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

pub(crate) fn mugration_output_settings(confidence: bool) -> OutputSettings {
  OutputSettings {
    trait_probabilities: confidence,
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
