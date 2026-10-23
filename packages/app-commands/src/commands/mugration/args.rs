use crate::commands::shared::ambiguity::{MugrationAmbiguityArgs, resolve_ambiguity_outputs};
use crate::commands::shared::config::ConfigArgs;
use crate::commands::shared::gtr_output::GtrOutputArgs;
use crate::commands::shared::jobs::JobsArgs;
use crate::commands::shared::metadata::{MetadataArgs, MetadataIdArgs, MetadataInput};
use crate::commands::shared::output_args::{
  MugrationOutputSelectionArgs, OutputCoreArgs, OutputRequest, TreeOutputPaths, selected_outputs,
};
use crate::commands::shared::required::missing_required_args;
use crate::commands::shared::topology_order_args::{TopologyOrder, TopologyOrderArgs};
use crate::commands::shared::tree_input::{TreeArgs, TreeDialectArgs};
#[cfg(feature = "clap")]
use crate::commands::shared::value_names::value_name_from_hint;
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

pub const MUGRATION_ABOUT: &str =
  "Reconstruct the ancestral states of a discrete trait, such as location or host, and fit a model of its changes";

#[derive(Debug, Clone)]
pub struct TreetimeMugrationArgs {
  pub(crate) tree: PathBuf,
  pub(crate) newick_dialect: NewickDialect,
  pub(crate) metadata: MetadataInput,
  pub(crate) attribute: String,
  pub(crate) weights: Option<PathBuf>,
  pub(crate) lat_longs: Option<PathBuf>,
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
      metadata: MetadataArgs { metadata },
      weights,
      lat_longs,
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
      lat_longs,
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
#[cfg_attr(feature = "clap", command(mut_args(value_name_from_hint)))]
#[schemars(rename = "MugrationConfig", description = MUGRATION_ABOUT)]
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

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub metadata: MetadataArgs,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub metadata_id: MetadataIdArgs,

  /// Metadata column with the discrete trait to reconstruct, for example `country`
  #[cfg_attr(feature = "clap", clap(long, value_name = "COLUMN", help_heading = "Input data"))]
  pub attribute: Option<String>,

  /// Table (CSV or TSV) with the equilibrium probability of each state
  ///
  /// The weight of a state is the probability that a randomly sampled sequence at equilibrium has the state, for
  /// example the population of each country. The first column holds the state and the second its weight, such as
  /// the row `micronesia,0.1` under the header `country,weight`.
  #[cfg_attr(feature = "clap", clap(long, help_heading = "Input data"))]
  #[cfg_attr(feature = "clap", clap(value_hint = ValueHint::FilePath))]
  #[schemars(extend("x-path" = "input"))]
  pub weights: Option<PathBuf>,
  /// Table with the coordinates of the trait values, which give the Auspice tree a map
  ///
  /// The format of augur's `--lat-longs`: one location per line, `<resolution> <name> <latitude> <longitude>`, with
  /// the fields separated by tabs, or by white space in a line without a tab. `<resolution>` is the `--attribute`
  /// name, and `<name>` a value of the trait; both match in lower case. A later line for the same location replaces
  /// an earlier one, and lines that start with `#` are skipped. Latitudes are in -90 to 90 degrees, longitudes in
  /// -180 to 180 degrees.
  #[cfg_attr(
    feature = "clap",
    clap(long, hide_short_help = true, value_hint = ValueHint::FilePath, help_heading = "Input data")
  )]
  #[schemars(extend("x-path" = "input"))]
  pub lat_longs: Option<PathBuf>,

  /// Value of the trait column that marks missing data
  #[cfg_attr(feature = "clap", clap(long, default_value = "?", help_heading = "Input data"))]
  #[default(_code = r#""?".to_owned()"#)]
  pub missing_data: String,

  /// Pseudo-count of the substitution model fit. Higher values give flatter models.
  #[cfg_attr(
    feature = "clap",
    clap(long, default_value_t = DEFAULT_PSEUDO_COUNT, help_heading = "Substitution model")
  )]
  #[default(DEFAULT_PSEUDO_COUNT)]
  pub pc: f64,

  /// Largest fraction of the trait values that may have no weight in the `--weights` file
  #[cfg_attr(
    feature = "clap",
    clap(
      long,
      hide_short_help = true,
      default_value_t = 0.5,
      help_heading = "Substitution model"
    )
  )]
  #[default = 0.5]
  pub missing_weights_threshold: f64,

  /// Number of iterations of the refinement of the substitution model from the data
  #[cfg_attr(
    feature = "clap",
    clap(
      long,
      value_name = "N",
      hide_short_help = true,
      default_value_t = 5,
      help_heading = "Substitution model"
    )
  )]
  #[default = 5]
  pub model_iterations: usize,

  /// Correction for biased sampling: a rough estimate of how many more state changes the tree would show if the
  /// samples were an even sample of the population
  ///
  /// The run multiplies the fitted rate of the substitution model by this factor, for example 2.5 for an estimate of
  /// 2.5 times more changes. Without this setting, the run applies no correction.
  #[cfg_attr(feature = "clap", clap(long, help_heading = "Substitution model"))]
  pub sampling_bias_correction: Option<f64>,

  /// Smooth the initial equilibrium frequencies with the pseudo-count before the first reconstruction pass
  ///
  /// By default, the initial model takes the frequencies of `--weights` as they are, and the pseudo-count applies only
  /// to the fits of the model. This setting flattens the prior of the first pass. It affects only models with
  /// `--weights`.
  #[cfg_attr(
    feature = "clap",
    clap(long, hide_short_help = true, help_heading = "Substitution model")
  )]
  pub smooth_initial_pi: bool,

  /// Exclude a near-uniform root posterior from the prior of the equilibrium frequencies
  ///
  /// By default, the most likely state of the root always enters the prior. This setting leaves out a root whose
  /// posterior carries no phylogenetic signal, which removes a bias at ambiguous roots that depends on the order of
  /// the states.
  #[cfg_attr(
    feature = "clap",
    clap(long, hide_short_help = true, help_heading = "Substitution model")
  )]
  pub filter_uninformative_root: bool,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub output: OutputCoreArgs,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub output_selection: MugrationOutputSelectionArgs,

  /// Write the probability of every state at every node
  ///
  /// Adds the `trait-probabilities` table to the outputs of `--output-all`, also when `--output-selection` lists
  /// other outputs. Without `--output-all`, give the file with `--output-trait-probabilities`.
  #[cfg_attr(feature = "clap", clap(long, help_heading = "Output"))]
  pub confidence: bool,

  /// Path to output augur-compatible node data JSON.
  ///
  /// Contains per-node discrete trait assignments, confidence profiles, entropy, the inferred substitution model, and
  /// branch state-change labels. The output is compatible with augur export v2 --node-data for Nextstrain pipeline
  /// integration.
  #[cfg_attr(feature = "clap", clap(long, hide_short_help = true, value_hint = ValueHint::FilePath, help_heading = "Output"))]
  #[schemars(extend("x-path" = "output"))]
  pub output_node_data: Option<PathBuf>,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub gtr_output: GtrOutputArgs,

  /// Path to output the inferred trait state of every node as a table.
  #[cfg_attr(feature = "clap", clap(long, hide_short_help = true, value_hint = ValueHint::FilePath, help_heading = "Output"))]
  #[schemars(extend("x-path" = "output"))]
  pub output_traits: Option<PathBuf>,

  /// Path to output the state probabilities of every node as a table.
  #[cfg_attr(feature = "clap", clap(long, hide_short_help = true, value_hint = ValueHint::FilePath, help_heading = "Output"))]
  #[schemars(extend("x-path" = "output"))]
  pub output_trait_probabilities: Option<PathBuf>,

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
