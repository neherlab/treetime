use crate::commands::ancestral::aa_model::AaModelName;
use crate::commands::shared::alignment::AlignmentArgs;
use crate::commands::shared::alphabet::AlphabetArgs;
use crate::commands::shared::ambiguity::{AncestralAmbiguityArgs, resolve_ambiguity_outputs};
use crate::commands::shared::config::ConfigArgs;
use crate::commands::shared::gap_fill::GapFillArgs;
use crate::commands::shared::gtr_output::GtrOutputArgs;
use crate::commands::shared::jobs::JobsArgs;
use crate::commands::shared::method_anc::{
  MarginalReconstruction, MethodAncArgs, MethodAncestralCli, SampleFromProfileArgs, check_model_with_method_anc,
};
use crate::commands::shared::model::{GtrModelNameCli, ModelArgs, resolve_model};
use crate::commands::shared::output_args::{
  AncestralOutputSelectionArgs, OutputCoreArgs, OutputRequest, TreeOutputPaths, selected_outputs,
};
use crate::commands::shared::reconstructed_fasta::ReconstructedFastaOutputArgs;
use crate::commands::shared::representation::{RepresentationArgs, RepresentationCli, representation_dense};
use crate::commands::shared::required::missing_required_args;
use crate::commands::shared::seed::SeedArgs;
use crate::commands::shared::sequence_inputs::{MissingSequenceArgs, SequenceInputs, missing_leaf_limit};
use crate::commands::shared::tip_sequences::{ImputeMissingDataArgs, TipSequenceArgs};
use crate::commands::shared::topology_order_args::{TopologyOrder, TopologyOrderArgs};
use crate::commands::shared::tree_input::{TreeArgs, TreeDialectArgs};
use crate::commands::shared::zero_based::ZeroBasedArgs;
use app_output::ambiguity_outputs::AmbiguityOutputs;
use app_output::command_outputs::CommandKind;
use app_output::output_conditions::OutputSettings;
use app_output::output_plan::OutputSelection;
#[cfg(feature = "clap")]
use clap::ValueHint;
use deser::{Deserialize, Serialize};
use eyre::Report;
use schemars::JsonSchema;
use smart_default::SmartDefault;
use std::path::{Path, PathBuf};
use treetime::ancestral::params::AncestralParams;
use treetime::ancestral::params::MethodAncestral;
use treetime::gtr::get_gtr::GtrModel;
use treetime::partition::marginal::posterior_states::AmbiguityThreshold;
use treetime::partition::marginal::sample::SampleMode;
use treetime::seq::mutation::PositionBase;
use treetime_io::nwk::NewickDialect;
use treetime_schema::{schema_defaults, skip_serializing_optionals};

pub fn ancestral_params(args: &TreetimeAncestralArgs, seed: u64) -> AncestralParams {
  AncestralParams {
    method: args.method_anc,
    model: args.model.clone(),
    dense: args.dense,
    include_leaves: args.include_leaves,
    report_ambiguous: args.report_ambiguous,
    impute_missing_data: args.impute_missing_data,
    seed,
    sample_from_profile: args.sample_from_profile,
  }
}

#[derive(Debug, Clone)]
pub struct TreetimeAncestralArgs {
  pub(crate) tree: PathBuf,
  pub(crate) newick_dialect: NewickDialect,
  pub(crate) sequences: SequenceInputs,
  pub(crate) allow_missing_sequences: bool,
  pub(crate) positions: PositionBase,
  pub(crate) method_anc: MethodAncestral,
  pub(crate) model: GtrModel,
  pub(crate) dense: Option<bool>,
  pub(crate) include_leaves: bool,
  pub(crate) report_ambiguous: bool,
  pub(crate) impute_missing_data: bool,
  pub(crate) sample_from_profile: SampleMode,
  pub(crate) seed: Option<u64>,
  pub(crate) translations: Option<String>,
  pub(crate) cds_selection: Vec<String>,
  pub(crate) annotation: Option<PathBuf>,
  pub(crate) aa_root_sequence: Option<PathBuf>,
  pub(crate) aa_model: AaModelName,
  pub(crate) ambiguity_outputs: AmbiguityOutputs,
  pub(crate) topology_order: TopologyOrder,
  pub(crate) outputs: OutputRequest,
}

impl TreetimeAncestralArgs {
  pub fn tree(&self) -> &Path {
    &self.tree
  }

  pub fn alignments(&self) -> &[PathBuf] {
    &self.sequences.alignments
  }
}

impl TryFrom<TreetimeAncestralArgsRaw> for TreetimeAncestralArgs {
  type Error = Report;

  #[deny(unused_variables)]
  #[expect(
    clippy::unneeded_field_pattern,
    reason = "the CLI reads its own settings before the conversion, and `..` would hide every setting that no run reads"
  )]
  fn try_from(raw: TreetimeAncestralArgsRaw) -> Result<Self, Report> {
    let TreetimeAncestralArgsRaw {
      config_args: _,
      alignment: AlignmentArgs { alignment },
      tree: TreeArgs { tree },
      newick_dialect: TreeDialectArgs { newick_dialect },
      alphabet_args,
      model_args: ModelArgs { model, model_params },
      method_anc: MethodAncArgs { method_anc },
      representation: RepresentationArgs { representation },
      sample_from_profile: SampleFromProfileArgs { sample_from_profile },
      gap_fill_args,
      zero_based: ZeroBasedArgs { zero_based },
      tip_sequences:
        TipSequenceArgs {
          include_leaves,
          reconstruct_tip_states,
          report_ambiguous,
        },
      impute_missing_data: ImputeMissingDataArgs { impute_missing_data },
      missing_sequences: MissingSequenceArgs {
        allow_missing_sequences,
      },
      output_node_data,
      gtr_output,
      reconstructed_fasta_output,
      translations,
      cds_selection,
      annotation,
      aa_root_sequence,
      aa_model,
      output_translations,
      output,
      output_selection,
      ambiguity: AncestralAmbiguityArgs {
        ambiguity_threshold,
        ambiguity_outputs,
      },
      topology_order,
      seed_args: SeedArgs { seed },
      jobs_args: _,
    } = raw;
    let tree = tree.ok_or_else(|| missing_required_args::<TreetimeAncestralArgsRaw>(&["tree"]))?;
    check_dependent_settings(&DependentSettings {
      method_anc,
      model,
      model_params: &model_params,
      translations: translations.is_some(),
      representation,
      ambiguity_threshold,
    })?;
    Ok(Self {
      tree,
      newick_dialect: newick_dialect.dialect(),
      outputs: ancestral_outputs(
        output,
        output_selection,
        alignment.len(),
        output_node_data,
        gtr_output,
        reconstructed_fasta_output,
        output_translations,
      ),
      sequences: ancestral_sequence_inputs(alignment, alphabet_args, gap_fill_args, allow_missing_sequences),
      allow_missing_sequences,
      positions: PositionBase::from_zero_based(zero_based),
      method_anc: method_anc.into(),
      model: resolve_model(model, &model_params)?,
      dense: representation_dense(representation),
      include_leaves: include_leaves || reconstruct_tip_states,
      report_ambiguous,
      impute_missing_data: impute_missing_data || reconstruct_tip_states,
      sample_from_profile: sample_from_profile.into(),
      seed,
      translations,
      cds_selection,
      annotation,
      aa_root_sequence,
      aa_model,
      ambiguity_outputs: resolve_ambiguity_outputs(ambiguity_threshold, &ambiguity_outputs)?,
      topology_order: ancestral_topology_order(topology_order)?,
    })
  }
}

#[deny(unused_variables)]
fn ancestral_outputs(
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
  AncestralOutputSelectionArgs { output_selection }: AncestralOutputSelectionArgs,
  alignment_count: usize,
  output_node_data: Option<PathBuf>,
  GtrOutputArgs { output_model }: GtrOutputArgs,
  ReconstructedFastaOutputArgs { output_sequences }: ReconstructedFastaOutputArgs,
  output_translations: Option<PathBuf>,
) -> OutputRequest {
  OutputRequest {
    command: CommandKind::Ancestral,
    alignment_count,
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
      (OutputSelection::NodeData, output_node_data),
      (OutputSelection::Model, output_model),
      (OutputSelection::Sequences, output_sequences),
      (OutputSelection::Translations, output_translations),
    ],
  }
}

#[deny(unused_variables)]
fn ancestral_sequence_inputs(
  alignment: Vec<PathBuf>,
  AlphabetArgs { alphabet }: AlphabetArgs,
  GapFillArgs { gap_fill }: GapFillArgs,
  allow_missing_sequences: bool,
) -> SequenceInputs {
  SequenceInputs {
    alignments: alignment,
    alphabet: alphabet.map(Into::into),
    gap_fill: gap_fill.into(),
    missing_leaves: missing_leaf_limit(allow_missing_sequences),
  }
}

#[deny(unused_variables)]
fn ancestral_topology_order(
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

struct DependentSettings<'a> {
  method_anc: MethodAncestralCli,
  model: Option<GtrModelNameCli>,
  model_params: &'a [String],
  translations: bool,
  representation: Option<RepresentationCli>,
  ambiguity_threshold: Option<AmbiguityThreshold>,
}

fn check_dependent_settings(settings: &DependentSettings<'_>) -> Result<(), Report> {
  check_model_with_method_anc(settings.method_anc, settings.model, settings.model_params)?;
  let marginal = MarginalReconstruction {
    runs: settings.method_anc.is_marginal() || settings.translations,
    condition: "--method-anc marginal, or the amino-acid reconstruction of --translations",
  };
  marginal.check("representation", settings.representation.is_some())?;
  marginal.check("ambiguity_threshold", settings.ambiguity_threshold.is_some())
}

#[derive(Debug, Clone, SmartDefault, JsonSchema, Serialize, Deserialize)]
#[deser(skip_serializing_optionals)]
#[schemars(transform = skip_serializing_optionals)]
#[schemars(default, deny_unknown_fields)]
#[schemars(transform = schema_defaults::<Self>)]
#[deser(default, deny_unknown_fields)]
#[cfg_attr(feature = "clap", derive(clap::Parser))]
#[schemars(rename = "AncestralConfig")]
pub struct TreetimeAncestralArgsRaw {
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
  pub alphabet_args: AlphabetArgs,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub model_args: ModelArgs,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub method_anc: MethodAncArgs,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub representation: RepresentationArgs,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub sample_from_profile: SampleFromProfileArgs,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub gap_fill_args: GapFillArgs,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub zero_based: ZeroBasedArgs,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub tip_sequences: TipSequenceArgs,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub impute_missing_data: ImputeMissingDataArgs,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub missing_sequences: MissingSequenceArgs,

  /// Path to output augur-compatible node data JSON.
  ///
  /// Contains per-node nucleotide mutations, reconstructed sequences, the alignment
  /// mask, genome annotations, and the reference (root) sequence. The output is
  /// compatible with augur export v2 --node-data for Nextstrain pipeline integration.
  ///
  /// Takes precedence over paths configured with `--output-all` and `--output-selection`.
  ///
  /// With several `--alignment` files, one file is written per alignment, and the file name
  /// must contain `{alignment}`, which is replaced by the number of the alignment (1, 2, ...).
  #[cfg_attr(feature = "clap", clap(long, value_hint = ValueHint::FilePath, help_heading = "Output"))]
  #[schemars(extend("x-path" = "output"))]
  pub output_node_data: Option<PathBuf>,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub gtr_output: GtrOutputArgs,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub reconstructed_fasta_output: ReconstructedFastaOutputArgs,

  /// Path template for per-CDS amino-acid FASTA alignments.
  ///
  /// The template must contain a CDS placeholder, replaced with each value from `--cds-selection` (or each
  /// CDS in `--annotation` when `--cds-selection` is omitted). Both `{cds}` (Nextclade
  /// `--output-translations`) and `%GENE` (augur) placeholders are accepted.
  #[cfg_attr(feature = "clap", clap(long, help_heading = "Input data"))]
  #[cfg_attr(feature = "clap", clap(value_hint = ValueHint::FilePath))]
  #[schemars(extend("x-path" = "input-template"))]
  pub translations: Option<String>,

  /// CDS names to reconstruct from `--translations`.
  ///
  /// When omitted, the CDS set is derived from `--annotation`.
  #[cfg_attr(
    feature = "clap",
    clap(
      long,
      value_name = "CDS",
      num_args = 1..,
      value_delimiter = ',',
      help_heading = "Input data"
    )
  )]
  pub cds_selection: Vec<String>,

  /// GFF3 file with CDS coordinates for Augur node data annotations.
  ///
  /// Also supplies the CDS set when `--cds-selection` is omitted.
  #[cfg_attr(feature = "clap", clap(long, help_heading = "Input data"))]
  #[cfg_attr(feature = "clap", clap(value_hint = ValueHint::FilePath))]
  #[schemars(extend("x-path" = "input"))]
  pub annotation: Option<PathBuf>,

  /// FASTA file with one amino-acid root/reference sequence per CDS.
  ///
  /// The root `aa_muts` of augur node data list the differences between this reference and the inferred root, and
  /// always omit `X` and gaps.
  #[cfg_attr(feature = "clap", clap(long, help_heading = "Ancestral reconstruction"))]
  #[cfg_attr(feature = "clap", clap(value_hint = ValueHint::FilePath))]
  #[schemars(extend("x-path" = "input"))]
  pub aa_root_sequence: Option<PathBuf>,

  /// Amino-acid substitution model. Mirrors the nucleotide `--model`; default `infer` matches augur.
  #[cfg_attr(feature = "clap", clap(long, value_enum, ignore_case = true, default_value_t = AaModelName::default(), help_heading = "Substitution model"))]
  pub aa_model: AaModelName,

  /// Path template for per-CDS reconstructed amino-acid FASTA output (including internal nodes).
  ///
  /// Off by default. When set, the reconstructed sequence of every node is written per CDS. The
  /// file name accepts the `{cds}`/`%GENE` placeholders of `--translations`, and needs one when
  /// more than one CDS is reconstructed; folder names are used as given.
  ///
  /// Takes precedence over paths configured with `--output-all` and `--output-selection`.
  #[cfg_attr(feature = "clap", clap(long, value_hint = ValueHint::FilePath, help_heading = "Output"))]
  #[schemars(extend("x-path" = "output"))]
  pub output_translations: Option<PathBuf>,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub output: OutputCoreArgs,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub output_selection: AncestralOutputSelectionArgs,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub ambiguity: AncestralAmbiguityArgs,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub topology_order: TopologyOrderArgs,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub seed_args: SeedArgs,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(skip)]
  #[deser(skip)]
  pub jobs_args: JobsArgs,
}
