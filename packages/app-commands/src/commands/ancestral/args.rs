use crate::commands::ancestral::aa_model::AaModelName;
use crate::commands::shared::alignment::AlignmentArgs;
use crate::commands::shared::alphabet::AlphabetArgs;
use crate::commands::shared::ambiguity::AncestralAmbiguityArgs;
use crate::commands::shared::config::ConfigArgs;
use crate::commands::shared::gap_fill::GapFillArgs;
use crate::commands::shared::gtr_output::GtrOutputArgs;
use crate::commands::shared::jobs::JobsArgs;
use crate::commands::shared::method_anc::{MethodAncArgs, SampleFromProfileArgs};
use crate::commands::shared::model::ModelArgs;
use crate::commands::shared::output_args::{AncestralOutputSelectionArgs, OutputCoreArgs};
use crate::commands::shared::reconstructed_fasta::ReconstructedFastaOutputArgs;
use crate::commands::shared::representation::RepresentationArgs;
use crate::commands::shared::required::missing_required_args;
use crate::commands::shared::seed::SeedArgs;
use crate::commands::shared::sequence_inputs::MissingSequenceArgs;
use crate::commands::shared::tip_sequences::{ImputeMissingDataArgs, TipSequenceArgs};
use crate::commands::shared::topology_order_args::TopologyOrderArgs;
use crate::commands::shared::tree_input::{TreeArgs, TreeDialectArgs};
use crate::commands::shared::zero_based::ZeroBasedArgs;
use app_output::ambiguity_outputs::AmbiguityOutputs;
#[cfg(feature = "clap")]
use clap::ValueHint;
use derive_more::{Deref, DerefMut};
use deser::{Deserialize, Serialize};
use eyre::Report;
use schemars::JsonSchema;
use smart_default::SmartDefault;
use std::path::{Path, PathBuf};
use treetime::ancestral::params::AncestralParams;
use treetime::ancestral::params::MethodAncestral;
use treetime::gtr::get_gtr::GtrModel;
use treetime::partition::marginal::sample::SampleMode;
use treetime_schema::{schema_defaults, skip_serializing_optionals};

pub fn ancestral_params(args: &TreetimeAncestralArgs, seed: u64) -> AncestralParams {
  AncestralParams {
    method: args.method_anc,
    model: args.model.clone(),
    dense: args.representation.dense(),
    include_leaves: args.include_leaves,
    report_ambiguous: args.tip_sequences.report_ambiguous,
    impute_missing_data: args.impute_missing_data,
    seed,
    sample_from_profile: args.sample_from_profile,
  }
}

#[derive(Debug, Clone, Deref, DerefMut)]
pub struct TreetimeAncestralArgs {
  #[deref]
  #[deref_mut]
  pub(crate) raw: TreetimeAncestralArgsRaw,
  pub(crate) tree: PathBuf,
  pub(crate) method_anc: MethodAncestral,
  pub(crate) include_leaves: bool,
  pub(crate) impute_missing_data: bool,
  pub(crate) sample_from_profile: SampleMode,
  pub(crate) ambiguity_outputs: AmbiguityOutputs,
  pub(crate) model: GtrModel,
}

impl TreetimeAncestralArgs {
  pub fn tree(&self) -> &Path {
    &self.tree
  }
}

impl TryFrom<TreetimeAncestralArgsRaw> for TreetimeAncestralArgs {
  type Error = Report;

  fn try_from(raw: TreetimeAncestralArgsRaw) -> Result<Self, Report> {
    let tree = raw
      .tree
      .tree
      .clone()
      .ok_or_else(|| missing_required_args::<TreetimeAncestralArgsRaw>(&["tree"]))?;
    Ok(Self {
      tree,
      method_anc: raw.method_anc.method_anc.into(),
      include_leaves: raw.tip_sequences.include_leaves || raw.tip_sequences.reconstruct_tip_states,
      impute_missing_data: raw.impute_missing_data.impute_missing_data || raw.tip_sequences.reconstruct_tip_states,
      sample_from_profile: raw.sample_from_profile.sample_from_profile.into(),
      ambiguity_outputs: raw.ambiguity.resolve()?,
      model: raw.model_args.model()?,
      raw,
    })
  }
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
  pub tree_dialect: TreeDialectArgs,

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
  /// must contain `{partition}`, which is replaced by the number of the alignment (1, 2, ...).
  #[cfg_attr(feature = "clap", clap(long, value_hint = ValueHint::FilePath, help_heading = "Output"))]
  #[schemars(extend("x-path" = "output"))]
  pub output_augur_node_data: Option<PathBuf>,

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
  /// The template must contain a CDS placeholder, replaced with each value from `--cdses` (or each
  /// CDS in `--annotation` when `--cdses` is omitted). Both `{cds}` (Nextclade
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
      long = "cdses",
      visible_alias = "genes",
      value_name = "CDS",
      num_args = 1..,
      value_delimiter = ',',
      help_heading = "Input data"
    )
  )]
  pub cdses: Vec<String>,

  /// GFF3 file with CDS coordinates for Augur node data annotations.
  ///
  /// Also supplies the CDS set when `--cdses` is omitted.
  #[cfg_attr(feature = "clap", clap(long, alias = "annotation-gff", help_heading = "Input data"))]
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
  #[cfg_attr(feature = "clap", clap(long, value_enum, default_value_t = AaModelName::default(), help_heading = "Substitution model"))]
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
  pub output_reconstructed_aa_fasta: Option<PathBuf>,

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
