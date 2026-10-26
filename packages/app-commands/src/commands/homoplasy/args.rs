use crate::commands::shared::alignment::{AlignmentArgs, check_required_alignment_paths};
use crate::commands::shared::alphabet::AlphabetArgs;
use crate::commands::shared::argument_error::argument_error;
use crate::commands::shared::config::ConfigArgs;
use crate::commands::shared::gap_fill::GapFillArgs;
use crate::commands::shared::jobs::JobsArgs;
use crate::commands::shared::method_anc::{
  MarginalReconstruction, MethodAncArgs, SampleFromProfileArgs, check_model_with_method_anc,
};
use crate::commands::shared::model::{ModelArgs, resolve_model};
use crate::commands::shared::output_args::{
  HomoplasyOutputSelectionArgs, OutputCoreArgs, OutputRequest, TreeOutputPaths, selected_outputs,
};
use crate::commands::shared::representation::{RepresentationArgs, representation_dense};
use crate::commands::shared::required::missing_required_args;
use crate::commands::shared::seed::SeedArgs;
use crate::commands::shared::sequence_inputs::{MissingSequenceArgs, SequenceInputs, missing_leaf_limit};
use crate::commands::shared::tip_sequences::ImputeMissingDataArgs;
use crate::commands::shared::topology_order_args::{TopologyOrder, TopologyOrderArgs};
use crate::commands::shared::tree_input::{TreeArgs, TreeDialectArgs};
#[cfg(feature = "clap")]
use crate::commands::shared::value_names::value_name_from_hint;
use crate::commands::shared::zero_based::ZeroBasedArgs;
use app_output::command_outputs::CommandKind;
use app_output::output_conditions::OutputSettings;
use app_output::output_plan::OutputSelection;
#[cfg(feature = "clap")]
use clap::ValueHint;
use deser::{Deserialize, Serialize};
use eyre::Report;
use schemars::JsonSchema;
use smart_default::SmartDefault;
use std::path::PathBuf;
use treetime::ancestral::params::{AncestralParams, MethodAncestral};
use treetime::gtr::get_gtr::GtrModel;
use treetime::partition::marginal::sample::SampleMode;
use treetime::seq::mutation::PositionBase;
use treetime_io::nwk::NewickDialect;
use treetime_schema::{schema_defaults, skip_serializing_optionals};

pub const HOMOPLASY_ABOUT: &str =
  "Reconstruct ancestral sequences and count homoplasies, mutations that occur more than once in a tree";

pub fn homoplasy_ancestral_params(args: &TreetimeHomoplasyArgs, seed: u64) -> AncestralParams {
  AncestralParams {
    method: args.method_anc,
    dense: args.dense,
    include_leaves: false,
    report_ambiguous: true,
    impute_missing_data: args.impute_missing_data,
    seed,
    sample_from_profile: args.sample_from_profile,
  }
}

#[derive(Debug, Clone)]
pub struct TreetimeHomoplasyArgs {
  pub(crate) tree: PathBuf,
  pub(crate) newick_dialect: NewickDialect,
  pub(crate) sequences: SequenceInputs,
  pub(crate) positions: PositionBase,
  pub(crate) method_anc: MethodAncestral,
  pub(crate) model: GtrModel,
  pub(crate) dense: Option<bool>,
  pub(crate) impute_missing_data: bool,
  pub(crate) sample_from_profile: SampleMode,
  pub(crate) seed: Option<u64>,
  pub(crate) constant_sites: usize,
  pub(crate) rescale: f64,
  pub(crate) detailed: bool,
  pub(crate) drms: Option<PathBuf>,
  pub(crate) num_mut: usize,
  pub(crate) topology_order: TopologyOrder,
  pub(crate) outputs: OutputRequest,
}

impl TreetimeHomoplasyArgs {
  pub fn alignments(&self) -> &[PathBuf] {
    &self.sequences.alignments
  }
}

impl TryFrom<TreetimeHomoplasyArgsRaw> for TreetimeHomoplasyArgs {
  type Error = Report;

  #[deny(unused_variables)]
  #[expect(
    clippy::unneeded_field_pattern,
    reason = "the CLI reads its own settings before the conversion, and `..` would hide every setting that no run reads"
  )]
  fn try_from(raw: TreetimeHomoplasyArgsRaw) -> Result<Self, Report> {
    let TreetimeHomoplasyArgsRaw {
      config_args: _,
      alignment: AlignmentArgs { alignment },
      tree: TreeArgs { tree },
      newick_dialect: TreeDialectArgs { newick_dialect },
      alphabet_args: AlphabetArgs { alphabet },
      model_args: ModelArgs { model, model_params },
      method_anc: MethodAncArgs { method_anc },
      representation: RepresentationArgs { representation },
      gap_fill_args: GapFillArgs { gap_fill },
      zero_based: ZeroBasedArgs { zero_based },
      impute_missing_data: ImputeMissingDataArgs { impute_missing_data },
      missing_sequences: MissingSequenceArgs {
        allow_missing_sequences,
      },
      sample_from_profile: SampleFromProfileArgs { sample_from_profile },
      seed_args: SeedArgs { seed },
      jobs_args: _,
      constant_sites,
      rescale,
      detailed,
      drms,
      num_mut,
      output,
      output_homoplasy_stats,
      output_homoplasy_report,
      output_selection,
      topology_order,
    } = raw;
    let tree = tree.ok_or_else(|| missing_required_args::<TreetimeHomoplasyArgsRaw>(&["tree"]))?;
    check_required_alignment_paths(&alignment)?;
    check_rescale(rescale)?;
    check_one_alignment(&alignment)?;
    check_model_with_method_anc(method_anc, model, &model_params)?;
    MarginalReconstruction {
      runs: method_anc.is_marginal(),
      condition: "--method-anc marginal",
    }
    .check("representation", representation.is_some())?;
    Ok(Self {
      tree,
      newick_dialect: newick_dialect.dialect(),
      sequences: SequenceInputs {
        alignments: alignment,
        alphabet: alphabet.map(Into::into),
        gap_fill: gap_fill.into(),
        missing_leaves: missing_leaf_limit(allow_missing_sequences),
        partition_entries: vec![],
      },
      positions: PositionBase::from_zero_based(zero_based),
      method_anc: method_anc.into(),
      model: resolve_model(model, &model_params)?,
      dense: representation_dense(representation),
      impute_missing_data,
      sample_from_profile: sample_from_profile.into(),
      seed,
      constant_sites,
      rescale,
      detailed,
      drms,
      num_mut,
      topology_order: homoplasy_topology_order(topology_order)?,
      outputs: homoplasy_outputs(
        output,
        output_selection,
        output_homoplasy_stats,
        output_homoplasy_report,
      ),
    })
  }
}

#[deny(unused_variables)]
fn homoplasy_outputs(
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
  HomoplasyOutputSelectionArgs { output_selection }: HomoplasyOutputSelectionArgs,
  output_homoplasy_stats: Option<PathBuf>,
  output_homoplasy_report: Option<PathBuf>,
) -> OutputRequest {
  OutputRequest {
    command: CommandKind::Homoplasy,
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
      (OutputSelection::HomoplasyStats, output_homoplasy_stats),
      (OutputSelection::HomoplasyReport, output_homoplasy_report),
    ],
  }
}

#[deny(unused_variables)]
fn homoplasy_topology_order(
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

fn check_rescale(rescale: f64) -> Result<(), Report> {
  if !rescale.is_finite() || rescale <= 0.0 {
    return Err(argument_error(
      &["rescale"],
      format!("--rescale must be a positive finite number, but got {rescale}"),
    ));
  }
  Ok(())
}

fn check_one_alignment(alignment: &[PathBuf]) -> Result<(), Report> {
  let n_alignments = alignment.len();
  if n_alignments > 1 {
    return Err(argument_error(
      &["alignment"],
      format!(
        "homoplasy reads exactly one --alignment file, but {n_alignments} were given: \
         the reported positions refer to the coordinates of one alignment"
      ),
    ));
  }
  Ok(())
}

#[derive(Debug, Clone, SmartDefault, JsonSchema, Serialize, Deserialize)]
#[deser(skip_serializing_optionals)]
#[schemars(transform = skip_serializing_optionals)]
#[schemars(default, deny_unknown_fields)]
#[schemars(transform = schema_defaults::<Self>)]
#[deser(default, deny_unknown_fields)]
#[cfg_attr(feature = "clap", derive(clap::Parser))]
#[cfg_attr(feature = "clap", command(mut_args(value_name_from_hint)))]
#[schemars(rename = "HomoplasyConfig", description = HOMOPLASY_ABOUT)]
pub struct TreetimeHomoplasyArgsRaw {
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
  pub missing_sequences: MissingSequenceArgs,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub alphabet_args: AlphabetArgs,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub zero_based: ZeroBasedArgs,

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
  pub impute_missing_data: ImputeMissingDataArgs,

  /// Number of constant sites that the alignment leaves out.
  ///
  /// Added to the alignment length to give the number of sites of the genome, which sets the
  /// number of sites without mutations and the rate of the Poisson comparison.
  #[cfg_attr(
    feature = "clap",
    clap(long = "const", value_name = "N", default_value_t = 0, help_heading = "Homoplasy")
  )]
  pub constant_sites: usize,

  /// Factor that multiplies every branch length of the input tree.
  ///
  /// Use it when the tree is not in substitutions per site, for example `--rescale=0.001` for a
  /// tree in substitutions per thousand sites. Scales the reconstruction, the tree outputs, and the
  /// reported tree lengths.
  #[cfg_attr(feature = "clap", clap(long, default_value_t = 1.0, help_heading = "Homoplasy"))]
  #[default = 1.0]
  pub rescale: f64,

  /// Add the mutations on terminal branches and the taxa that carry recurrent mutations to the
  /// report.
  ///
  /// The statistics JSON always contains them.
  #[cfg_attr(feature = "clap", clap(long, help_heading = "Homoplasy"))]
  pub detailed: bool,

  /// TSV file of drug resistance mutations (DRM) that annotates the report.
  ///
  /// Columns: `GENOMIC_POSITION` (counted from 1), `ALT_BASE`, `DRUG`, `GENE`, `SUBSTITUTION`. A mutation at a listed
  /// position gets the gene and the drug, and the substitution when its derived base is a listed `ALT_BASE`. A path
  /// that ends in `.gz`, `.bz2`, `.xz`, or `.zst` reads a compressed file.
  #[cfg_attr(feature = "clap", clap(long, help_heading = "Homoplasy"))]
  #[cfg_attr(feature = "clap", clap(value_hint = ValueHint::FilePath))]
  #[schemars(extend("x-path" = "input"))]
  pub drms: Option<PathBuf>,

  /// Number of rows of each list in the report.
  ///
  /// The statistics JSON always contains the complete lists.
  #[cfg_attr(
    feature = "clap",
    clap(
      long,
      value_name = "N",
      short = 'n',
      default_value_t = 10,
      help_heading = "Homoplasy"
    )
  )]
  #[default = 10]
  pub num_mut: usize,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub output: OutputCoreArgs,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub output_selection: HomoplasyOutputSelectionArgs,

  /// Path to output homoplasy statistics JSON.
  ///
  /// Contains the counts, histograms, Poisson comparison, and complete ranked lists of the report, for substitutions,
  /// changes involving ambiguous characters, and insertions and deletions.
  #[cfg_attr(feature = "clap", clap(long, hide_short_help = true, value_hint = ValueHint::FilePath, help_heading = "Output"))]
  #[schemars(extend("x-path" = "output"))]
  pub output_homoplasy_stats: Option<PathBuf>,

  /// Path to output homoplasy report text.
  ///
  /// The run also logs the report at info level (`-v`).
  #[cfg_attr(feature = "clap", clap(long, hide_short_help = true, value_hint = ValueHint::FilePath, help_heading = "Output"))]
  #[schemars(extend("x-path" = "output"))]
  pub output_homoplasy_report: Option<PathBuf>,

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
