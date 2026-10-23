use crate::commands::shared::alignment::AlignmentArgs;
use crate::commands::shared::alphabet::AlphabetArgs;
use crate::commands::shared::ambiguity::{TimetreeAmbiguityArgs, resolve_ambiguity_outputs};
use crate::commands::shared::argument_error::argument_error;
use crate::commands::shared::branch_lengths::{BranchLengthArgs, indel_policy};
use crate::commands::shared::clock_filter::ClockFilterArgs;
use crate::commands::shared::config::ConfigArgs;
use crate::commands::shared::dates_input::DatesInput;
use crate::commands::shared::dependent_settings::{check_applies, setting_flag};
use crate::commands::shared::documented_default::documented_default;
use crate::commands::shared::gap_fill::GapFillArgs;
use crate::commands::shared::gtr_output::GtrOutputArgs;
use crate::commands::shared::jobs::JobsArgs;
use crate::commands::shared::metadata::{DateColumnArgs, MetadataArgs, MetadataIdArgs, MetadataInput};
use crate::commands::shared::model::{
  ModelArgs, PartitionRatesArgs, check_partition_rates, resolve_model, resolve_partition_rates,
};
use crate::commands::shared::output_args::{
  ClockOutputArgs, DivergenceUnits, OutputCoreArgs, OutputRequest, TimetreeOutputSelectionArgs, TreeOutputPaths,
  selected_outputs,
};
use crate::commands::shared::partition::{PartitionArgs, PartitionSpec, partition_entries};
use crate::commands::shared::reconstructed_fasta::ReconstructedFastaOutputArgs;
use crate::commands::shared::representation::{RepresentationArgs, representation_dense};
use crate::commands::shared::required::missing_required_args;
use crate::commands::shared::reroot::{RerootArgs, RerootTarget, RootingArgs};
use crate::commands::shared::seed::SeedArgs;
use crate::commands::shared::sequence_inputs::{MissingSequenceArgs, SequenceInputs, missing_leaf_limit};
use crate::commands::shared::tip_sequences::{ImputeMissingDataArgs, TipSequenceArgs};
use crate::commands::shared::topology_order_args::{TopologyOrder, TopologyOrderArgs};
use crate::commands::shared::tree_input::{TreeArgs, TreeDialectArgs};
#[cfg(feature = "clap")]
use crate::commands::shared::value_names::value_name_from_hint;
use crate::commands::shared::zero_based::ZeroBasedArgs;
use crate::commands::timetree::branch_length_mode::{BranchLengthModeCli, ResolvedBranchLengthMode};
use crate::commands::timetree::coalescent::CoalescentArg;
#[cfg(feature = "clap")]
use crate::commands::timetree::coalescent::parse_coalescent;
use app_output::ambiguity_outputs::AmbiguityOutputs;
use app_output::command_outputs::CommandKind;
use app_output::output_conditions::{NotReconstructed, OutputSettings, SequenceOutputs};
use app_output::output_plan::OutputSelection;
#[cfg(feature = "clap")]
use clap::ValueHint;
use deser::{Deserialize, Serialize};
use eyre::Report;
use itertools::Itertools;
use schemars::JsonSchema;
use smart_default::SmartDefault;
use std::path::PathBuf;
use treetime::gtr::get_gtr::GtrModel;
use treetime::optimize::indel::IndelPolicy;
use treetime::optimize::partition_rates::PartitionRates;
use treetime::seq::mutation::PositionBase;
use treetime::timetree::coalescent_timescale::CoalescentMode;
use treetime_grid::MaxGridPoints;
use treetime_io::nwk::NewickDialect;
use treetime_schema::{schema_defaults, skip_serializing_optionals};
use treetime_utils::o;

pub const TIMETREE_ABOUT: &str =
  "Infer a time-scaled tree from a tree, the sampling dates of its tips, and optionally an alignment";

pub const N_SKYLINE_DEFAULT: usize = 20;

const N_SKYLINE_MIN: usize = 2;

pub const SKYLINE_STIFFNESS_DEFAULT: f64 = 2.0;

pub const COALESCENT_CONFIDENCE_DEFAULT: f64 = 2.0;

pub const GEN_PER_YEAR_DEFAULT: f64 = 50.0;

#[derive(Debug, Clone)]
pub struct TreetimeTimetreeArgs {
  pub(crate) tree: PathBuf,
  pub(crate) newick_dialect: NewickDialect,
  pub(crate) dates: DatesInput,
  pub(crate) sequences: SequenceInputs,
  pub(crate) branch_lengths: TimetreeBranchLengths,
  pub(crate) sequence_model: TimetreeSequenceModel,
  pub(crate) marginal_only_settings: Vec<&'static str>,
  pub(crate) positions: PositionBase,
  pub(crate) clock: TimetreeClock,
  pub(crate) reroot: RerootTarget,
  pub(crate) keep_root: bool,
  pub(crate) dating: TimetreeDating,
  pub(crate) tip_sequences: TimetreeTipSequences,
  pub(crate) seed: Option<u64>,
  pub(crate) ambiguity_outputs: AmbiguityOutputs,
  pub(crate) topology_order: TopologyOrder,
  pub(crate) outputs: OutputRequest,
}

impl TreetimeTimetreeArgs {
  pub fn alignments(&self) -> &[PathBuf] {
    &self.sequences.alignments
  }
}

impl TryFrom<TreetimeTimetreeArgsRaw> for TreetimeTimetreeArgs {
  type Error = Report;

  #[deny(unused_variables)]
  #[expect(
    clippy::unneeded_field_pattern,
    reason = "the CLI reads its own settings before the conversion, and `..` would hide every setting that no run reads"
  )]
  fn try_from(raw: TreetimeTimetreeArgsRaw) -> Result<Self, Report> {
    let TreetimeTimetreeArgsRaw {
      config_args: _,
      alignment: AlignmentArgs { alignment },
      partition: PartitionArgs { partition },
      tree: TreeArgs { tree },
      newick_dialect: TreeDialectArgs { newick_dialect },
      missing_sequences,
      metadata,
      metadata_id,
      date_column,
      sequence_length,
      clock_rate,
      clock_std_dev,
      branch_length_mode,
      confidence,
      max_grid_points,
      keep_polytomies,
      relax,
      max_iter,
      coalescent,
      n_skyline,
      skyline_stiffness,
      coalescent_confidence,
      clock_filter,
      reroot: RerootArgs { reroot },
      rooting: RootingArgs { outgroup, keep_root },
      model_args,
      partition_rates_args,
      alphabet_args,
      representation,
      gap_fill_args,
      zero_based: ZeroBasedArgs { zero_based },
      tip_sequences,
      impute_missing_data,
      branch_lengths,
      output_node_data,
      gtr_output: GtrOutputArgs { output_model },
      reconstructed_fasta_output: ReconstructedFastaOutputArgs { output_sequences },
      clock_output: ClockOutputArgs {
        output_clock_model,
        output_rtt,
      },
      output_dates,
      output_tracelog,
      output_coalescent,
      output_coalescent_json,
      output,
      output_selection,
      ambiguity,
      topology_order,
      seed_args: SeedArgs { seed },
      jobs_args: _,
      gen_per_year,
    } = raw;
    let (tree, dates) = timetree_inputs(tree, metadata, metadata_id, date_column)?;
    let alignment_count = alignment.len();
    let marginal_only_settings = timetree_marginal_only_settings(&representation, &partition_rates_args, &ambiguity);
    check_marginal_only_settings(&marginal_only_settings, branch_length_mode, alignment_count > 0)?;
    let output_settings = timetree_output_settings(branch_length_mode, alignment_count > 0, coalescent.is_some());
    let coalescent = coalescent_settings(
      coalescent,
      n_skyline,
      skyline_stiffness,
      coalescent_confidence,
      gen_per_year,
    )?;
    let sequences = timetree_sequence_inputs(alignment, &partition, alphabet_args, gap_fill_args, missing_sequences)?;
    let file_paths = vec![
      (OutputSelection::NodeData, output_node_data),
      (OutputSelection::Model, output_model),
      (OutputSelection::Sequences, output_sequences),
      (OutputSelection::ClockModel, output_clock_model),
      (OutputSelection::Rtt, output_rtt),
      (OutputSelection::Dates, output_dates),
      (OutputSelection::Tracelog, output_tracelog),
      (OutputSelection::Coalescent, output_coalescent),
      (OutputSelection::CoalescentJson, output_coalescent_json),
    ];
    Ok(Self {
      tree,
      newick_dialect: newick_dialect.dialect(),
      dates,
      branch_lengths: timetree_branch_lengths(branch_length_mode, alignment_count, sequence_length, branch_lengths)?,
      sequence_model: timetree_sequence_model(model_args, representation, partition_rates_args, &sequences)?,
      marginal_only_settings,
      clock: timetree_clock(clock_rate, clock_std_dev, confidence, clock_filter, relax)?,
      reroot: RerootTarget::new(reroot.map(Into::into), outgroup),
      keep_root,
      dating: timetree_dating(max_iter, keep_polytomies, confidence, max_grid_points, coalescent),
      outputs: timetree_outputs(output, output_selection, alignment_count, output_settings, file_paths),
      sequences,
      positions: PositionBase::from_zero_based(zero_based),
      tip_sequences: timetree_tip_sequences(tip_sequences, impute_missing_data),
      seed,
      ambiguity_outputs: timetree_ambiguity_outputs(ambiguity)?,
      topology_order: timetree_topology_order(topology_order)?,
    })
  }
}

#[deny(unused_variables)]
fn timetree_inputs(
  tree: Option<PathBuf>,
  MetadataArgs { metadata }: MetadataArgs,
  MetadataIdArgs {
    metadata_id_columns,
    metadata_delimiters,
  }: MetadataIdArgs,
  DateColumnArgs {
    date_column,
    date_format,
  }: DateColumnArgs,
) -> Result<(PathBuf, DatesInput), Report> {
  let (tree, metadata) = match (tree, metadata) {
    (Some(tree), Some(metadata)) => (tree, metadata),
    (tree, metadata) => {
      let missing = [("tree", tree.is_none()), ("metadata", metadata.is_none())]
        .into_iter()
        .filter_map(|(key, missing)| missing.then_some(key))
        .collect_vec();
      return Err(missing_required_args::<TreetimeTimetreeArgsRaw>(&missing));
    },
  };
  let dates = DatesInput {
    metadata: MetadataInput {
      path: metadata,
      id_columns: metadata_id_columns,
      delimiters: metadata_delimiters,
    },
    date_column,
    date_format,
  };
  Ok((tree, dates))
}

fn timetree_marginal_only_settings(
  representation: &RepresentationArgs,
  partition_rates: &PartitionRatesArgs,
  ambiguity: &TimetreeAmbiguityArgs,
) -> Vec<&'static str> {
  [
    ("representation", representation.representation.is_some()),
    ("partition_rates", partition_rates.partition_rates.is_some()),
    ("ambiguity_threshold", ambiguity.ambiguity_threshold.is_some()),
  ]
  .into_iter()
  .filter_map(|(setting, given)| given.then_some(setting))
  .collect_vec()
}

fn check_marginal_only_settings(
  marginal_only_settings: &[&str],
  branch_length_mode: BranchLengthModeCli,
  has_alignment: bool,
) -> Result<(), Report> {
  if branch_length_mode == BranchLengthModeCli::Auto {
    return Ok(());
  }
  check_branch_length_settings(
    marginal_only_settings,
    &ResolvedBranchLengthMode::new(branch_length_mode, has_alignment, None),
  )
}

#[deny(unused_variables)]
fn timetree_branch_lengths(
  branch_length_mode: BranchLengthModeCli,
  alignment_count: usize,
  sequence_length: Option<usize>,
  BranchLengthArgs {
    no_indels,
    divergence_units,
  }: BranchLengthArgs,
) -> Result<TimetreeBranchLengths, Report> {
  check_sequence_source(alignment_count > 0, sequence_length, branch_length_mode)?;
  Ok(TimetreeBranchLengths {
    mode: branch_length_mode,
    sequence_length,
    indel_policy: indel_policy(no_indels),
    mutation_units: divergence_units == DivergenceUnits::Mutations,
  })
}

#[derive(Debug, Clone, Copy)]
pub(crate) struct TimetreeBranchLengths {
  pub(crate) mode: BranchLengthModeCli,
  pub(crate) sequence_length: Option<usize>,
  pub(crate) indel_policy: IndelPolicy,
  pub(crate) mutation_units: bool,
}

#[deny(unused_variables)]
fn timetree_sequence_model(
  ModelArgs { model, model_params }: ModelArgs,
  RepresentationArgs { representation }: RepresentationArgs,
  PartitionRatesArgs { partition_rates }: PartitionRatesArgs,
  sequences: &SequenceInputs,
) -> Result<TimetreeSequenceModel, Report> {
  check_partition_rates(
    partition_rates,
    sequences.alignments.len(),
    sequences.partition_entries.len(),
  )?;
  Ok(TimetreeSequenceModel {
    model: resolve_model(model, &model_params)?,
    dense: representation_dense(representation),
    partition_rates: resolve_partition_rates(partition_rates),
  })
}

#[derive(Debug, Clone)]
pub(crate) struct TimetreeSequenceModel {
  pub(crate) model: GtrModel,
  pub(crate) dense: Option<bool>,
  pub(crate) partition_rates: PartitionRates,
}

#[deny(unused_variables)]
fn timetree_clock(
  clock_rate: Option<f64>,
  clock_std_dev: Option<f64>,
  confidence: bool,
  ClockFilterArgs {
    clock_filter,
    allow_negative_rate,
  }: ClockFilterArgs,
  relax: Vec<f64>,
) -> Result<TimetreeClock, Report> {
  check_applies(
    "clock_std_dev",
    clock_std_dev.is_some(),
    "confidence",
    confidence,
    || o!("applies only with --confidence: only the date intervals use the clock rate uncertainty"),
  )?;
  Ok(TimetreeClock {
    clock_rate,
    clock_std_dev,
    clock_filter,
    allow_negative_rate,
    relax,
  })
}

#[derive(Debug, Clone)]
pub(crate) struct TimetreeClock {
  pub(crate) clock_rate: Option<f64>,
  pub(crate) clock_std_dev: Option<f64>,
  pub(crate) clock_filter: f64,
  pub(crate) allow_negative_rate: bool,
  pub(crate) relax: Vec<f64>,
}

#[deny(unused_variables)]
fn timetree_dating(
  max_iter: usize,
  keep_polytomies: bool,
  confidence: bool,
  max_grid_points: Option<MaxGridPoints>,
  coalescent: CoalescentSettings,
) -> TimetreeDating {
  TimetreeDating {
    max_iter,
    keep_polytomies,
    confidence,
    max_grid_points: max_grid_points.unwrap_or_default(),
    coalescent,
  }
}

#[derive(Debug, Clone, Copy)]
pub(crate) struct TimetreeDating {
  pub(crate) max_iter: usize,
  pub(crate) keep_polytomies: bool,
  pub(crate) confidence: bool,
  pub(crate) max_grid_points: MaxGridPoints,
  pub(crate) coalescent: CoalescentSettings,
}

#[deny(unused_variables)]
fn timetree_outputs(
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
  TimetreeOutputSelectionArgs { output_selection }: TimetreeOutputSelectionArgs,
  alignment_count: usize,
  settings: OutputSettings,
  file_paths: Vec<(OutputSelection, Option<PathBuf>)>,
) -> OutputRequest {
  OutputRequest {
    command: CommandKind::Timetree,
    alignment_count,
    settings,
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
    file_paths,
  }
}

#[deny(unused_variables)]
fn timetree_sequence_inputs(
  alignment: Vec<PathBuf>,
  partition: &[PartitionSpec],
  AlphabetArgs { alphabet }: AlphabetArgs,
  GapFillArgs { gap_fill }: GapFillArgs,
  MissingSequenceArgs {
    allow_missing_sequences,
  }: MissingSequenceArgs,
) -> Result<SequenceInputs, Report> {
  Ok(SequenceInputs {
    alignments: alignment,
    alphabet: alphabet.map(Into::into),
    gap_fill: gap_fill.into(),
    missing_leaves: missing_leaf_limit(allow_missing_sequences),
    partition_entries: partition_entries(partition)?,
  })
}

#[deny(unused_variables)]
fn timetree_tip_sequences(
  TipSequenceArgs {
    include_leaves,
    reconstruct_tip_states,
    report_ambiguous,
  }: TipSequenceArgs,
  ImputeMissingDataArgs { impute_missing_data }: ImputeMissingDataArgs,
) -> TimetreeTipSequences {
  TimetreeTipSequences {
    include_leaves: include_leaves || reconstruct_tip_states,
    report_ambiguous,
    impute_missing_data: impute_missing_data || reconstruct_tip_states,
  }
}

#[derive(Debug, Clone, Copy)]
pub(crate) struct TimetreeTipSequences {
  pub(crate) include_leaves: bool,
  pub(crate) report_ambiguous: bool,
  pub(crate) impute_missing_data: bool,
}

#[deny(unused_variables)]
fn timetree_ambiguity_outputs(
  TimetreeAmbiguityArgs {
    ambiguity_threshold,
    ambiguity_outputs,
  }: TimetreeAmbiguityArgs,
) -> Result<AmbiguityOutputs, Report> {
  resolve_ambiguity_outputs(ambiguity_threshold, &ambiguity_outputs)
}

#[deny(unused_variables)]
fn timetree_topology_order(
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

pub(crate) fn check_branch_length_settings(
  marginal_only_settings: &[&str],
  resolved: &ResolvedBranchLengthMode,
) -> Result<(), Report> {
  for setting in marginal_only_settings {
    resolved.require_marginal(&[setting], &setting_flag(setting))?;
  }
  Ok(())
}

pub(crate) fn timetree_output_settings(
  branch_length_mode: BranchLengthModeCli,
  has_alignment: bool,
  coalescent: bool,
) -> OutputSettings {
  let sequences = match (branch_length_mode, has_alignment) {
    (BranchLengthModeCli::Input, _) => SequenceOutputs::NotReconstructed(NotReconstructed::InputBranchLengths),
    (_, false) => SequenceOutputs::NotReconstructed(NotReconstructed::NoAlignment),
    (BranchLengthModeCli::Marginal, true) => SequenceOutputs::Reconstructed,
    (BranchLengthModeCli::Auto, true) => SequenceOutputs::WhenBranchLengthsReestimated,
  };
  OutputSettings {
    coalescent,
    trait_probabilities: false,
    sequences,
  }
}

#[derive(Debug, Clone, Copy)]
pub(crate) struct CoalescentSettings {
  pub(crate) mode: CoalescentMode,
  pub(crate) n_skyline: usize,
  pub(crate) skyline_stiffness: f64,
  pub(crate) confidence: f64,
  pub(crate) gen_per_year: f64,
}

fn check_sequence_source(
  has_alignment: bool,
  sequence_length: Option<usize>,
  branch_length_mode: BranchLengthModeCli,
) -> Result<(), Report> {
  if has_alignment {
    return check_applies("sequence_length", sequence_length.is_some(), "alignment", false, || {
      o!(
        "applies only without --alignment: the run takes the sequence length from the alignment, and an \
            alignment of variable sites with a separate genome length is not supported"
      )
    });
  }
  if branch_length_mode == BranchLengthModeCli::Marginal {
    return Err(argument_error(
      &["alignment", "branch_length_mode"],
      "--branch-length-mode=marginal requires an alignment: add --alignment or use --branch-length-mode=input",
    ));
  }
  if sequence_length.is_none_or(|length| length == 0) {
    return Err(argument_error(
      &["alignment", "sequence_length"],
      "--sequence-length is required without an alignment: the time inference needs the sequence length",
    ));
  }
  Ok(())
}

fn coalescent_settings(
  coalescent: Option<CoalescentArg>,
  n_skyline: Option<usize>,
  skyline_stiffness: Option<f64>,
  coalescent_confidence: Option<f64>,
  gen_per_year: Option<f64>,
) -> Result<CoalescentSettings, Report> {
  let actual = || coalescent.map_or_else(|| o!("not set"), |coalescent| coalescent.to_string());
  let skyline = matches!(coalescent, Some(CoalescentArg::Skyline));
  for (setting, given) in [
    ("n_skyline", n_skyline.is_some()),
    ("skyline_stiffness", skyline_stiffness.is_some()),
  ] {
    check_applies(setting, given, "coalescent", skyline, || {
      format!("applies only to --coalescent skyline, but --coalescent is {}", actual())
    })?;
  }
  if let Some(n_skyline) = n_skyline.filter(|&n_skyline| n_skyline < N_SKYLINE_MIN) {
    return Err(argument_error(
      &["n_skyline"],
      format!("--n-skyline is {n_skyline}, but the skyline needs at least {N_SKYLINE_MIN} grid points"),
    ));
  }
  let fitted = matches!(coalescent, Some(CoalescentArg::Opt | CoalescentArg::Skyline));
  check_applies(
    "coalescent_confidence",
    coalescent_confidence.is_some(),
    "coalescent",
    fitted,
    || {
      format!(
        "applies only to --coalescent opt and --coalescent skyline, but --coalescent is {}",
        actual()
      )
    },
  )?;
  check_applies(
    "gen_per_year",
    gen_per_year.is_some(),
    "coalescent",
    coalescent.is_some(),
    || o!("applies only with --coalescent, which is not set"),
  )?;
  Ok(CoalescentSettings {
    mode: coalescent.map_or(CoalescentMode::Disabled, CoalescentArg::mode),
    n_skyline: n_skyline.unwrap_or(N_SKYLINE_DEFAULT),
    skyline_stiffness: skyline_stiffness.unwrap_or(SKYLINE_STIFFNESS_DEFAULT),
    confidence: coalescent_confidence.unwrap_or(COALESCENT_CONFIDENCE_DEFAULT),
    gen_per_year: gen_per_year.unwrap_or(GEN_PER_YEAR_DEFAULT),
  })
}

impl TreetimeTimetreeArgsRaw {
  pub fn n_skyline(&self) -> usize {
    self.n_skyline.unwrap_or(N_SKYLINE_DEFAULT)
  }

  pub fn skyline_stiffness(&self) -> f64 {
    self.skyline_stiffness.unwrap_or(SKYLINE_STIFFNESS_DEFAULT)
  }

  pub fn coalescent_confidence(&self) -> f64 {
    self.coalescent_confidence.unwrap_or(COALESCENT_CONFIDENCE_DEFAULT)
  }

  pub fn gen_per_year(&self) -> f64 {
    self.gen_per_year.unwrap_or(GEN_PER_YEAR_DEFAULT)
  }
}

#[derive(Debug, Clone, SmartDefault, JsonSchema, Serialize, Deserialize)]
#[deser(skip_serializing_optionals)]
#[schemars(transform = skip_serializing_optionals)]
#[schemars(default, deny_unknown_fields)]
#[schemars(transform = schema_defaults::<Self>)]
#[deser(default, deny_unknown_fields)]
#[cfg_attr(feature = "clap", derive(clap::Parser))]
#[cfg_attr(
  feature = "clap",
  command(
    mut_arg("n_skyline", |arg| documented_default(&N_SKYLINE_DEFAULT).help(arg)),
    mut_arg("skyline_stiffness", |arg| documented_default(&SKYLINE_STIFFNESS_DEFAULT).help(arg)),
    mut_arg("coalescent_confidence", |arg| documented_default(&COALESCENT_CONFIDENCE_DEFAULT).help(arg)),
    mut_arg("gen_per_year", |arg| documented_default(&GEN_PER_YEAR_DEFAULT).help(arg)),
    mut_args(value_name_from_hint),
  )
)]
#[schemars(rename = "TimetreeConfig", description = TIMETREE_ABOUT)]
pub struct TreetimeTimetreeArgsRaw {
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
  pub partition: PartitionArgs,

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
  pub metadata: MetadataArgs,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub metadata_id: MetadataIdArgs,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub date_column: DateColumnArgs,

  /// Length of the sequences, used to calculate the expected variation of the branch lengths
  ///
  /// Required without `--alignment`, and applies only then: with an alignment, the run takes the length from the
  /// alignment.
  #[cfg_attr(feature = "clap", clap(long, value_name = "N", help_heading = "Input data"))]
  pub sequence_length: Option<usize>,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub alphabet_args: AlphabetArgs,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub zero_based: ZeroBasedArgs,

  /// Fixed rate of the molecular clock, in substitutions per site per year
  ///
  /// Without this setting, the run estimates the rate.
  #[schemars(example = 0.001)]
  #[cfg_attr(feature = "clap", clap(long, help_heading = "Molecular clock"))]
  pub clock_rate: Option<f64>,

  /// Standard deviation of the provided clock rate estimate
  ///
  /// Applies only with `--confidence`, whose date intervals include the clock rate uncertainty.
  #[cfg_attr(feature = "clap", clap(long, help_heading = "Molecular clock"))]
  pub clock_std_dev: Option<f64>,

  /// Use an autocorrelated relaxed clock with the given prior strengths, for example `--relax 1.0 0.5`
  ///
  /// SLACK is the strength of the Gaussian prior on the rate deviation of each branch, and COUPLING the strength of
  /// the coupling of the rates of parent and child branches. Values around 1.0 give weak priors; larger values
  /// constrain the rate deviations more strongly. A COUPLING of 0 (`--relax 1.0 0`) gives an uncorrelated clock.
  #[schemars(example = [1.0, 0.0])]
  #[cfg_attr(feature = "clap", clap(long, num_args = 2, value_names = ["SLACK", "COUPLING"], help_heading = "Molecular clock"))]
  pub relax: Vec<f64>,

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

  /// How the branch lengths of the time inference are obtained
  ///
  /// `auto` keeps the input branch lengths without an alignment, and also when the longest input branch is longer
  /// than 0.1 substitutions per site, because the re-estimated branch lengths are accurate only at short evolutionary
  /// distances; otherwise it re-estimates them as `marginal` does. When no input branch has a length, `auto` chooses
  /// `marginal`. The run logs the choice.
  #[cfg_attr(feature = "clap", clap(long, value_enum, default_value_t = BranchLengthModeCli::default(), help_heading = "Branch lengths"))]
  pub branch_length_mode: BranchLengthModeCli,

  /// Maximum number of iterations of the inference loop
  ///
  /// Polytomy resolution and coalescent models of the time inference need at least 2 iterations.
  #[default = 2]
  #[cfg_attr(feature = "clap", clap(long, value_name = "N", default_value_t = TreetimeTimetreeArgsRaw::default().max_iter, help_heading = "Branch lengths"))]
  pub max_iter: usize,

  /// Estimate the 90% intervals of the node dates
  ///
  /// Runs a final marginal round of the time inference after the iterations and writes the 90% interval of every
  /// dated node. The intervals include the clock rate uncertainty when one is known: `--clock-std-dev` for a given
  /// rate, or the Gaussian estimate of an estimated rate. The rate uncertainty is added by inferring the dates again
  /// at the rate plus and minus one standard deviation and combining both parts in quadrature.
  #[cfg_attr(feature = "clap", clap(long, help_heading = "Dating"))]
  pub confidence: bool,

  /// Largest number of points of one probability grid during time inference
  ///
  /// A run stops with an error when a grid would need more points. The value bounds single grids, not the total
  /// memory of a run. When unset, the run takes the limit of the app or server that runs it, otherwise 1000000.
  #[schemars(example = 1000000)]
  #[cfg_attr(
    feature = "clap",
    clap(long, value_name = "N", hide_short_help = true, help_heading = "Dating")
  )]
  pub max_grid_points: Option<MaxGridPoints>,

  /// Keep the polytomies of the input tree
  ///
  /// Without this flag, each polytomy is resolved into binary mergers by the node dates. The resolution is
  /// stochastic: the run logs the seed it draws, and `--seed` reproduces the run.
  #[cfg_attr(feature = "clap", clap(long, help_heading = "Polytomies"))]
  pub keep_polytomies: bool,

  /// Coalescent prior: a fixed time scale in substitutions per site (as in TreeTime 0.x and augur, for example 0.01)
  /// or as a time (10y, 26w, 90d); opt fits a constant time scale; skyline fits a piecewise-constant one
  ///
  /// A time is a positive number directly followed by the unit y (years), w (weeks) or d (days). A time scale in
  /// substitutions per site is converted to years with the clock rate of each iteration. Without this setting the node
  /// times have no coalescent prior.
  #[schemars(example = 0.01)]
  #[cfg_attr(
    feature = "clap",
    clap(long, value_parser = parse_coalescent, value_name = "TC|opt|skyline", help_heading = "Coalescent prior")
  )]
  pub coalescent: Option<CoalescentArg>,

  /// Number of grid points of the skyline coalescent model
  ///
  /// Applies only to `--coalescent skyline`. Sets how many piecewise linear segments model Tc(t) over time. Must be at
  /// least 2.
  #[cfg_attr(feature = "clap", clap(long, value_name = "N", help_heading = "Coalescent prior"))]
  #[schemars(range(min = N_SKYLINE_MIN), transform = documented_default(&N_SKYLINE_DEFAULT))]
  pub n_skyline: Option<usize>,

  /// Smoothing stiffness of the skyline coalescent model
  ///
  /// Penalizes log-fold-changes of the coalescent time scale Tc between adjacent
  /// skyline segments: with z = ln Tc, the objective adds
  /// `(stiffness/2) * Σ (ln(Tc_{i+1}/Tc_i))^2`. Because it acts on log Tc, the
  /// stiffness is dimensionless and scale-independent. Larger values enforce a
  /// smoother Tc(t). Applies only to `--coalescent skyline`.
  #[cfg_attr(
    feature = "clap",
    clap(long, hide_short_help = true, help_heading = "Coalescent prior")
  )]
  #[schemars(transform = documented_default(&SKYLINE_STIFFNESS_DEFAULT))]
  pub skyline_stiffness: Option<f64>,

  /// Width of the bands of the coalescent time scale Tc, in standard deviations
  ///
  /// Applies only to the inferred coalescent priors (`--coalescent opt` and `--coalescent skyline`).
  /// The band spans `Tc * exp(±confidence * σ)`, where `σ` is
  /// the standard deviation of `ln Tc` from the coalescent likelihood curvature. A fixed
  /// time scale is not inferred and therefore has no band.
  #[cfg_attr(
    feature = "clap",
    clap(long, hide_short_help = true, help_heading = "Coalescent prior")
  )]
  #[schemars(transform = documented_default(&COALESCENT_CONFIDENCE_DEFAULT))]
  pub coalescent_confidence: Option<f64>,

  /// Generations per year, to convert the coalescent time scale Tc into an effective population size
  ///
  /// The coalescent output reports an effective population size `N_e = Tc * gen_per_year`. Tc is
  /// already in calendar years, so this factor rescales it into generation units, the standard
  /// axis of a skyline plot. Only affects the reported `N_e`; it does not enter the inference. Applies only with
  /// `--coalescent`.
  #[cfg_attr(
    feature = "clap",
    clap(long, hide_short_help = true, help_heading = "Coalescent prior")
  )]
  #[schemars(transform = documented_default(&GEN_PER_YEAR_DEFAULT))]
  pub gen_per_year: Option<f64>,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub model_args: ModelArgs,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub partition_rates_args: PartitionRatesArgs,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub representation: RepresentationArgs,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub gap_fill_args: GapFillArgs,

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
  pub output: OutputCoreArgs,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub output_selection: TimetreeOutputSelectionArgs,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub branch_lengths: BranchLengthArgs,

  /// Path to output augur-compatible node data JSON.
  ///
  /// Contains per-node dates, branch lengths, clock model parameters, confidence intervals, and divergence metrics. The
  /// output is compatible with augur export v2 --node-data for Nextstrain pipeline integration.
  #[cfg_attr(feature = "clap", clap(long, hide_short_help = true, value_hint = ValueHint::FilePath, help_heading = "Output"))]
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

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub clock_output: ClockOutputArgs,

  /// Path to output node dates TSV: calendar date and numeric date of every node, with the bounds of the 90% max
  /// posterior region with `--confidence`.
  #[cfg_attr(feature = "clap", clap(long, hide_short_help = true, value_hint = ValueHint::FilePath, help_heading = "Output"))]
  #[schemars(extend("x-path" = "output"))]
  pub output_dates: Option<PathBuf>,

  /// Path to output the convergence values of every iteration as a table
  #[cfg_attr(feature = "clap", clap(long, hide_short_help = true, value_hint = ValueHint::FilePath, help_heading = "Output"))]
  #[schemars(extend("x-path" = "output"))]
  pub output_tracelog: Option<PathBuf>,

  /// Path to output the coalescent time scale as a table (one row per skyline segment).
  ///
  /// Written when a coalescent prior is set with `--coalescent`. A fixed time scale writes one band-less segment over
  /// the tree span.
  #[cfg_attr(feature = "clap", clap(long, hide_short_help = true, value_hint = ValueHint::FilePath, help_heading = "Output"))]
  #[schemars(extend("x-path" = "output"))]
  pub output_coalescent: Option<PathBuf>,

  /// Path to output the coalescent time scale as a rich JSON document (inputs + segments).
  ///
  /// Written when a coalescent model is set.
  #[cfg_attr(feature = "clap", clap(long, hide_short_help = true, value_hint = ValueHint::FilePath, help_heading = "Output"))]
  #[schemars(extend("x-path" = "output"))]
  pub output_coalescent_json: Option<PathBuf>,

  #[cfg_attr(feature = "clap", clap(flatten))]
  #[schemars(flatten)]
  #[deser(flatten)]
  pub ambiguity: TimetreeAmbiguityArgs,

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
