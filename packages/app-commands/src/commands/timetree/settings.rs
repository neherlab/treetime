use crate::commands::shared::argument_error::argument_error;
use crate::commands::shared::branch_lengths::{BranchLengthArgs, indel_policy};
use crate::commands::shared::clock_filter::{ClockFilterArgs, TipPruning};
use crate::commands::shared::dependent_settings::check_applies;
use crate::commands::shared::model::{PartitionRatesArgs, check_partition_rates, resolve_partition_rates};
use crate::commands::shared::output_args::DivergenceUnits;
use crate::commands::shared::representation::{RepresentationArgs, representation_dense};
use crate::commands::shared::sequence_inputs::SequenceInputs;
use crate::commands::shared::tip_sequences::{ImputeMissingDataArgs, TipSequenceArgs};
use crate::commands::timetree::args::{
  COALESCENT_CONFIDENCE_DEFAULT, GEN_PER_YEAR_DEFAULT, N_SKYLINE_DEFAULT, N_SKYLINE_MIN, SKYLINE_STIFFNESS_DEFAULT,
};
use crate::commands::timetree::branch_length_mode::BranchLengthModeCli;
use crate::commands::timetree::coalescent::CoalescentArg;
use eyre::Report;
use treetime::gtr::get_gtr::GtrModel;
use treetime::optimize::indel::IndelPolicy;
use treetime::optimize::partition_rates::PartitionRates;
use treetime::timetree::coalescent_timescale::CoalescentMode;
use treetime_grid::MaxGridPoints;
use treetime_utils::o;

#[deny(unused_variables)]
pub(crate) fn timetree_branch_lengths(
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
pub(crate) fn timetree_sequence_model(
  model: GtrModel,
  RepresentationArgs { representation }: RepresentationArgs,
  PartitionRatesArgs { partition_rates }: PartitionRatesArgs,
  sequences: &SequenceInputs,
) -> Result<TimetreeSequenceModel, Report> {
  check_partition_rates(
    partition_rates,
    sequences.alignments.len(),
    sequences.partitions.entry_count(),
    sequences.partitions.is_scheme(),
  )?;
  Ok(TimetreeSequenceModel {
    model,
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
pub(crate) fn timetree_clock(
  clock_rate: Option<f64>,
  clock_std_dev: Option<f64>,
  confidence: bool,
  ClockFilterArgs {
    clock_filter,
    allow_negative_rate,
    prune_filtered_tips,
    prune_undated_tips,
    keep_tips,
    keep_tips_file,
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
    prune_tips: TipPruning::new(
      clock_filter,
      prune_filtered_tips,
      prune_undated_tips,
      keep_tips,
      keep_tips_file,
    )?,
    relax,
  })
}

#[derive(Debug, Clone)]
pub(crate) struct TimetreeClock {
  pub(crate) clock_rate: Option<f64>,
  pub(crate) clock_std_dev: Option<f64>,
  pub(crate) clock_filter: f64,
  pub(crate) allow_negative_rate: bool,
  pub(crate) prune_tips: TipPruning,
  pub(crate) relax: Vec<f64>,
}

#[deny(unused_variables)]
pub(crate) fn timetree_dating(
  max_iter: usize,
  keep_polytomies: bool,
  no_collapse_short_branches: bool,
  confidence: bool,
  max_grid_points: Option<MaxGridPoints>,
  coalescent: CoalescentSettings,
) -> TimetreeDating {
  TimetreeDating {
    max_iter,
    keep_polytomies,
    collapse_short_branches: !no_collapse_short_branches,
    confidence,
    max_grid_points: max_grid_points.unwrap_or_default(),
    coalescent,
  }
}

#[derive(Debug, Clone, Copy)]
pub(crate) struct TimetreeDating {
  pub(crate) max_iter: usize,
  pub(crate) keep_polytomies: bool,
  pub(crate) collapse_short_branches: bool,
  pub(crate) confidence: bool,
  pub(crate) max_grid_points: MaxGridPoints,
  pub(crate) coalescent: CoalescentSettings,
}

#[deny(unused_variables)]
pub(crate) fn timetree_tip_sequences(
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

pub(crate) fn coalescent_settings(
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
