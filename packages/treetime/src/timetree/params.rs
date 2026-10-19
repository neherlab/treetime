use crate::clock::date_constraints::DateConstraints;
use crate::clock::find_best_root::params::RerootSpec;
use crate::gtr::get_gtr::GtrModelName;
use crate::optimize::indel::IndelPolicy;
use crate::optimize::params::BranchLengthMode;
use crate::optimize::partition_rates::PartitionRates;
use crate::progress::LogSink;
use crate::progress_info;
use smart_default::SmartDefault;
use std::fmt::Debug;
use treetime_grid::MaxGridPoints;

pub(crate) fn compute_effective_time_marginal(
  time_marginal: TimeMarginalMode,
  confidence: bool,
  log: &dyn LogSink,
) -> TimeMarginalMode {
  if confidence && time_marginal == TimeMarginalMode::Never {
    progress_info!(
      log,
      "--confidence: promoting time-marginal from never to only-final for CI estimation"
    );
    TimeMarginalMode::OnlyFinal
  } else {
    time_marginal
  }
}

pub struct TimetreeParams {
  pub model: GtrModelName,
  pub dense: Option<bool>,
  pub branch_length_mode: BranchLengthMode,
  pub indel_policy: IndelPolicy,
  pub partition_rates: PartitionRates,
  pub sequence_length: Option<usize>,
  pub clock_rate: Option<f64>,
  pub clock_std_dev: Option<f64>,
  pub keep_root: bool,
  pub reroot_spec: RerootSpec,
  pub allow_negative_rate: bool,
  pub clock_filter: f64,
  pub max_iter: usize,
  pub resolve_polytomies: bool,
  pub relax: Vec<f64>,
  pub coalescent: Option<f64>,
  pub coalescent_opt: bool,
  pub coalescent_skyline: bool,
  pub skyline_n_points: usize,
  pub skyline_stiffness: f64,
  pub coalescent_confidence: f64,
  pub gen_per_year: f64,
  pub n_branches_posterior: Option<usize>,
  pub time_marginal: TimeMarginalMode,
  pub confidence: bool,
  pub include_leaves: bool,
  pub report_ambiguous: bool,
  pub impute_missing_data: bool,
  pub sequence_outputs_requested: bool,
  pub seed: u64,
  pub max_grid_points: MaxGridPoints,
}

pub(crate) struct TimetreeContext {
  pub final_sequences: FinalSequences,
  pub time_marginal: TimeMarginalMode,
  pub date_constraints: DateConstraints,
  pub sequence_length: usize,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum FinalSequences {
  Skip,
  FromLastRound,
  AfterMarginalUpdate,
}

#[derive(Copy, Debug, Clone, PartialEq, Eq, PartialOrd, Ord, SmartDefault)]
pub enum TimeMarginalMode {
  #[default]
  Never,
  Always,
  OnlyFinal,
}

impl TimeMarginalMode {
  pub(crate) fn runs_final_round(self) -> bool {
    self == Self::OnlyFinal
  }

  pub(crate) fn gives_confidence_intervals(self) -> bool {
    matches!(self, Self::OnlyFinal | Self::Always)
  }
}
