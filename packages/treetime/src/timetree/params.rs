use crate::clock::date_constraints::DateConstraints;
use crate::clock::prune_filtered::PruneFilteredTips;
use crate::optimize::indel::IndelPolicy;
use crate::optimize::partition_rates::PartitionRates;
use crate::reroot::params::RerootSpec;
use crate::timetree::coalescent_timescale::CoalescentMode;
use std::fmt::Debug;
use strum_macros::Display;
use treetime_grid::MaxGridPoints;

pub struct TimetreeParams {
  pub dense: Option<bool>,
  pub branch_length_mode: BranchLengthMode,
  pub indel_policy: IndelPolicy,
  pub partition_rates: PartitionRates,
  pub clock_rate: Option<f64>,
  pub clock_std_dev: Option<f64>,
  pub reroot_spec: Option<RerootSpec>,
  pub allow_negative_rate: bool,
  pub clock_filter: f64,
  pub prune_filtered_tips: Option<PruneFilteredTips>,
  pub max_iter: usize,
  pub resolve_polytomies: bool,
  pub relax: Vec<f64>,
  pub coalescent: CoalescentMode,
  pub n_skyline: usize,
  pub skyline_stiffness: f64,
  pub coalescent_confidence: f64,
  pub gen_per_year: f64,
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
  pub date_constraints: DateConstraints,
  pub sequence_length: usize,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum FinalSequences {
  Skip,
  FromLastRound,
  AfterMarginalUpdate,
}

#[derive(Copy, Debug, Clone, Display, PartialEq, Eq, PartialOrd, Ord)]
#[strum(serialize_all = "kebab-case")]
pub enum BranchLengthMode {
  Input,
  Marginal,
}
