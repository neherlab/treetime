use crate::clock::clock_regression::ClockVarianceParams;
use crate::clock::date_constraints::DateConstraints;
use crate::clock::find_best_root::params::{BranchPointOptimizationParams, RerootSpec};
use crate::gtr::get_gtr::GtrModelName;
use crate::make_report;
use crate::optimize::params::BranchLengthMode;
use crate::progress::LogSink;
use crate::seq::alignment::get_common_length;
use crate::{progress_info, progress_warn};
use eyre::Report;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use smart_default::SmartDefault;
use std::fmt::Debug;
use treetime_primitives::AlignmentRecord;

pub(crate) fn compute_effective_time_marginal(
  time_marginal: TimeMarginalMode,
  confidence: bool,
  clock_std_dev: Option<f64>,
  covariation: bool,
  log: &dyn LogSink,
) -> TimeMarginalMode {
  if confidence && time_marginal == TimeMarginalMode::Never {
    if clock_std_dev.is_some() || covariation {
      progress_info!(
        log,
        "--confidence: promoting time-marginal from never to only-final for CI estimation"
      );
      TimeMarginalMode::OnlyFinal
    } else {
      progress_warn!(
        log,
        "Cannot estimate confidence intervals without clock rate uncertainty. \
         Specify --clock-std-dev or rerun with --covariation. \
         Proceeding without confidence estimation."
      );
      TimeMarginalMode::Never
    }
  } else {
    time_marginal
  }
}

#[derive(Copy, Debug, Clone, PartialEq, Eq, PartialOrd, Ord, SmartDefault, Serialize, Deserialize, JsonSchema)]
#[serde(rename_all = "kebab-case")]
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
}

#[expect(
  clippy::as_conversions,
  reason = "count/index numeric cast is exact for the domain range"
)]
pub(crate) fn build_covariation_clock_params(
  covariation: bool,
  sequence_length: Option<usize>,
  tip_slack: Option<f64>,
  aln: Option<&[AlignmentRecord]>,
  log: &dyn LogSink,
) -> Result<Option<ClockVarianceParams>, Report> {
  if !covariation {
    return Ok(None);
  }

  let seq_len = if let Some(aln_data) = aln {
    get_common_length(aln_data)? as f64
  } else {
    sequence_length.ok_or_else(|| make_report!("--sequence-length required for --covariation without alignment"))?
      as f64
  };

  let tip_slack = tip_slack.unwrap_or(10.0);

  progress_info!(
    log,
    "Covariation-aware clock regression: seq_len={seq_len}, tip_slack={tip_slack}"
  );

  Ok(Some(ClockVarianceParams {
    variance_factor: 1.0 / seq_len,
    variance_offset: 0.0,
    variance_offset_leaf: tip_slack * tip_slack / (seq_len * seq_len),
  }))
}

pub struct TimetreeParams {
  pub model: GtrModelName,
  pub dense: Option<bool>,
  pub branch_length_mode: BranchLengthMode,
  pub no_indels: bool,
  pub sequence_length: Option<usize>,
  pub clock_rate: Option<f64>,
  pub clock_std_dev: Option<f64>,
  pub keep_root: bool,
  pub reroot_spec: RerootSpec,
  pub allow_negative_rate: bool,
  pub clock_filter: f64,
  pub covariation: bool,
  pub tip_slack: Option<f64>,
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
  pub impute_missing_data: bool,
  pub sequence_outputs_requested: bool,
  pub seed: Option<u64>,
}

pub(crate) struct TimetreeContext {
  pub final_sequences: bool,
  pub final_marginal_update: bool,
  pub time_marginal: TimeMarginalMode,
  pub date_constraints: DateConstraints,
  pub covariation_clock_params: ClockVarianceParams,
  pub branch_params: BranchPointOptimizationParams,
}
