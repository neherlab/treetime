use crate::job::JobId;
use crate::results::clades::matched_ancestors;
use crate::results::run_results::{CommandResults, RunResults, run_results};
use crate::results::timetree::TimetreeEstimates;
use crate::results::tree::ResultTree;
use crate::runs::manager::RunManager;
use crate::runs::record::RunStatus;
use eyre::Report;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use treetime_utils::datetime::year_fraction::year_fraction_days_between;

/// Comparison of the results of two runs.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
pub struct RunComparison {
  /// Estimates side by side; present when both runs are finished time-tree runs with a tree.
  pub estimates: Option<EstimateComparison>,
  /// Date shifts of the ancestors both trees share; present when both runs are finished time-tree runs with a tree.
  pub ancestors: Option<AncestorComparison>,
}

/// Estimates of two time-tree runs and their differences, second minus first.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
pub struct EstimateComparison {
  /// Estimates of the first run.
  pub first: TimetreeEstimates,
  /// Estimates of the second run.
  pub second: TimetreeEstimates,
  /// Shift of the root date, in days.
  pub root_shift_days: Option<f64>,
  /// Change of the width of the root-date interval, in days.
  pub root_interval_change_days: Option<f64>,
  /// Change of the clock rate, in percent of the first rate.
  pub clock_rate_change_percent: Option<f64>,
  /// Change of the number of samples the clock model left out.
  pub excluded_samples_change: i64,
  /// Change of the total log likelihood; absent when either value is missing or not finite.
  pub log_likelihood_change: Option<f64>,
}

/// Date shifts of the ancestors two trees share, matched by their set of samples.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
pub struct AncestorComparison {
  /// Ancestors dated in both trees, with the shift of their date.
  pub shifts: Vec<AncestorShift>,
  /// Number of ancestors in the first tree.
  pub ancestors: usize,
  /// Mean absolute shift, in days.
  pub mean_absolute_shift_days: Option<f64>,
}

/// Date shift of one ancestor between two trees.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
pub struct AncestorShift {
  /// Name of the ancestor in the first tree.
  pub name: String,
  /// Number of samples below the ancestor.
  pub tips: usize,
  /// Date in the first tree, as a decimal year.
  pub date_first: f64,
  /// Date in the second tree minus the date in the first, in days.
  pub shift_days: f64,
}

pub fn compare_runs(manager: &RunManager, first: &JobId, second: &JobId) -> Result<RunComparison, Report> {
  let finished_timetree = |id: &JobId| -> Result<Option<RunResults>, Report> {
    let record = manager.get(id)?;
    if record.status == RunStatus::Ok {
      run_results(manager, id).map(Some)
    } else {
      Ok(None)
    }
  };
  match (finished_timetree(first)?, finished_timetree(second)?) {
    (Some(first), Some(second)) => Ok(compare_results(&first, &second)),
    _ => Ok(RunComparison {
      estimates: None,
      ancestors: None,
    }),
  }
}

pub fn compare_results(first: &RunResults, second: &RunResults) -> RunComparison {
  let estimates = match (&first.results, &second.results) {
    (CommandResults::Timetree(a), CommandResults::Timetree(b)) => a.estimates.clone().zip(b.estimates.clone()),
    _ => None,
  };
  let timetrees = matches!(
    (&first.results, &second.results),
    (CommandResults::Timetree(_), CommandResults::Timetree(_))
  );
  RunComparison {
    ancestors: first
      .tree
      .as_ref()
      .zip(second.tree.as_ref())
      .filter(|_| timetrees)
      .map(|(a, b)| compare_ancestors(a, b)),
    estimates: estimates.map(|(a, b)| compare_estimates(a, b)),
  }
}

pub fn compare_estimates(first: TimetreeEstimates, second: TimetreeEstimates) -> EstimateComparison {
  let finite = |estimates: &TimetreeEstimates| {
    estimates
      .log_likelihood
      .map(|value| value.0)
      .filter(|value| value.is_finite())
  };
  EstimateComparison {
    root_shift_days: first
      .root_date
      .zip(second.root_date)
      .map(|(a, b)| year_fraction_days_between(a, b)),
    root_interval_change_days: first
      .root_interval
      .zip(second.root_interval)
      .map(|(a, b)| b.days - a.days),
    clock_rate_change_percent: first
      .clock_rate
      .zip(second.clock_rate)
      .filter(|(a, _)| *a != 0.0)
      .map(|(a, b)| (b - a) / a * 100.0),
    excluded_samples_change: count(second.excluded_samples) - count(first.excluded_samples),
    log_likelihood_change: finite(&first).zip(finite(&second)).map(|(a, b)| b - a),
    first,
    second,
  }
}

pub fn compare_ancestors(first: &ResultTree, second: &ResultTree) -> AncestorComparison {
  let shifts = matched_ancestors(first, second)
    .into_iter()
    .filter_map(|(a, b)| {
      let node = &first.nodes[a];
      let date_first = node.date?;
      let date_second = second.nodes[b].date?;
      Some(AncestorShift {
        name: node.name.clone(),
        tips: node.tips,
        date_first,
        shift_days: year_fraction_days_between(date_first, date_second),
      })
    })
    .collect::<Vec<_>>();
  #[allow(
    clippy::as_conversions,
    clippy::cast_precision_loss,
    reason = "an ancestor count is far below 2^52"
  )]
  let mean_absolute_shift_days =
    (!shifts.is_empty()).then(|| shifts.iter().map(|shift| shift.shift_days.abs()).sum::<f64>() / shifts.len() as f64);
  AncestorComparison {
    ancestors: first.nodes.iter().filter(|node| !node.is_tip()).count(),
    mean_absolute_shift_days,
    shifts,
  }
}

#[allow(
  clippy::as_conversions,
  clippy::cast_possible_wrap,
  reason = "a sample count is far below i64::MAX"
)]
fn count(value: usize) -> i64 {
  value as i64
}
