use crate::config::catalog::{SettingRole, command_settings};
use crate::job::JobId;
use crate::results::clades::matched_ancestors;
use crate::results::run_results::{CommandResults, RunResults, run_results};
use crate::results::timetree::TimetreeEstimates;
use crate::results::tree::ResultTree;
use crate::results::year_date::YearDate;
use crate::runs::manager::RunManager;
use crate::runs::record::{RunRecord, RunStatus};
use crate::runs::setting_differences::{SettingDifference, setting_differences};
use deser::Serialize;
use eyre::Report;
use schemars::JsonSchema;
use treetime_schema::skip_serializing_optionals;
use treetime_utils::datetime::year_fraction::year_fraction_days_between;

/// Comparison of two runs: their settings and their results.
#[derive(Clone, Debug, PartialEq, JsonSchema, Serialize)]
#[deser(skip_serializing_optionals)]
#[schemars(transform = skip_serializing_optionals)]
pub struct RunComparison {
  /// Settings and inputs that differ; absent when the runs execute different commands.
  pub settings: Option<SettingsComparison>,
  /// Estimates side by side; present when both runs are finished time-tree runs with a tree.
  pub estimates: Option<EstimateComparison>,
  /// Date shifts of the ancestors both trees share; present when both runs are finished time-tree runs with a tree.
  pub ancestors: Option<AncestorComparison>,
}

/// Settings and inputs that differ between two runs of the same command.
#[derive(Clone, Debug, PartialEq, Eq, JsonSchema, Serialize)]
pub struct SettingsComparison {
  /// Settings and inputs whose values differ.
  pub differences: Vec<SettingDifference>,
  /// Number of settings and inputs compared.
  pub compared: usize,
  /// Whether both runs have the same configuration hash: the same settings on the same input contents.
  pub same_config_hash: bool,
}

/// Estimates of two time-tree runs and their differences, second minus first.
#[derive(Clone, Debug, PartialEq, JsonSchema, Serialize)]
#[deser(skip_serializing_optionals)]
#[schemars(transform = skip_serializing_optionals)]
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
#[derive(Clone, Debug, PartialEq, JsonSchema, Serialize)]
#[deser(skip_serializing_optionals)]
#[schemars(transform = skip_serializing_optionals)]
pub struct AncestorComparison {
  /// Ancestors dated in both trees, with the shift of their date.
  pub shifts: Vec<AncestorShift>,
  /// Number of ancestors in the first tree.
  pub ancestors: usize,
  /// Mean absolute shift, in days.
  pub mean_absolute_shift_days: Option<f64>,
}

/// Date shift of one ancestor between two trees.
#[derive(Clone, Debug, PartialEq, JsonSchema, Serialize)]
pub struct AncestorShift {
  /// Name of the ancestor in the first tree.
  pub name: String,
  /// Number of samples below the ancestor.
  pub tips: usize,
  /// Date in the first tree.
  pub date_first: YearDate,
  /// Date in the second tree minus the date in the first, in days.
  pub shift_days: f64,
}

pub fn compare_runs(manager: &RunManager, first: &JobId, second: &JobId) -> Result<RunComparison, Report> {
  let first_record = manager.get(first)?;
  let second_record = manager.get(second)?;
  let settings = compare_settings(&first_record, &second_record)?;
  if first_record.status == RunStatus::Ok && second_record.status == RunStatus::Ok {
    Ok(compare_results(
      &run_results(manager, first)?,
      &run_results(manager, second)?,
      settings,
    ))
  } else {
    Ok(RunComparison {
      settings,
      estimates: None,
      ancestors: None,
    })
  }
}

fn compare_settings(first: &RunRecord, second: &RunRecord) -> Result<Option<SettingsComparison>, Report> {
  if first.config.command() != second.config.command() {
    return Ok(None);
  }
  let compared = command_settings(first.config.command())?
    .settings
    .iter()
    .filter(|spec| spec.role != SettingRole::Output)
    .count();
  Ok(Some(SettingsComparison {
    differences: setting_differences(first, second)?,
    compared,
    same_config_hash: first.config_hash.is_some() && first.config_hash == second.config_hash,
  }))
}

fn compare_results(first: &RunResults, second: &RunResults, settings: Option<SettingsComparison>) -> RunComparison {
  let timetrees = match (&first.results, &second.results) {
    (CommandResults::Timetree(a), CommandResults::Timetree(b)) => Some((a, b)),
    _ => None,
  };
  RunComparison {
    settings,
    ancestors: timetrees
      .and_then(|_| first.tree.as_ref().zip(second.tree.as_ref()))
      .map(|(a, b)| compare_ancestors(a, b)),
    estimates: timetrees
      .and_then(|(a, b)| a.estimates.clone().zip(b.estimates.clone()))
      .map(|(a, b)| compare_estimates(a, b)),
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
      .as_ref()
      .zip(second.root_date.as_ref())
      .map(|(a, b)| year_fraction_days_between(a.year, b.year)),
    root_interval_change_days: first
      .root_interval
      .as_ref()
      .zip(second.root_interval.as_ref())
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
      let date_first = node.date.as_ref()?;
      let date_second = second.nodes[b].date.as_ref()?;
      Some(AncestorShift {
        name: node.name.clone(),
        tips: node.tips,
        shift_days: year_fraction_days_between(date_first.year, date_second.year),
        date_first: date_first.clone(),
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
