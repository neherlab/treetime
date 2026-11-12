use crate::config::catalog::{SettingPurpose, command_settings};
use crate::job::JobId;
use crate::results::clades::{matched_ancestors, same_topology};
use crate::results::run_results::{CommandResults, RunResults, run_results};
use crate::results::timetree::TimetreeEstimates;
use crate::results::tree::ResultTree;
use crate::results::year_date::YearDate;
use crate::runs::manager::RunManager;
use crate::runs::record::{RunRecord, RunStatus};
use crate::runs::setting_differences::{
  ComparedConfig, SettingDifference, compared_inputs, compared_settings, setting_differences,
};
use deser::Serialize;
use eyre::Report;
use schemars::JsonSchema;
use treetime_schema::skip_serializing_optionals;
use treetime_utils::cast::ToF64Exact;
use treetime_utils::datetime::year_fraction::year_fraction_days_between;

const SUBSTITUTION_MODEL_GROUP: &str = "Substitution model";

pub fn compare_runs(manager: &RunManager, first: &JobId, second: &JobId) -> Result<RunComparison, Report> {
  let first_record = manager.get(first)?;
  let second_record = manager.get(second)?;
  let same_analysis = first_record.config.command() == second_record.config.command();
  let settings = compare_settings(&first_record, &second_record)?;
  if first_record.status == RunStatus::Ok && second_record.status == RunStatus::Ok {
    Ok(compare_results(
      &run_results(manager, first)?,
      &run_results(manager, second)?,
      same_analysis,
      settings,
    ))
  } else {
    Ok(RunComparison {
      same_analysis,
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
  let catalog = command_settings(first.config.command())?;
  let first_inputs = compared_inputs(&first.inputs);
  let second_inputs = compared_inputs(&second.inputs);
  let differences = setting_differences(
    &catalog,
    &ComparedConfig {
      settings: &first.config.settings()?,
      inputs: &first_inputs,
    },
    &ComparedConfig {
      settings: &second.config.settings()?,
      inputs: &second_inputs,
    },
  );
  Ok(Some(SettingsComparison {
    differences,
    compared: compared_settings(&catalog).count(),
    same_config_hash: first.config_hash.is_some() && first.config_hash == second.config_hash,
  }))
}

fn compare_results(
  first: &RunResults,
  second: &RunResults,
  same_analysis: bool,
  settings: Option<SettingsComparison>,
) -> RunComparison {
  let timetrees = match (&first.results, &second.results) {
    (CommandResults::Timetree(a), CommandResults::Timetree(b)) => Some((a, b)),
    _ => None,
  };
  let trees = first.tree.as_ref().zip(second.tree.as_ref());
  let differences = settings
    .as_ref()
    .map(|settings| settings.differences.as_slice())
    .unwrap_or_default();
  RunComparison {
    same_analysis,
    ancestors: timetrees.and(trees).map(|(a, b)| compare_ancestors(a, b)),
    estimates: timetrees
      .and_then(|(a, b)| a.estimates.clone().zip(b.estimates.clone()))
      .map(|(a, b)| {
        let mismatches = likelihood_mismatches(&a, &b, differences, trees);
        compare_estimates(a, b, mismatches)
      }),
    settings,
  }
}

/// Comparison of two runs: their settings and their results.
#[derive(Clone, Debug, PartialEq, JsonSchema, Serialize)]
#[deser(skip_serializing_optionals)]
#[schemars(transform = skip_serializing_optionals)]
pub struct RunComparison {
  /// Whether both runs are the same analysis, that is, they execute the same command. Runs of different analyses have
  /// no settings comparison and no estimates.
  pub same_analysis: bool,
  /// Settings and inputs that differ; absent when the runs are different analyses.
  pub settings: Option<SettingsComparison>,
  /// Estimates side by side; present when both runs are finished time-tree runs with a tree.
  pub estimates: Option<EstimateComparison>,
  /// Date shifts of the ancestors both trees share; present when both runs are finished time-tree runs with a tree.
  pub ancestors: Option<AncestorComparison>,
}

/// Settings and inputs that differ between two runs of the same command. The settings that the app sets for every
/// run, the output files and their paths, are not compared.
#[derive(Clone, Debug, PartialEq, Eq, JsonSchema, Serialize)]
pub struct SettingsComparison {
  /// Settings whose resolved values differ, and inputs whose file contents differ.
  pub differences: Vec<SettingDifference>,
  /// Number of settings and inputs compared.
  pub compared: usize,
  /// Whether both runs have the same configuration hash: the same settings on the same input contents.
  pub same_config_hash: bool,
}

pub fn compare_estimates(
  first: TimetreeEstimates,
  second: TimetreeEstimates,
  log_likelihood_mismatches: Vec<LikelihoodMismatch>,
) -> EstimateComparison {
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
      .and_then(|(a, b)| year_fraction_days_between(a.year, b.year)),
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
    log_likelihood_change: finite(&first)
      .zip(finite(&second))
      .filter(|_| log_likelihood_mismatches.is_empty())
      .map(|(a, b)| b - a),
    log_likelihood_mismatches,
    first,
    second,
  }
}

pub fn likelihood_mismatches(
  first: &TimetreeEstimates,
  second: &TimetreeEstimates,
  differences: &[SettingDifference],
  trees: Option<(&ResultTree, &ResultTree)>,
) -> Vec<LikelihoodMismatch> {
  let differ_in = |group: &str| differences.iter().any(|difference| difference.group == group);
  [
    (
      differences
        .iter()
        .any(|difference| difference.purpose == SettingPurpose::Data),
      LikelihoodMismatch::Data,
    ),
    (
      differ_in(SUBSTITUTION_MODEL_GROUP),
      LikelihoodMismatch::SubstitutionModel,
    ),
    (
      first.coalescent_prior != second.coalescent_prior,
      LikelihoodMismatch::CoalescentPrior,
    ),
    (
      trees.is_none_or(|(a, b)| !same_topology(a, b)),
      LikelihoodMismatch::Topology,
    ),
  ]
  .into_iter()
  .filter_map(|(differs, mismatch)| differs.then_some(mismatch))
  .collect()
}

/// A difference between two time-tree runs that changes what their log likelihoods measure, so that their difference
/// is no measure of fit.
#[derive(Clone, Copy, Debug, PartialEq, Eq, JsonSchema, Serialize)]
#[schemars(rename_all = "kebab-case")]
#[deser(rename_all = "kebab-case")]
pub enum LikelihoodMismatch {
  /// The runs read different input data: different file contents, or different settings of how the files are read.
  Data,
  /// The runs use different substitution models.
  SubstitutionModel,
  /// The runs use different coalescent priors.
  CoalescentPrior,
  /// The trees of the runs differ in topology, for example because the runs resolved polytomies with different
  /// random seeds or polytomy settings: each node that the resolution adds lowers the log likelihood. Also when a run
  /// wrote no tree.
  Topology,
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
  /// Change of the total log likelihood; absent when either value is missing or not finite, or when
  /// `log_likelihood_mismatches` is not empty.
  pub log_likelihood_change: Option<f64>,
  /// Why the log likelihoods of the two runs cannot be compared; empty when they can.
  pub log_likelihood_mismatches: Vec<LikelihoodMismatch>,
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
        shift_days: year_fraction_days_between(date_first.year, date_second.year)?,
        date_first: date_first.clone(),
      })
    })
    .collect::<Vec<_>>();
  let mean_absolute_shift_days = (!shifts.is_empty())
    .then(|| shifts.iter().map(|shift| shift.shift_days.abs()).sum::<f64>() / shifts.len().to_f64_exact());
  AncestorComparison {
    ancestors: first.nodes.iter().filter(|node| !node.is_tip()).count(),
    mean_absolute_shift_days,
    shifts,
  }
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

#[expect(
  clippy::as_conversions,
  clippy::cast_possible_wrap,
  reason = "a sample count is far below i64::MAX"
)]
fn count(value: usize) -> i64 {
  value as i64
}
