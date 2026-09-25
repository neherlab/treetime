use crate::results::tree::ResultTree;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use std::collections::BTreeSet;
use treetime::clock::clock_model::{ClockModel, ClockModelStats};
use treetime::clock::rtt::{ClockDateSource, ClockRegressionResult};
use treetime_utils::datetime::year_fraction::year_fraction_days_between;

/// Results of a `clock` run.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
pub struct ClockResults {
  /// Estimates of the clock model.
  pub estimates: ClockEstimates,
  /// Samples and line of the root-to-tip regression; absent when the run wrote no clock regression table.
  pub root_to_tip: Option<RootToTip>,
}

/// Estimates of a `clock` run.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
pub struct ClockEstimates {
  /// Clock rate in substitutions per site per year; absent when the run wrote no clock model.
  pub clock_rate: Option<f64>,
  /// Whether the clock rate was fixed by the user instead of estimated.
  pub clock_rate_fixed: bool,
  /// Correlation coefficient of the root-to-tip regression; absent for a fixed rate.
  pub r: Option<f64>,
  /// Coefficient of determination of the root-to-tip regression.
  pub r_squared: Option<f64>,
  /// Number of samples with a date.
  pub dated_samples: usize,
  /// Number of dated samples the clock filter flagged as outliers.
  pub outliers: usize,
}

/// The points and line of a root-to-tip regression, as TreeTime fitted it.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
pub struct RootToTip {
  /// Samples as the regression saw them.
  pub points: Vec<RootToTipPoint>,
  /// Line of the clock model; absent when the run wrote no clock model.
  pub line: Option<ClockLine>,
}

/// Line of a clock model: divergence = rate * date + intercept.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
pub struct ClockLine {
  /// Clock rate in substitutions per site per year.
  pub rate: f64,
  /// Divergence at year 0.
  pub intercept: f64,
}

/// One sample of a root-to-tip regression.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
pub struct RootToTipPoint {
  /// Name of the sample.
  pub name: String,
  /// Date the regression used, as a decimal year.
  pub date: Option<f64>,
  /// Where the date came from; absent when the table does not say.
  pub date_source: Option<ClockDateSource>,
  /// Root-to-tip divergence the regression used.
  pub div: f64,
  /// Date the clock model predicts from the divergence.
  pub predicted_date: f64,
  /// Sampling date minus the predicted date, in days.
  pub residual_days: Option<f64>,
  /// Whether the clock filter flagged the sample as an outlier.
  pub outlier: bool,
}

pub fn clock_results(
  tree: Option<&ResultTree>,
  model: Option<&ClockModel>,
  rows: Option<&[ClockRegressionResult]>,
) -> ClockResults {
  let root_to_tip = rows.map(|rows| root_to_tip(tree, model, rows));
  let dated = root_to_tip
    .iter()
    .flat_map(|regression| &regression.points)
    .filter(|point| point.date.is_some())
    .collect::<Vec<_>>();
  ClockResults {
    estimates: ClockEstimates {
      clock_rate: model.map(ClockModel::clock_rate),
      clock_rate_fixed: model.is_some_and(is_fixed),
      r: model.and_then(ClockModel::r_val),
      r_squared: model.and_then(ClockModel::r_val).map(|r| r * r),
      dated_samples: dated.len(),
      outliers: dated.iter().filter(|point| point.outlier).count(),
    },
    root_to_tip,
  }
}

pub fn root_to_tip(tree: Option<&ResultTree>, model: Option<&ClockModel>, rows: &[ClockRegressionResult]) -> RootToTip {
  let tips: Option<BTreeSet<&str>> = tree.map(|tree| tree.tips().map(|tip| tip.name.as_str()).collect());
  let points = rows
    .iter()
    .filter_map(|row| {
      let name = row.name.as_deref()?;
      let listed = tips
        .as_ref()
        .map_or_else(|| row.date.is_some(), |tips| tips.contains(name));
      listed.then(|| RootToTipPoint {
        name: name.to_owned(),
        date: row.date,
        date_source: row.date_source,
        div: row.div,
        predicted_date: row.predicted_date,
        residual_days: row
          .date
          .map(|date| year_fraction_days_between(row.predicted_date, date)),
        outlier: row.is_outlier,
      })
    })
    .collect();
  RootToTip {
    points,
    line: model.map(|model| ClockLine {
      rate: model.clock_rate(),
      intercept: model.intercept(),
    }),
  }
}

pub fn is_fixed(model: &ClockModel) -> bool {
  matches!(model.stats(), ClockModelStats::Fixed)
}
