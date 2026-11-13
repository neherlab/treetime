use crate::results::clock::{RootToTipPoint, has_fixed_rate};
use crate::results::year_date::YearDate;
use deser::Serialize;
use schemars::JsonSchema;
use treetime::clock::clock_model::ClockModel;
use treetime_schema::skip_serializing_optionals;
use treetime_utils::cast::ToF64Exact;

const CLOCK_RATE_UNIT: &str = "substitutions per site per year";

const RESIDUAL_MEAN_SQUARE_UNIT: &str = "(substitutions per site)²";

pub fn regression_statistics(model: &ClockModel, points: &[RootToTipPoint]) -> Vec<RegressionStatistic> {
  let number = |value: f64| StatisticValue::Number { value };
  [
    Some((RegressionStatisticKind::ClockRate, number(model.clock_rate()))),
    Some((
      RegressionStatisticKind::RootDate,
      StatisticValue::Date {
        date: YearDate::new(model.t_mrca()),
      },
    )),
    model
      .r_val()
      .map(|r| (RegressionStatisticKind::CorrelationCoefficient, number(r))),
    model
      .r_squared()
      .map(|r_squared| (RegressionStatisticKind::RSquared, number(r_squared))),
    residual_mean_square(model, points).map(|value| (RegressionStatisticKind::ResidualMeanSquare, number(value))),
  ]
  .into_iter()
  .flatten()
  .map(|(statistic, value)| RegressionStatistic {
    statistic,
    label: statistic.label().to_owned(),
    definition: statistic.definition().to_owned(),
    unit: statistic.unit().map(str::to_owned),
    value,
  })
  .collect()
}

/// One statistic of a root-to-tip regression, as Tempest lists them.
#[derive(Clone, Debug, PartialEq, JsonSchema, Serialize)]
#[deser(skip_serializing_optionals)]
#[schemars(transform = skip_serializing_optionals)]
pub struct RegressionStatistic {
  /// Which statistic this is.
  pub statistic: RegressionStatisticKind,
  /// Name of the statistic, for example `Correlation coefficient r`.
  pub label: String,
  /// What the statistic measures, in one or two sentences.
  pub definition: String,
  /// Unit of a numeric value; absent for a value without a unit.
  pub unit: Option<String>,
  /// Value of the statistic.
  pub value: StatisticValue,
}

/// A statistic of a root-to-tip regression.
#[derive(Clone, Copy, Debug, PartialEq, Eq, JsonSchema, Serialize)]
#[schemars(rename_all = "kebab-case")]
#[deser(rename_all = "kebab-case")]
pub enum RegressionStatisticKind {
  /// Slope of the regression line.
  ClockRate,
  /// Date at which the regression line reaches zero divergence.
  RootDate,
  /// Correlation coefficient of sampling date and divergence.
  CorrelationCoefficient,
  /// Square of the correlation coefficient.
  RSquared,
  /// Mean square of the divergence residuals.
  ResidualMeanSquare,
}

impl RegressionStatisticKind {
  const fn label(self) -> &'static str {
    self.facts().0
  }

  const fn definition(self) -> &'static str {
    self.facts().1
  }

  const fn unit(self) -> Option<&'static str> {
    self.facts().2
  }

  #[rustfmt::skip]
  const fn facts(self) -> (&'static str, &'static str, Option<&'static str>) {
    match self {
      Self::ClockRate              => ("Clock rate",                "Slope of the root-to-tip regression line: the divergence that accumulates per year. It is the estimated rate, or the rate given in the settings when the rate is fixed.", Some(CLOCK_RATE_UNIT)),
      Self::RootDate               => ("Root date (x-intercept)",   "Date at which the regression line reaches zero root-to-tip divergence: the age of the root that the regression implies, before the dates of the nodes are inferred.", None),
      Self::CorrelationCoefficient => ("Correlation coefficient r", "Pearson correlation of sampling date and root-to-tip divergence over the samples of the regression, from -1 to 1. Positive values mean that later samples have diverged more.", None),
      Self::RSquared               => ("R²",                        "Square of the correlation coefficient: the fraction of the variance of the root-to-tip divergence that the sampling dates explain.", None),
      Self::ResidualMeanSquare     => ("Residual mean square",      "Sum of the squared differences between the divergence of each sample and the regression line, divided by the residual degrees of freedom: the number of samples in the regression minus the fitted parameters (two, or one when the rate is fixed).", Some(RESIDUAL_MEAN_SQUARE_UNIT)),
    }
  }
}

/// Value of a regression statistic.
#[derive(Clone, Debug, PartialEq, JsonSchema, Serialize)]
#[schemars(tag = "kind", rename_all = "kebab-case")]
#[deser(tag = "kind", rename_all = "kebab-case")]
pub enum StatisticValue {
  /// A number.
  Number {
    /// The number.
    value: f64,
  },
  /// A date.
  Date {
    /// The date.
    date: YearDate,
  },
}

fn residual_mean_square(model: &ClockModel, points: &[RootToTipPoint]) -> Option<f64> {
  let residuals = points
    .iter()
    .filter(|point| !point.outlier)
    .filter_map(|point| Some(point.div - model.div(point.date.as_ref()?.year)))
    .collect::<Vec<_>>();
  let fitted_parameters = if has_fixed_rate(model) { 1 } else { 2 };
  let degrees_of_freedom = residuals.len().checked_sub(fitted_parameters).filter(|dof| *dof > 0)?;
  Some(residuals.iter().map(|residual| residual * residual).sum::<f64>() / degrees_of_freedom.to_f64_exact())
}
