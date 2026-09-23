use crate::clock::clock_set::ClockSet;
use crate::make_error;
use eyre::Report;
use getset::Getters;
use log::{debug, warn};
use ndarray::Array2;
use serde::{Deserialize, Serialize};
use treetime_utils::array::serde::{array2_as_vec, array2_from_vec};
use treetime_utils::fmt::float::float_to_significant_digits;
use treetime_utils::io::json::{JsonPretty, json_write_str};

pub trait ClockLine {
  fn clock_rate(&self) -> f64;
  fn intercept(&self) -> f64;
  fn clock_deviation(&self, date: f64, div: f64) -> f64;
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RegressionStats {
  pub chisq: f64,
  pub r_val: f64,
  #[serde(serialize_with = "array2_as_vec", deserialize_with = "array2_from_vec")]
  pub hessian: Array2<f64>,
  #[serde(serialize_with = "array2_as_vec", deserialize_with = "array2_from_vec")]
  pub cov: Array2<f64>,
}

#[must_use]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ClockRegression {
  clock_rate: f64,
  intercept: f64,
  chisq: f64,
  r_val: f64,
  #[serde(serialize_with = "array2_as_vec", deserialize_with = "array2_from_vec")]
  hessian: Array2<f64>,
  #[serde(serialize_with = "array2_as_vec", deserialize_with = "array2_from_vec")]
  cov: Array2<f64>,
}

impl ClockLine for ClockRegression {
  fn clock_rate(&self) -> f64 {
    self.clock_rate
  }

  fn intercept(&self) -> f64 {
    self.intercept
  }

  fn clock_deviation(&self, date: f64, div: f64) -> f64 {
    date * self.clock_rate + self.intercept - div
  }
}

#[expect(
  clippy::same_name_method,
  reason = "the inherent accessor and the trait method return the same field"
)]
impl ClockRegression {
  pub(crate) fn clock_rate(&self) -> f64 {
    self.clock_rate
  }

  pub(crate) fn intercept(&self) -> f64 {
    self.intercept
  }

  pub(crate) fn chisq(&self) -> f64 {
    self.chisq
  }

  pub(crate) fn r_val(&self) -> f64 {
    self.r_val
  }

  pub(crate) fn from_clock_set(clock_set: &ClockSet) -> Result<Self, Report> {
    let det = clock_set.determinant();
    if det <= 0.0 {
      debug!("ClockSet: {}", json_write_str(clock_set, JsonPretty(true))?);
      debug!("ClockSet determinant: {det}");
      return make_error!("No variation in sampling dates! Please specify your clock rate explicitly.");
    }

    let clock_rate = clock_set.clock_rate(det);
    Ok(Self {
      clock_rate,
      intercept: clock_set.intercept(clock_rate),
      chisq: clock_set.chisq(),
      r_val: clock_set.r_val(),
      hessian: clock_set.hessian(),
      cov: clock_set.cov(),
    })
  }

  pub(crate) fn hessian(&self) -> &Array2<f64> {
    &self.hessian
  }

  pub fn cov(&self) -> &Array2<f64> {
    &self.cov
  }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum ClockModelStats {
  Estimated(RegressionStats),
  Fixed,
}

#[must_use]
#[derive(Debug, Clone, Serialize, Deserialize, Getters)]
pub struct ClockModel {
  clock_rate: f64,

  intercept: f64,

  #[getset(get = "pub")]
  stats: ClockModelStats,
}

#[expect(
  clippy::same_name_method,
  reason = "the inherent accessor and the trait method return the same field"
)]
impl ClockModel {
  pub fn clock_rate(&self) -> f64 {
    self.clock_rate
  }

  pub fn intercept(&self) -> f64 {
    self.intercept
  }

  pub fn from_regression(regression: &ClockRegression) -> Result<Self, Report> {
    if regression.clock_rate <= 0.0 {
      return make_error!(
        "Estimated clock rate is non-positive ({:.6e}).\n\n\
         This means the root-to-tip regression found no positive correlation between \
         sampling dates and genetic divergence, which prevents time-scaled analysis.\n\n\
         Suggestions:\n\
         - Specify a known substitution rate with --clock-rate\n\
         - Verify that sampling dates are correct and span a sufficient time range\n\
         - Check that the alignment has enough informative sites",
        regression.clock_rate
      );
    }

    Ok(Self::from_regression_unchecked(regression))
  }

  pub(crate) fn from_regression_allow_negative(regression: &ClockRegression) -> Self {
    if regression.clock_rate <= 0.0 {
      warn!(
        "Estimated clock rate is non-positive ({:.6e}). The root-to-tip regression found no positive \
         correlation between sampling dates and genetic divergence. Continuing, but the dates lack a \
         reliable temporal signal; interpret the clock results with caution or specify a known rate with \
         --clock-rate.",
        regression.clock_rate
      );
    }
    Self::from_regression_unchecked(regression)
  }

  fn from_regression_unchecked(regression: &ClockRegression) -> Self {
    Self {
      clock_rate: regression.clock_rate,
      intercept: regression.intercept,
      stats: ClockModelStats::Estimated(RegressionStats {
        chisq: regression.chisq,
        r_val: regression.r_val,
        hessian: regression.hessian.clone(),
        cov: regression.cov.clone(),
      }),
    }
  }

  pub(crate) fn with_fixed_rate(clock_set: &ClockSet, clock_rate: f64) -> Result<Self, Report> {
    if clock_rate <= 0.0 {
      return make_error!(
        "Specified clock rate must be positive, got {clock_rate:.6e}.\n\n\
         The clock rate is the expected number of substitutions per site per year.\n\n\
         Suggestions:\n\
         - Provide a positive value, e.g. --clock-rate=0.001\n\
         - Omit --clock-rate to let the rate be estimated from the data"
      );
    }
    Ok(Self {
      clock_rate,
      intercept: clock_set.intercept(clock_rate),
      stats: ClockModelStats::Fixed,
    })
  }

  fn get_regression_stat<T>(&self, f: impl FnOnce(&RegressionStats) -> T) -> Option<T> {
    match &self.stats {
      ClockModelStats::Estimated(stats) => Some(f(stats)),
      ClockModelStats::Fixed => None,
    }
  }

  pub fn chisq(&self) -> Option<f64> {
    self.get_regression_stat(|s| s.chisq)
  }

  pub fn r_val(&self) -> Option<f64> {
    self.get_regression_stat(|s| s.r_val)
  }

  pub fn hessian(&self) -> Option<&Array2<f64>> {
    match &self.stats {
      ClockModelStats::Estimated(stats) => Some(&stats.hessian),
      ClockModelStats::Fixed => None,
    }
  }

  pub fn cov(&self) -> Option<&Array2<f64>> {
    match &self.stats {
      ClockModelStats::Estimated(stats) => Some(&stats.cov),
      ClockModelStats::Fixed => None,
    }
  }

  pub(crate) fn date(&self, div: f64) -> f64 {
    (div - self.intercept()) / self.clock_rate()
  }

  pub fn div(&self, date: f64) -> f64 {
    date * self.clock_rate() + self.intercept()
  }

  pub fn t_mrca(&self) -> f64 {
    self.date(0.0)
  }

  pub fn equation_str(&self) -> String {
    format!(
      "div = {:}t {:} {:}",
      float_to_significant_digits(self.clock_rate(), 3),
      if self.intercept() < 0.0 { "-" } else { "+" },
      float_to_significant_digits(self.intercept().abs(), 3)
    )
  }

  #[cfg(test)]
  pub(crate) fn for_testing(clock_rate: f64, intercept: f64) -> Self {
    Self {
      clock_rate,
      intercept,
      stats: ClockModelStats::Fixed,
    }
  }

  #[cfg(test)]
  pub(crate) fn for_testing_with_stats(clock_rate: f64, intercept: f64, stats: ClockModelStats) -> Self {
    Self {
      clock_rate,
      intercept,
      stats,
    }
  }
}

impl ClockLine for ClockModel {
  fn clock_rate(&self) -> f64 {
    self.clock_rate
  }

  fn intercept(&self) -> f64 {
    self.intercept
  }

  fn clock_deviation(&self, date: f64, div: f64) -> f64 {
    date * self.clock_rate + self.intercept - div
  }
}
