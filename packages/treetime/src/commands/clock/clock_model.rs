use crate::commands::clock::clock_set::ClockSet;
use crate::make_error;
use eyre::Report;
use getset::{CopyGetters, Getters};
use log::debug;
use ndarray::Array2;
use serde::{Deserialize, Serialize};
use treetime_io::json::{JsonPretty, json_write_str};
use treetime_utils::fmt::float::float_to_significant_digits;

/// Regression statistics from clock model estimation
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RegressionStats {
  pub chisq: f64,
  pub r_val: f64,
  pub hessian: Array2<f64>,
  pub cov: Array2<f64>,
}

/// Clock model statistics - either estimated from data or fixed by user
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum ClockModelStats {
  /// Clock rate estimated from root-to-tip regression with full statistics
  Estimated(RegressionStats),
  /// Clock rate fixed by user input (no regression statistics available)
  Fixed,
}

#[must_use]
#[derive(Debug, Clone, Serialize, Deserialize, CopyGetters, Getters)]
pub struct ClockModel {
  #[getset(get_copy = "pub")]
  clock_rate: f64,

  #[getset(get_copy = "pub")]
  intercept: f64,

  #[getset(get = "pub")]
  stats: ClockModelStats,
}

impl ClockModel {
  pub fn new(clock_set: &ClockSet) -> Result<Self, Report> {
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
      stats: ClockModelStats::Estimated(RegressionStats {
        chisq: clock_set.chisq(),
        r_val: clock_set.r_val(),
        hessian: clock_set.hessian(),
        cov: clock_set.cov(),
      }),
    })
  }

  pub fn with_fixed_rate(clock_set: &ClockSet, clock_rate: f64) -> ClockModel {
    ClockModel {
      clock_rate,
      intercept: clock_set.intercept(clock_rate),
      stats: ClockModelStats::Fixed,
    }
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

  pub fn date(&self, div: f64) -> f64 {
    (div - self.intercept()) / self.clock_rate()
  }

  pub fn div(&self, date: f64) -> f64 {
    date * self.clock_rate() + self.intercept()
  }

  pub fn clock_deviation(&self, date: f64, div: f64) -> f64 {
    self.div(date) - div
  }

  /// Time of root (most recent common ancestor)
  pub fn t_mrca(&self) -> f64 {
    self.date(0.0)
  }

  /// String showing line equation (for display)
  pub fn equation_str(&self) -> String {
    format!(
      "div = {:}t {:} {:}",
      float_to_significant_digits(self.clock_rate(), 3),
      if self.intercept() < 0.0 { "-" } else { "+" },
      float_to_significant_digits(self.intercept().abs(), 3)
    )
  }
}
