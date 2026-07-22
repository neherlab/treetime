use crate::make_error;
use crate::payload::clock_set::ClockSet;
use eyre::Report;
use getset::Getters;
use log::{debug, warn};
use ndarray::Array2;
use serde::{Deserialize, Serialize};
use treetime_utils::fmt::float::float_to_significant_digits;
use treetime_utils::io::json::{JsonPretty, json_write_str};

/// Shared interface for types that define a clock line: `div = rate * date + intercept`.
pub trait ClockLine {
  fn clock_rate(&self) -> f64;
  fn intercept(&self) -> f64;
  fn clock_deviation(&self, date: f64, div: f64) -> f64;
}

/// Regression statistics from clock model estimation
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RegressionStats {
  pub chisq: f64,
  pub r_val: f64,
  pub hessian: Array2<f64>,
  pub cov: Array2<f64>,
}

/// Raw regression result from root-to-tip analysis. Clock rate can be any sign.
///
/// Used by the estimation pipeline (clock filtering, root search diagnostics)
/// where negative rates are acceptable. Consumers that require positive rates
/// convert to `ClockModel` at their boundary.
#[must_use]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ClockRegression {
  clock_rate: f64,
  intercept: f64,
  chisq: f64,
  r_val: f64,
  hessian: Array2<f64>,
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

#[allow(clippy::same_name_method)]
impl ClockRegression {
  pub fn clock_rate(&self) -> f64 {
    self.clock_rate
  }

  pub fn intercept(&self) -> f64 {
    self.intercept
  }

  pub fn chisq(&self) -> f64 {
    self.chisq
  }

  pub fn r_val(&self) -> f64 {
    self.r_val
  }

  pub fn from_clock_set(clock_set: &ClockSet) -> Result<Self, Report> {
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

  pub fn hessian(&self) -> &Array2<f64> {
    &self.hessian
  }

  pub fn cov(&self) -> &Array2<f64> {
    &self.cov
  }
}

/// Clock model statistics - either estimated from data or fixed by user
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum ClockModelStats {
  /// Clock rate estimated from root-to-tip regression with full statistics
  Estimated(RegressionStats),
  /// Clock rate fixed by user input (no regression statistics available)
  Fixed,
}

/// Fitted clock model: a root-to-tip line `div = rate * date + intercept`.
///
/// The rate sign depends on how the model is constructed:
/// - `from_regression` / `with_fixed_rate` guarantee a positive rate. Time-scaled
///   analysis (timetree) requires this because `time = divergence / rate` is only
///   well-defined for a positive rate; a non-positive rate is rejected as an error.
/// - `from_regression_allow_negative` permits a non-positive rate (with a warning)
///   for the clock command, which only reports the root-to-tip regression and does
///   not perform time inference. This matches v0, which continues with negative
///   rates (e.g. under `--keep-root` or `--allow-negative-rate`).
#[must_use]
#[derive(Debug, Clone, Serialize, Deserialize, Getters)]
pub struct ClockModel {
  clock_rate: f64,

  intercept: f64,

  #[getset(get = "pub")]
  stats: ClockModelStats,
}

#[allow(clippy::same_name_method)]
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

  /// Builds a clock model from regression, permitting a non-positive rate.
  ///
  /// Emits a warning when the rate is non-positive but proceeds, unlike
  /// `from_regression` which errors. Used by the clock command, which reports the
  /// root-to-tip regression without performing time inference. A negative slope is
  /// a legitimate (if temporally uninformative) regression result; rejecting it
  /// would defeat `--allow-negative-rate` and `--keep-root`. See
  /// `kb/decisions/timetree-rejects-negative-clock-rate.md`.
  pub fn from_regression_allow_negative(regression: &ClockRegression) -> Self {
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

  pub fn with_fixed_rate(clock_set: &ClockSet, clock_rate: f64) -> Result<Self, Report> {
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

  pub fn date(&self, div: f64) -> f64 {
    (div - self.intercept()) / self.clock_rate()
  }

  pub fn div(&self, date: f64) -> f64 {
    date * self.clock_rate() + self.intercept()
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

  /// Create a clock model with specified rate and intercept (for testing)
  #[cfg(test)]
  pub fn for_testing(clock_rate: f64, intercept: f64) -> Self {
    Self {
      clock_rate,
      intercept,
      stats: ClockModelStats::Fixed,
    }
  }

  /// Create a clock model with specified rate, intercept, and stats (for testing)
  #[cfg(test)]
  pub fn for_testing_with_stats(clock_rate: f64, intercept: f64, stats: ClockModelStats) -> Self {
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
