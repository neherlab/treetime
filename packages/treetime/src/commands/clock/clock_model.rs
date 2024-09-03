use crate::commands::clock::clock_set::ClockSet;
use crate::make_error;
use crate::utils::float_fmt::float_to_significant_digits;
use eyre::Report;
use getset::{CopyGetters, Getters};
use ndarray::Array2;
use serde::{Deserialize, Serialize};

#[must_use]
#[derive(Debug, Default, Serialize, Deserialize, CopyGetters, Getters)]
pub struct ClockModel {
  #[getset(get_copy = "pub")]
  clock_rate: f64,

  #[getset(get_copy = "pub")]
  intercept: f64,

  #[getset(get_copy = "pub")]
  chisq: f64,

  #[getset(get_copy = "pub")]
  r_val: f64,

  #[getset(get = "pub")]
  hessian: Array2<f64>,

  #[getset(get = "pub")]
  cov: Array2<f64>,
}

impl ClockModel {
  pub fn new(clock_set: &ClockSet) -> Result<Self, Report> {
    let det = clock_set.determinant();

    if det <= 0.0 {
      return make_error!("No variation in sampling dates! Please specify your clock rate explicitly.");
    }

    Ok(Self {
      clock_rate: clock_set.clock_rate(det),
      intercept: clock_set.intercept(clock_set.clock_rate(det)),
      chisq: clock_set.chisq(),
      r_val: clock_set.r_val(),
      hessian: clock_set.hessian(),
      cov: clock_set.cov(),
    })
  }

  pub fn clock_model_fixed_rate(clock_set: &ClockSet, clock_rate: f64) -> ClockModel {
    let intercept = clock_set.intercept(clock_rate);
    let hessian = clock_set.hessian();
    ClockModel {
      clock_rate,
      intercept,
      hessian,
      chisq: 0.0,
      r_val: 0.0,
      cov: clock_set.cov(),
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
