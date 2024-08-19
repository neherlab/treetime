use crate::make_error;
use auto_ops::impl_op_ex;
use eyre::Report;
use getset::{CopyGetters, Getters};
use ndarray::Array2;
use serde::{Deserialize, Serialize};

#[must_use]
#[derive(Debug, Default, Serialize, Deserialize, CopyGetters, Getters)]
pub struct ClockModel {
  #[getset(get_copy = "pub")]
  rate: f64,

  #[getset(get_copy = "pub")]
  intercept: f64,

  #[getset(get_copy = "pub")]
  chisq: f64,

  #[getset(get = "pub")]
  hessian: Array2<f64>,
}

#[must_use]
#[derive(Debug, Default, Clone, Serialize, Deserialize, CopyGetters)]
#[getset(get_copy = "pub")]
pub struct ClockSet {
  t_sum: f64,
  tsq_sum: f64,
  d_sum: f64,
  dsq_sum: f64,
  dt_sum: f64,
  norm: f64,
}

impl ClockSet {
  pub fn leaf_contribution(tip_value: f64) -> Self {
    Self {
      t_sum: tip_value,
      tsq_sum: tip_value.powi(2),
      d_sum: 0.0,
      dt_sum: 0.0,
      dsq_sum: 0.0,
      norm: 1.0,
    }
  }

  pub fn propagate_averages(&self, branch_value: f64, branch_variance: f64) -> Self {
    let denom = 1.0 / (1.0 + branch_variance * self.norm);

    // Eq 11 in Neher 2018 -- contribution of children doesn't change
    let t_sum = self.t_sum * denom;

    // Eq. 13
    let tsq_sum = self.tsq_sum - branch_variance * self.t_sum.powi(2) * denom;

    // Eq. 12 -- add branch_value norm times
    let d_sum = (self.d_sum + branch_value * self.norm) * denom;

    // Eq. 14
    let dt_sum = self.dt_sum + branch_value * self.t_sum
      - branch_variance * self.t_sum * (self.d_sum + self.norm * branch_value) * denom;

    // Eq. A.2
    let dsq_sum = self.dsq_sum + 2.0 * branch_value * self.d_sum + branch_value.powi(2) * self.norm
      - branch_variance
        * (self.d_sum.powi(2) + 2.0 * branch_value * self.d_sum * self.norm + branch_value.powi(2) * self.norm.powi(2))
        * denom;

    let norm = self.norm * denom;

    Self {
      t_sum,
      tsq_sum,
      d_sum,
      dsq_sum,
      dt_sum,
      norm,
    }
  }

  pub fn clock_model(&self) -> Result<ClockModel, Report> {
    let det = self.determinant();

    if det <= 0.0 {
      return make_error!("No variation in sampling dates! Please specify your clock rate explicitly.");
    }

    let rate = self.rate(det);
    let intercept = self.intercept(rate);
    let hessian = self.hessian();
    let chisq = self.chisq(det);

    Ok(ClockModel {
      rate,
      intercept,
      chisq,
      hessian,
    })
  }

  pub fn clock_model_fixed_rate(&self, rate: f64) -> ClockModel {
    let intercept = self.intercept(rate);
    let hessian = self.hessian();

    ClockModel {
      rate,
      intercept,
      hessian,
      chisq: 0.0,
    }
  }

  fn determinant(&self) -> f64 {
    self.tsq_sum() * self.norm() - self.t_sum().powi(2)
  }

  fn rate(&self, det: f64) -> f64 {
    (self.dt_sum() * self.norm() - self.t_sum() * self.d_sum()) / det
  }

  fn intercept(&self, rate: f64) -> f64 {
    (self.d_sum() - self.t_sum() * rate) / self.norm()
  }

  fn hessian(&self) -> Array2<f64> {
    Array2::from_shape_vec((2, 2), vec![self.tsq_sum(), self.t_sum(), self.t_sum(), self.norm()]).unwrap()
  }

  fn chisq(&self, det: f64) -> f64 {
    0.5
      * (self.dsq_sum() * self.norm()
        - self.d_sum().powi(2)
        - (self.dt_sum() * self.norm() - self.d_sum() * self.t_sum()).powi(2) / det)
      / self.norm()
  }
}

impl_op_ex!(+ |a: &ClockSet, b: &ClockSet| -> ClockSet {
    ClockSet {
        t_sum: a.t_sum + b.t_sum,
        tsq_sum: a.tsq_sum + b.tsq_sum,
        d_sum: a.d_sum + b.d_sum,
        dsq_sum: a.dsq_sum + b.dsq_sum,
        dt_sum: a.dt_sum + b.dt_sum,
        norm: a.norm + b.norm,
    }
});

impl_op_ex!(-|a: &ClockSet, b: &ClockSet| -> ClockSet {
  ClockSet {
    t_sum: a.t_sum - b.t_sum,
    tsq_sum: a.tsq_sum - b.tsq_sum,
    d_sum: a.d_sum - b.d_sum,
    dsq_sum: a.dsq_sum - b.dsq_sum,
    dt_sum: a.dt_sum - b.dt_sum,
    norm: a.norm - b.norm,
  }
});

impl_op_ex!(+= |a: &mut ClockSet, b: &ClockSet| {
    a.t_sum += b.t_sum;
    a.tsq_sum += b.tsq_sum;
    a.d_sum += b.d_sum;
    a.dsq_sum += b.dsq_sum;
    a.dt_sum += b.dt_sum;
    a.norm += b.norm;
});

impl_op_ex!(-= |a: &mut ClockSet, b: &ClockSet| {
    a.t_sum -= b.t_sum;
    a.tsq_sum -= b.tsq_sum;
    a.d_sum -= b.d_sum;
    a.dsq_sum -= b.dsq_sum;
    a.dt_sum -= b.dt_sum;
    a.norm -= b.norm;
});
