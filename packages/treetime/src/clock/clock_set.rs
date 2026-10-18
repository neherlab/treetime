use auto_ops::impl_op_ex;
use deser::Serialize;
use getset::CopyGetters;
use ndarray::{Array2, array};

#[must_use]
#[derive(Debug, Default, Clone, PartialEq, CopyGetters, Serialize)]
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
  pub(crate) fn leaf_contribution(date: Option<f64>) -> Self {
    Self {
      t_sum: date.unwrap_or(0.0),
      tsq_sum: date.unwrap_or(0.0).powi(2),
      d_sum: 0.0,
      dt_sum: 0.0,
      dsq_sum: 0.0,
      norm: if date.is_some() { 1.0 } else { 0.0 },
    }
  }

  pub(crate) fn outlier_contribution() -> Self {
    Self {
      t_sum: 0.0,
      tsq_sum: 0.0,
      d_sum: 0.0,
      dt_sum: 0.0,
      dsq_sum: 0.0,
      norm: 0.0,
    }
  }

  /// Sums with every divergence in the set lengthened by `branch_length`, as seen from the parent end of a branch:
  /// `d -> d + l`, so `sum d += l n`, `sum t d += l sum t`, `sum d^2 += 2 l sum d + l^2 n`; the date sums stay.
  pub(crate) fn shift(&self, branch_length: f64) -> Self {
    let l = branch_length;
    Self {
      d_sum: self.d_sum + l * self.norm,
      dt_sum: self.dt_sum + l * self.t_sum,
      dsq_sum: self.dsq_sum + 2.0 * l * self.d_sum + l * l * self.norm,
      ..self.clone()
    }
  }

  pub(crate) fn determinant(&self) -> f64 {
    self.tsq_sum() * self.norm() - self.t_sum().powi(2)
  }

  pub(crate) fn clock_rate(&self, det: f64) -> f64 {
    (self.dt_sum() * self.norm() - self.t_sum() * self.d_sum()) / det
  }

  pub(crate) fn intercept(&self, rate: f64) -> f64 {
    (self.d_sum() - self.t_sum() * rate) / self.norm()
  }

  pub(crate) fn hessian(&self) -> Array2<f64> {
    array![[self.tsq_sum(), self.t_sum()], [self.t_sum(), self.norm()]]
  }

  pub(crate) fn cov(&self) -> Array2<f64> {
    let det_inv = 1.0 / self.determinant();
    array![
      [self.norm() * det_inv, -self.t_sum() * det_inv],
      [-self.t_sum() * det_inv, self.tsq_sum() * det_inv],
    ]
  }

  pub(crate) fn chisq(&self) -> f64 {
    let det = self.determinant();
    0.5
      * (self.dsq_sum() * self.norm()
        - self.d_sum().powi(2)
        - (self.dt_sum() * self.norm() - self.d_sum() * self.t_sum()).powi(2) / det)
      / self.norm()
  }

  pub(crate) fn chisq_fixed_rate(&self, rate: f64) -> f64 {
    0.5
      * (self.dsq_sum() - 2.0 * rate * self.dt_sum() + rate.powi(2) * self.tsq_sum()
        - (self.d_sum() - rate * self.t_sum()).powi(2) / self.norm())
      / self.norm()
  }

  pub(crate) fn r_val(&self) -> f64 {
    (self.dt_sum * self.norm() - self.t_sum * self.d_sum)
      / ((self.dsq_sum * self.norm() - self.d_sum.powi(2)) * (self.tsq_sum * self.norm() - self.t_sum.powi(2))).sqrt()
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
