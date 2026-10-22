use crate::test_utils::AnalyticBranch;
use statrs::function::factorial::ln_factorial;
use std::sync::Arc;
use treetime_utils::cast::ToF64Exact;

#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct PoissonCherry {
  pub counts: (u64, u64),
  pub lambda: f64,
  pub tc: Option<f64>,
}

impl PoissonCherry {
  pub(crate) fn log_z(&self) -> f64 {
    let (k1, k2) = self.counts;
    let k = k1 + k2;
    let multinomial = ln_factorial(k) - ln_factorial(k1) - ln_factorial(k2);
    let k = k.to_f64_exact();
    match self.tc {
      None => multinomial - self.lambda.ln() - (k + 1.0) * 2.0_f64.ln(),
      Some(tc) => multinomial + k * self.lambda.ln() - tc.ln() - (k + 1.0) * self.root_age_rate().ln(),
    }
  }

  pub(crate) fn root_age(&self) -> GammaPosterior {
    let (k1, k2) = self.counts;
    GammaPosterior {
      shape: (k1 + k2).to_f64_exact() + 1.0,
      rate: self.root_age_rate(),
    }
  }

  pub(crate) fn log_z_slope(&self) -> f64 {
    let (k1, k2) = self.counts;
    let k = (k1 + k2).to_f64_exact();
    k - (k + 1.0) * 2.0 * self.lambda / self.root_age_rate()
  }

  pub(crate) fn branch(&self, count: u64) -> AnalyticBranch {
    let lambda = self.lambda;
    let log_count_factorial = ln_factorial(count);
    let count = count.to_f64_exact();
    let neg_log = move |tau: f64| {
      let mean = lambda * tau;
      if count > 0.0 {
        mean - count * mean.ln() + log_count_factorial
      } else {
        mean
      }
    };
    AnalyticBranch {
      neg_log: Arc::new(neg_log),
      monotone_after: count / lambda,
    }
  }

  pub(crate) fn root_neg_log_prior(&self, tip_date: f64, time: f64) -> f64 {
    self.tc.map_or(0.0, |tc| tc.ln() + (tip_date - time) / tc)
  }

  fn root_age_rate(&self) -> f64 {
    2.0 * self.lambda + self.tc.map_or(0.0, |tc| 1.0 / tc)
  }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct GammaPosterior {
  pub shape: f64,
  pub rate: f64,
}

impl GammaPosterior {
  pub(crate) fn mean(&self) -> f64 {
    self.shape / self.rate
  }

  pub(crate) fn variance(&self) -> f64 {
    self.shape / self.rate.powi(2)
  }
}
