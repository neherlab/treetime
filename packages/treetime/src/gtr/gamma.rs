use crate::make_error;
use eyre::Report;
use getset::{CopyGetters, Getters};
use ndarray::Array1;
use statrs::distribution::{ContinuousCDF, Gamma};
use statrs::function::gamma::checked_gamma_lr;
use treetime_utils::cast::ToF64Exact;

pub const GAMMA_SHAPE_MIN: f64 = 0.02;

pub const GAMMA_SHAPE_MAX: f64 = 100.0;

pub const INITIAL_GAMMA_SHAPE: f64 = 1.0;

pub const NO_RATE_VARIATION: usize = 1;

#[derive(Clone, Debug, PartialEq, Getters, CopyGetters)]
pub struct GammaRates {
  #[getset(get_copy = "pub")]
  shape: Option<f64>,
  #[getset(get = "pub")]
  rates: Array1<f64>,
  #[getset(get = "pub")]
  weights: Array1<f64>,
}

impl GammaRates {
  pub fn new(categories: usize, shape: f64) -> Result<Self, Report> {
    if categories == 0 {
      return make_error!("Number of gamma rate categories must be at least 1, got 0");
    }
    if categories == 1 {
      return Ok(Self::default());
    }
    if !(GAMMA_SHAPE_MIN..=GAMMA_SHAPE_MAX).contains(&shape) {
      return make_error!("Gamma shape must lie between {GAMMA_SHAPE_MIN} and {GAMMA_SHAPE_MAX}, got {shape}");
    }
    let n = categories.to_f64_exact();
    let distribution = Gamma::new(shape, shape)?;
    let bounds = (0..=categories)
      .map(|k| match k {
        0 => Ok(0.0),
        k if k == categories => Ok(1.0),
        k => upper_mass(shape, distribution.inverse_cdf(k.to_f64_exact() / n)),
      })
      .collect::<Result<Vec<f64>, Report>>()?;
    let rates = bounds
      .windows(2)
      .map(|pair| n * (pair[1] - pair[0]))
      .collect::<Array1<f64>>();
    Ok(Self {
      shape: Some(shape),
      rates,
      weights: Array1::from_elem(categories, 1.0 / n),
    })
  }

  pub fn categories(&self) -> usize {
    self.rates.len()
  }
}

impl Default for GammaRates {
  fn default() -> Self {
    Self {
      shape: None,
      rates: Array1::ones(1),
      weights: Array1::ones(1),
    }
  }
}

fn upper_mass(shape: f64, cut: f64) -> Result<f64, Report> {
  if cut <= 0.0 {
    return Ok(0.0);
  }
  Ok(checked_gamma_lr(shape + 1.0, shape * cut)?)
}
