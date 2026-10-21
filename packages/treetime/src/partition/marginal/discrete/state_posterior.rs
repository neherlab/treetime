use crate::make_internal_report;
use eyre::Report;
use getset::{CopyGetters, Getters};
use ndarray::Array1;
use treetime_utils::array::ndarray::argmax_first;

const ENTROPY_PROBABILITY_OFFSET: f64 = 1e-12;

#[derive(Clone, Debug, PartialEq, CopyGetters, Getters)]
pub struct StatePosterior {
  #[getset(get_copy = "pub")]
  state: usize,
  #[getset(get = "pub")]
  profile: Array1<f64>,
  #[getset(get_copy = "pub")]
  entropy: f64,
}

impl StatePosterior {
  pub fn new(profile: Array1<f64>) -> Result<Self, Report> {
    let state = argmax_first(&profile.view())
      .ok_or_else(|| make_internal_report!("A state posterior needs a profile with at least one state"))?;
    let entropy = -profile
      .iter()
      .map(|&p| p * (p + ENTROPY_PROBABILITY_OFFSET).ln())
      .sum::<f64>();
    Ok(Self {
      state,
      profile,
      entropy,
    })
  }
}
