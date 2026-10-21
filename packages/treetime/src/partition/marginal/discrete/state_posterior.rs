use crate::make_internal_report;
use crate::partition::marginal::posterior_states::most_likely_state;
use eyre::Report;
use getset::{CopyGetters, Getters};
use ndarray::Array1;
use ndarray_stats::EntropyExt;

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
    let state = most_likely_state(&profile)
      .ok_or_else(|| make_internal_report!("A state posterior needs a profile with at least one state"))?;
    let entropy = profile
      .entropy()
      .map_err(|err| make_internal_report!("The entropy of a state posterior failed: {err}"))?;
    Ok(Self {
      state,
      profile,
      entropy: without_negative_zero(entropy),
    })
  }
}

fn without_negative_zero(value: f64) -> f64 {
  value + 0.0
}
