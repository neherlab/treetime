use itertools::izip;
use ndarray::Array1;
use std::collections::BTreeMap;
use treetime::partition::storage::discrete::DiscreteStates;

pub(crate) fn build_confidence_map(states: &DiscreteStates, profile: &Array1<f64>) -> BTreeMap<String, f64> {
  izip!(states.iter(), profile)
    .filter(|&(_, &probability)| probability > 0.001)
    .map(|(state, &probability)| (state.to_owned(), probability))
    .collect()
}

pub(crate) fn compute_entropy(profile: &Array1<f64>) -> f64 {
  const TINY: f64 = 1e-12;
  -profile.iter().map(|&p| p * (p + TINY).ln()).sum::<f64>()
}
