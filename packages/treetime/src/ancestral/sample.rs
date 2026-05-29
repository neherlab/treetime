use ndarray::ArrayView1;
use rand::Rng;
use serde::{Deserialize, Serialize};
use smart_default::SmartDefault;
use treetime_utils::array::ndarray::argmax_first;

#[derive(Clone, Copy, Debug, PartialEq, Eq, SmartDefault, Serialize, Deserialize)]
#[cfg_attr(feature = "clap", derive(clap::ValueEnum))]
#[cfg_attr(feature = "clap", value(rename_all = "kebab-case"))]
pub enum SampleMode {
  #[default]
  Argmax,
  Root,
  All,
}

pub fn sample_from_profile(profile: ArrayView1<f64>, rng: &mut impl Rng) -> usize {
  let cumsum: Vec<f64> = profile
    .iter()
    .scan(0.0, |acc, &x| {
      *acc += x;
      Some(*acc)
    })
    .collect();

  let total = cumsum.last().copied().unwrap_or(0.0);
  if total <= 0.0 {
    return 0;
  }

  let threshold = rng.r#gen::<f64>() * total;
  cumsum.iter().position(|&c| c >= threshold).unwrap_or(0)
}

pub fn resolve_profile(profile: ArrayView1<f64>, sample: bool, rng: &mut impl Rng) -> usize {
  if sample {
    sample_from_profile(profile, rng)
  } else {
    argmax_first(&profile).unwrap_or(0)
  }
}
