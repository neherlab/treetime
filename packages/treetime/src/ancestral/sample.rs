use ndarray::ArrayView1;
use rand::{Rng, RngCore};
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use smart_default::SmartDefault;
use treetime_utils::array::ndarray::argmax_first;

#[derive(Clone, Copy, Debug, PartialEq, Eq, SmartDefault, Serialize, Deserialize, JsonSchema)]
#[serde(rename_all = "kebab-case")]
pub enum SampleMode {
  #[default]
  Argmax,
  Root,
  All,
}

impl SampleMode {
  pub(crate) fn samples_node(self, is_root: bool) -> bool {
    match self {
      SampleMode::Argmax => false,
      SampleMode::Root => is_root,
      SampleMode::All => true,
    }
  }
}

pub(crate) fn sample_from_profile(profile: ArrayView1<f64>, rng: &mut dyn RngCore) -> usize {
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

pub(crate) enum Resolve<'r> {
  Argmax,
  Sample(&'r mut dyn RngCore),
}

pub(crate) fn resolve_profile(profile: ArrayView1<f64>, resolve: &mut Resolve) -> usize {
  match resolve {
    Resolve::Argmax => argmax_first(&profile).unwrap_or(0),
    Resolve::Sample(rng) => sample_from_profile(profile, &mut **rng),
  }
}
