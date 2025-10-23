use crate::distribution::distribution::Distribution;
use crate::distribution::distribution_convolution::distribution_convolution;

pub fn convolve_safe(a: &Distribution, b: &Distribution) -> Distribution {
  distribution_convolution(a, b).unwrap_or_else(|_| Distribution::empty())
}

pub fn multiply_distributions(existing: Option<Distribution>, new: Distribution) -> Distribution {
  if let Some(existing) = existing {
    convolve_safe(&existing, &new)
  } else {
    new
  }
}
