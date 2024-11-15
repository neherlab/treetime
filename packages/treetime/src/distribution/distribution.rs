use crate::distribution::distribution_function::DistributionFunction;
use crate::distribution::distribution_point::DistributionPoint;
use crate::distribution::distribution_range::DistributionRange;
use serde::{Deserialize, Serialize};

#[allow(variant_size_differences)]
#[derive(Default, Serialize, Deserialize)]
pub enum Distribution {
  #[default]
  Empty,
  Point(DistributionPoint<f64>),
  Range(DistributionRange<f64>),

  #[serde(skip)]
  Function(DistributionFunction<f64>),
}
