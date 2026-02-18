pub(crate) mod distribution_core;
pub(crate) mod distribution_ops;
pub(crate) mod distribution_scaled;
pub(crate) mod policy;

mod __tests__;

pub use distribution_core::distribution::{
  Distribution, DistributionNegLog, DistributionPlain, TIME_EPSILON, TIME_LIMIT,
};
pub use distribution_core::formula::DistributionFormula;
pub use distribution_core::function::DistributionFunction;
pub use distribution_core::point::DistributionPoint;
pub use distribution_core::range::DistributionRange;
pub use distribution_ops::convolve::distribution_convolution;
pub use distribution_ops::divide::distribution_division;
pub use distribution_ops::map::distribution_map;
pub use distribution_ops::multiply::distribution_multiplication;
pub use distribution_ops::negate::{distribution_negation, distribution_negation_inplace};
pub use distribution_ops::scalar_multiply::distribution_scalar_multiplication;
pub use distribution_ops::subtract::distribution_subtraction;
pub use distribution_ops::time_bounds::{
  distribution_time_bounds_contains, distribution_time_bounds_intersection, distribution_time_bounds_overlaps,
  distribution_time_bounds_union,
};
pub use distribution_scaled::convolve::scaled_distribution_convolution;
pub use distribution_scaled::divide::scaled_distribution_division;
pub use distribution_scaled::multiply::{scaled_distribution_multiplication, scaled_distribution_multiply_many};
pub use distribution_scaled::scaled::ScaledDistribution;
pub use policy::{NegLog, Plain, PolicyMarker, SupportsConvolution, SupportsSubtraction, YAxisPolicy};
