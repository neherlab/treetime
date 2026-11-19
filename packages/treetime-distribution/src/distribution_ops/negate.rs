use crate::Distribution;
use crate::distribution_core::function::DistributionFunction;
use crate::distribution_core::point::DistributionPoint;
use crate::distribution_core::range::DistributionRange;
use eyre::Report;

pub(crate) fn distribution_negation(dist: &Distribution) -> Result<Distribution, Report> {
  match dist {
    Distribution::Empty => Ok(Distribution::empty()),
    Distribution::Point(p) => Ok(negate_point(p)),
    Distribution::Range(r) => Ok(negate_range(r)),
    Distribution::Function(f) => negate_function(f),
  }
}

fn negate_point(point: &DistributionPoint) -> Distribution {
  Distribution::point(-point.t(), point.amplitude())
}

fn negate_range(range: &DistributionRange) -> Distribution {
  Distribution::range((-range.end(), -range.start()), range.amplitude())
}

fn negate_function(func: &DistributionFunction) -> Result<Distribution, Report> {
  let mut result = func.clone();
  result.negate_arg_inplace()?;
  Ok(Distribution::Function(result))
}
