use crate::Distribution;
use crate::distribution_core::function::DistributionFunction;
use crate::distribution_core::point::DistributionPoint;
use crate::distribution_core::range::DistributionRange;
use crate::policy::YAxisPolicy;
use eyre::Report;
use treetime_utils::make_error;

pub(crate) fn distribution_negation<Y: YAxisPolicy>(dist: &Distribution<Y>) -> Result<Distribution<Y>, Report> {
  match dist {
    Distribution::Empty => Ok(Distribution::empty()),
    Distribution::Point(p) => Ok(negate_point(p)),
    Distribution::Range(r) => Ok(negate_range(r)),
    Distribution::Function(f) => negate_function(f),
    Distribution::Formula(_) => make_error!("Negation not implemented for Formula distributions"),
  }
}

fn negate_point<Y: YAxisPolicy>(point: &DistributionPoint<f64, Y>) -> Distribution<Y> {
  Distribution::point(-point.t(), point.amplitude())
}

fn negate_range<Y: YAxisPolicy>(range: &DistributionRange<f64, Y>) -> Distribution<Y> {
  Distribution::range((-range.end(), -range.start()), range.amplitude())
}

fn negate_function<Y: YAxisPolicy>(func: &DistributionFunction<f64, Y>) -> Result<Distribution<Y>, Report> {
  let mut result = func.clone();
  result.negate_arg_inplace()?;
  Ok(Distribution::Function(result))
}
