use crate::Distribution;
use crate::distribution_core::function::DistributionFunction;
use crate::distribution_core::point::DistributionPoint;
use crate::distribution_core::range::DistributionRange;
use crate::policy::YAxisPolicy;
use eyre::Report;
use treetime_utils::make_error;

/// Negate a distribution by reflecting it across the time axis: f(x) -> f(-x).
pub fn distribution_negation<Y: YAxisPolicy>(dist: &Distribution<Y>) -> Result<Distribution<Y>, Report> {
  match dist {
    Distribution::Empty => Ok(Distribution::empty()),
    Distribution::Point(p) => Ok(negate_point(p)),
    Distribution::Range(r) => Ok(negate_range(r)),
    Distribution::Function(f) => negate_function(f),
    Distribution::Formula(_) => make_error!("Negation not implemented for Formula distributions"),
  }
}

/// Negate a distribution in-place by reflecting it across the time axis: f(x) -> f(-x).
pub fn distribution_negation_inplace<Y: YAxisPolicy>(dist: &mut Distribution<Y>) -> Result<(), Report> {
  match dist {
    Distribution::Empty => Ok(()),
    Distribution::Point(p) => {
      negate_point_inplace(p);
      Ok(())
    },
    Distribution::Range(r) => {
      negate_range_inplace(r);
      Ok(())
    },
    Distribution::Function(f) => negate_function_inplace(f),
    Distribution::Formula(_) => make_error!("Negation in place not implemented for Formula distributions"),
  }
}

fn negate_point<Y: YAxisPolicy>(point: &DistributionPoint<f64, Y>) -> Distribution<Y> {
  Distribution::point(-point.t(), point.amplitude())
}

fn negate_point_inplace<Y: YAxisPolicy>(point: &mut DistributionPoint<f64, Y>) {
  *point = DistributionPoint::new(-point.t(), point.amplitude());
}

fn negate_range<Y: YAxisPolicy>(range: &DistributionRange<f64, Y>) -> Distribution<Y> {
  Distribution::range((-range.end(), -range.start()), range.amplitude())
}

fn negate_range_inplace<Y: YAxisPolicy>(range: &mut DistributionRange<f64, Y>) {
  *range = DistributionRange::new((-range.end(), -range.start()), range.amplitude());
}

fn negate_function<Y: YAxisPolicy>(func: &DistributionFunction<f64, Y>) -> Result<Distribution<Y>, Report> {
  let mut result = func.clone();
  result.negate_arg_inplace()?;
  Ok(Distribution::Function(result))
}

fn negate_function_inplace<Y: YAxisPolicy>(func: &mut DistributionFunction<f64, Y>) -> Result<(), Report> {
  func.negate_arg_inplace()
}
