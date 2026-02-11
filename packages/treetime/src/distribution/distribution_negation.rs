use crate::distribution::distribution::Distribution;
use crate::distribution::distribution_function::DistributionFunction;
use crate::distribution::distribution_point::DistributionPoint;
use crate::distribution::distribution_range::DistributionRange;
use crate::distribution::y_axis_policy::YAxisPolicy;
use eyre::Report;

/// Negate a distribution by reflecting it across the time axis: f(x) -> f(-x).
pub fn distribution_negation<Y: YAxisPolicy>(dist: &Distribution<Y>) -> Distribution<Y> {
  match dist {
    Distribution::Empty => Distribution::empty(),
    Distribution::Point(p) => negate_point(p),
    Distribution::Range(r) => negate_range(r),
    Distribution::Function(f) => negate_function(f).unwrap(),
    Distribution::Formula(_) => panic!("Negation not implemented for Formula distributions"),
  }
}

/// Negate a distribution in-place by reflecting it across the time axis: f(x) -> f(-x).
pub fn distribution_negation_inplace<Y: YAxisPolicy>(dist: &mut Distribution<Y>) {
  match dist {
    Distribution::Empty => {},
    Distribution::Point(p) => negate_point_inplace(p),
    Distribution::Range(r) => negate_range_inplace(r),
    Distribution::Function(f) => negate_function_inplace(f),
    Distribution::Formula(_) => panic!("Negation in place not implemented for Formula distributions"),
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
  result.negate_arg_inplace();
  Ok(Distribution::Function(result))
}

fn negate_function_inplace<Y: YAxisPolicy>(func: &mut DistributionFunction<f64, Y>) {
  func.negate_arg_inplace();
}
