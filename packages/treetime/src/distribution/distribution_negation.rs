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

#[cfg(test)]
mod tests {
  use super::{distribution_negation, distribution_negation_inplace};
  use crate::distribution::distribution::DistributionPlain as Distribution;
  use ndarray::array;

  #[test]
  fn test_negate_empty() {
    let dist: Distribution = Distribution::empty();
    let actual: Distribution = distribution_negation(&dist);
    let expected: Distribution = Distribution::empty();
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_negate_point() {
    let dist: Distribution = Distribution::point(2.0, 3.0);
    let actual: Distribution = distribution_negation(&dist);
    let expected: Distribution = Distribution::point(-2.0, 3.0);
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_negate_point_zero() {
    let dist: Distribution = Distribution::point(0.0, 5.0);
    let actual: Distribution = distribution_negation(&dist);
    let expected: Distribution = Distribution::point(0.0, 5.0);
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_negate_range() {
    let dist: Distribution = Distribution::range((1.0, 4.0), 2.0);
    let actual: Distribution = distribution_negation(&dist);
    let expected: Distribution = Distribution::range((-4.0, -1.0), 2.0);
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_negate_range_symmetric() {
    let dist: Distribution = Distribution::range((-3.0, 3.0), 1.0);
    let actual: Distribution = distribution_negation(&dist);
    let expected: Distribution = Distribution::range((-3.0, 3.0), 1.0);
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_negate_function() {
    let t = array![0.0, 1.0, 2.0];
    let y = array![1.0, 2.0, 3.0];
    let dist: Distribution = Distribution::function(t, y).unwrap();

    let actual: Distribution = distribution_negation(&dist);

    let expected_t = array![-2.0, -1.0, 0.0];
    let expected_y = array![3.0, 2.0, 1.0];
    let expected: Distribution = Distribution::function(expected_t, expected_y).unwrap();
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_negate_inplace_point() {
    let mut actual: Distribution = Distribution::point(2.0, 3.0);
    distribution_negation_inplace(&mut actual);
    let expected: Distribution = Distribution::point(-2.0, 3.0);
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_negate_inplace_range() {
    let mut actual: Distribution = Distribution::range((1.0, 4.0), 2.0);
    distribution_negation_inplace(&mut actual);
    let expected: Distribution = Distribution::range((-4.0, -1.0), 2.0);
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_negate_inplace_function() {
    let t = array![0.0, 1.0, 2.0];
    let y = array![1.0, 2.0, 3.0];
    let mut actual: Distribution = Distribution::function(t, y).unwrap();

    distribution_negation_inplace(&mut actual);

    let expected_t = array![-2.0, -1.0, 0.0];
    let expected_y = array![3.0, 2.0, 1.0];
    let expected: Distribution = Distribution::function(expected_t, expected_y).unwrap();
    assert_eq!(expected, actual);
  }
}
