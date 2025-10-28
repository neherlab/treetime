use crate::distribution::distribution::Distribution;
use crate::distribution::distribution_function::DistributionFunction;
use crate::distribution::distribution_point::DistributionPoint;
use crate::distribution::distribution_range::DistributionRange;
use eyre::Report;
use ndarray::Array1;

/// Negates a distribution by reflecting it across the time axis: f(x) → f(-x).
pub fn distribution_negation(dist: &Distribution) -> Distribution {
  match dist {
    Distribution::Empty => Distribution::empty(),
    Distribution::Point(p) => negate_point(p),
    Distribution::Range(r) => negate_range(r),
    Distribution::Function(f) => negate_function(f).unwrap(),
  }
}

/// Negates a distribution in-place by reflecting it across the time axis: f(x) → f(-x).
pub fn distribution_negation_inplace(dist: &mut Distribution) {
  match dist {
    Distribution::Empty => {},
    Distribution::Point(p) => negate_point_inplace(p),
    Distribution::Range(r) => negate_range_inplace(r),
    Distribution::Function(f) => negate_function_inplace(f),
  }
}

fn negate_point(point: &DistributionPoint<f64>) -> Distribution {
  Distribution::point(-point.t(), point.amplitude())
}

fn negate_point_inplace(point: &mut DistributionPoint<f64>) {
  *point = DistributionPoint::new(-point.t(), point.amplitude());
}

fn negate_range(range: &DistributionRange<f64>) -> Distribution {
  Distribution::range((-range.end(), -range.start()), range.amplitude())
}

fn negate_range_inplace(range: &mut DistributionRange<f64>) {
  *range = DistributionRange::new((-range.end(), -range.start()), range.amplitude());
}

fn negate_function(func: &DistributionFunction<f64>) -> Result<Distribution, Report> {
  let t = func.t();
  let y = func.y();

  let negated_t = t.mapv(|x| -x);
  let mut t_vec: Vec<f64> = negated_t.to_vec();
  let mut y_vec: Vec<f64> = y.to_vec();

  t_vec.reverse();
  y_vec.reverse();

  Distribution::function(Array1::from_vec(t_vec), Array1::from_vec(y_vec))
}

fn negate_function_inplace(func: &mut DistributionFunction<f64>) {
  let t = func.t_mut();
  t.mapv_inplace(|x| -x);

  let t_len = t.len();
  for i in 0..t_len / 2 {
    t.swap(i, t_len - 1 - i);
  }

  let y = func.y_mut();
  let y_len = y.len();
  for i in 0..y_len / 2 {
    y.swap(i, y_len - 1 - i);
  }
}

#[cfg(test)]
mod tests {
  use super::*;
  use approx::assert_ulps_eq;

  #[test]
  fn test_negate_empty() {
    let dist = Distribution::empty();
    let negated = distribution_negation(&dist);
    assert!(matches!(negated, Distribution::Empty));
  }

  #[test]
  fn test_negate_point() {
    let dist = Distribution::point(2.0, 3.0);
    let negated = distribution_negation(&dist);
    match negated {
      Distribution::Point(p) => {
        assert_ulps_eq!(p.t(), -2.0);
        assert_ulps_eq!(p.amplitude(), 3.0);
      },
      _ => panic!("Expected Point variant"),
    }
  }

  #[test]
  fn test_negate_point_zero() {
    let dist = Distribution::point(0.0, 5.0);
    let negated = distribution_negation(&dist);
    match negated {
      Distribution::Point(p) => {
        assert_ulps_eq!(p.t(), 0.0);
        assert_ulps_eq!(p.amplitude(), 5.0);
      },
      _ => panic!("Expected Point variant"),
    }
  }

  #[test]
  fn test_negate_range() {
    let dist = Distribution::range((1.0, 4.0), 2.0);
    let negated = distribution_negation(&dist);
    match negated {
      Distribution::Range(r) => {
        assert_ulps_eq!(r.start(), -4.0);
        assert_ulps_eq!(r.end(), -1.0);
        assert_ulps_eq!(r.amplitude(), 2.0);
      },
      _ => panic!("Expected Range variant"),
    }
  }

  #[test]
  fn test_negate_range_symmetric() {
    let dist = Distribution::range((-3.0, 3.0), 1.0);
    let negated = distribution_negation(&dist);
    match negated {
      Distribution::Range(r) => {
        assert_ulps_eq!(r.start(), -3.0);
        assert_ulps_eq!(r.end(), 3.0);
        assert_ulps_eq!(r.amplitude(), 1.0);
      },
      _ => panic!("Expected Range variant"),
    }
  }

  #[test]
  fn test_negate_function() {
    let t = Array1::from_vec(vec![0.0, 1.0, 2.0]);
    let y = Array1::from_vec(vec![1.0, 2.0, 3.0]);
    let dist = Distribution::function(t, y).unwrap();

    let negated = distribution_negation(&dist);
    match negated {
      Distribution::Function(f) => {
        assert_eq!(f.t().len(), 3);
        assert_ulps_eq!(f.t()[0], -2.0);
        assert_ulps_eq!(f.t()[1], -1.0);
        assert_ulps_eq!(f.t()[2], 0.0);
        assert_ulps_eq!(f.y()[0], 3.0);
        assert_ulps_eq!(f.y()[1], 2.0);
        assert_ulps_eq!(f.y()[2], 1.0);
      },
      _ => panic!("Expected Function variant"),
    }
  }

  #[test]
  fn test_negate_inplace_point() {
    let mut dist = Distribution::point(2.0, 3.0);
    distribution_negation_inplace(&mut dist);
    match dist {
      Distribution::Point(p) => {
        assert_ulps_eq!(p.t(), -2.0);
        assert_ulps_eq!(p.amplitude(), 3.0);
      },
      _ => panic!("Expected Point variant"),
    }
  }

  #[test]
  fn test_negate_inplace_range() {
    let mut dist = Distribution::range((1.0, 4.0), 2.0);
    distribution_negation_inplace(&mut dist);
    match dist {
      Distribution::Range(r) => {
        assert_ulps_eq!(r.start(), -4.0);
        assert_ulps_eq!(r.end(), -1.0);
        assert_ulps_eq!(r.amplitude(), 2.0);
      },
      _ => panic!("Expected Range variant"),
    }
  }

  #[test]
  fn test_negate_inplace_function() {
    let t = Array1::from_vec(vec![0.0, 1.0, 2.0]);
    let y = Array1::from_vec(vec![1.0, 2.0, 3.0]);
    let mut dist = Distribution::function(t, y).unwrap();

    distribution_negation_inplace(&mut dist);
    match dist {
      Distribution::Function(f) => {
        assert_eq!(f.t().len(), 3);
        assert_ulps_eq!(f.t()[0], -2.0);
        assert_ulps_eq!(f.t()[1], -1.0);
        assert_ulps_eq!(f.t()[2], 0.0);
        assert_ulps_eq!(f.y()[0], 3.0);
        assert_ulps_eq!(f.y()[1], 2.0);
        assert_ulps_eq!(f.y()[2], 1.0);
      },
      _ => panic!("Expected Function variant"),
    }
  }
}
