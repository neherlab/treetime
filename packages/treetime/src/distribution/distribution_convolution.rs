use crate::distribution::distribution::Distribution;
use crate::distribution::distribution_function::DistributionFunction;
use crate::distribution::distribution_point::DistributionPoint;
use crate::distribution::distribution_range::DistributionRange;
use approx::ulps_eq;
use eyre::Report;
use itertools::{Itertools, izip};
use ndarray::{Array1, array};
use ordered_float::OrderedFloat;

pub fn distribution_convolution(a: &Distribution, b: &Distribution) -> Result<Distribution, Report> {
  match (a, b) {
    (Distribution::Empty, _) | (_, Distribution::Empty) => {
      Ok(Distribution::Empty) //
    },
    (Distribution::Point(a), Distribution::Point(b)) => {
      Ok(convolution_point_point(a, b)) //
    },
    (Distribution::Point(a), Distribution::Range(b)) | (Distribution::Range(b), Distribution::Point(a)) => {
      Ok(convolution_point_range(a, b)) //
    },
    (Distribution::Range(a), Distribution::Range(b)) => {
      convolution_range_range(a, b) //
    },
    (Distribution::Point(a), Distribution::Function(b)) | (Distribution::Function(b), Distribution::Point(a)) => {
      Ok(convolution_point_function(a, b)?) //
    },
    (Distribution::Range(a), Distribution::Function(b)) | (Distribution::Function(b), Distribution::Range(a)) => {
      Ok(convolution_range_function(a, b)) //
    },
    (Distribution::Function(a), Distribution::Function(b)) => {
      Ok(convolution_function_function(a, b)) //
    },
  }
}

fn convolution_point_point(a: &DistributionPoint<f64>, b: &DistributionPoint<f64>) -> Distribution {
  let x = a.t() + b.t();
  let y = a.amplitude() * b.amplitude();
  Distribution::point(x, y)
}

fn convolution_range_range(a: &DistributionRange<f64>, b: &DistributionRange<f64>) -> Result<Distribution, Report> {
  let start = a.start() + b.start();
  let end = a.end() + b.end();

  let peak_start = f64::max(a.start() + b.start(), a.end() + b.start());
  let peak_end = f64::min(a.end() + b.end(), a.start() + b.end());

  let peak_amplitude = a.amplitude() * b.amplitude();

  if ulps_eq!(&peak_start, &peak_end, max_ulps = 10) {
    let x = array![start, peak_start, end];
    let y = array![0.0, peak_amplitude, 0.0];
    Distribution::function(x, y)
  } else {
    let x = array![start, peak_start, peak_end, end];
    let y = array![0.0, peak_amplitude, peak_amplitude, 0.0];
    Distribution::function(x, y)
  }
}

fn convolution_point_range(p: &DistributionPoint<f64>, r: &DistributionRange<f64>) -> Distribution {
  let begin = r.start() + p.t();
  let end = r.end() + p.t();
  let amplitude = p.amplitude() * r.amplitude();
  Distribution::range((begin, end), amplitude)
}

fn convolution_point_function(
  p: &DistributionPoint<f64>,
  f: &DistributionFunction<f64>,
) -> Result<Distribution, Report> {
  let t = f.t().map(|t| t + p.t());
  let y = f.y().map(|y| y * p.amplitude());
  Distribution::function(t, y)
}

fn convolution_range_function(r: &DistributionRange<f64>, f: &DistributionFunction<f64>) -> Distribution {
  let t_out = f.t().clone();
  let mut y_out = Array1::zeros(f.y().len());

  for (i, &ti) in f.t().iter().enumerate() {
    let mask = f.t().mapv(|x| (x >= ti - r.end()) && (x <= ti - r.start()));
    let filtered_y = f.y() * &mask.mapv(|x| if x { 1.0 } else { 0.0 });
    y_out[i] = r.amplitude() * filtered_y.sum();
  }

  Distribution::function(t_out, y_out).unwrap()
}

/// Discrete convolution of two functions using functional ndarray operations
/// (a * b)(t) = sum_s a(s) * b(t - s)
fn convolution_function_function(a: &DistributionFunction<f64>, b: &DistributionFunction<f64>) -> Distribution {
  // Pair times with corresponding amplitudes for both functions
  let a = izip!(a.t().iter(), a.y().iter());
  let b = izip!(b.t().iter(), b.y().iter());

  let (times, amplitudes): (Vec<f64>, Vec<f64>) =
    // Generate all combinations of (time_a, amplitude_a) × (time_b, amplitude_b)
    a.cartesian_product(b)
    // Apply convolution: t = t_a + t_b, y = y_a * y_b
    .map(|((&ta, &ya), (&tb, &yb))| (OrderedFloat(ta + tb), ya * yb))
    // Group by output time: algorithmically builds HashMap<OrderedFloat, Vec<f64>> where
    // each unique time becomes a key, and all amplitudes mapping to that time are collected into a vector
    .into_group_map()
    .into_iter()
    // Sum amplitudes for each unique output time
    .map(|(time, amplitudes)| (time.into_inner(), amplitudes.into_iter().sum::<f64>()))
    // Sort results by time for proper distribution ordering
    .sorted_by_key(|(time, _)| OrderedFloat(*time))
    // Unzip into separate vectors for times and amplitudes
    .unzip();

  let t_out = Array1::from(times);
  let y_out = Array1::from(amplitudes);

  Distribution::function(t_out, y_out).unwrap()
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::pretty_assert_ulps_eq;
  use pretty_assertions::assert_eq;

  #[test]
  fn test_convolution_empty() {
    let a = Distribution::empty();
    let b = Distribution::function(array![], array![]).unwrap();
    let actual = distribution_convolution(&a, &b).unwrap();
    let expected = Distribution::empty();
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_convolution_point_point() {
    let a = Distribution::point(2.0, 3.0);
    let b = Distribution::point(5.0, 4.0);
    let actual = distribution_convolution(&a, &b).unwrap();
    let expected = Distribution::point(7.0, 12.0);
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_convolution_range_range_same_length() {
    let a = Distribution::range((2.0, 4.0), 3.0);
    let b = Distribution::range((6.0, 8.0), 2.0);
    let actual = distribution_convolution(&a, &b).unwrap();
    let expected = {
      let x = array![8.0, 10.0, 12.0];
      let y = array![0.0, 6.0, 0.0];
      Distribution::function(x, y).unwrap()
    };
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_convolution_range_range_different_length() {
    let a = Distribution::range((2.0, 4.0), 3.0);
    let b = Distribution::range((6.0, 9.0), 2.0);
    let actual = distribution_convolution(&a, &b).unwrap();
    let expected = {
      let x = array![8.0, 10.0, 11.0, 13.0];
      let y = array![0.0, 6.0, 6.0, 0.0];
      Distribution::function(x, y).unwrap()
    };
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_convolution_point_function() {
    let p = Distribution::point(3.0, 2.0);

    let x = array![0.0, 1.0, 2.0, 3.0, 4.0];
    let y = array![1.0, 2.0, 3.0, 4.0, 5.0];
    let f = Distribution::function(x, y).unwrap();
    let actual = distribution_convolution(&p, &f).unwrap();

    let x = array![3.0, 4.0, 5.0, 6.0, 7.0];
    let y = array![2.0, 4.0, 6.0, 8.0, 10.0];
    let expected = Distribution::function(x, y).unwrap();

    assert_eq!(expected, actual);
  }

  #[test]
  fn test_convolution_range_function() {
    let r = Distribution::range((1.0, 3.0), 2.0);

    let x = array![0.0, 1.0, 2.0, 3.0, 4.0, 5.0];
    let y = array![0.0, 1.0, 0.0, 2.0, 1.0, 0.0];
    let f = Distribution::function(x, y).unwrap();
    let actual = distribution_convolution(&r, &f).unwrap();

    let x = array![0.0, 1.0, 2.0, 3.0, 4.0, 5.0];
    let y = array![0.0, 0.0, 2.0, 2.0, 6.0, 6.0];
    let expected = Distribution::function(x, y).unwrap();

    assert_eq!(expected, actual);
  }

  #[test]
  fn test_convolution_point_range() {
    let p = Distribution::point(3.0, 2.0);
    let r = Distribution::range((1.0, 4.0), 1.5);
    let actual = distribution_convolution(&p, &r).unwrap();
    let expected = Distribution::range((4.0, 7.0), 3.0);
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_convolution_range_point() {
    let r = Distribution::range((1.0, 4.0), 1.5);
    let p = Distribution::point(3.0, 2.0);
    let actual = distribution_convolution(&r, &p).unwrap();
    let expected = Distribution::range((4.0, 7.0), 3.0);
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_convolution_function_function() {
    // Test 1: Two-point function convolving with three-point function
    // f1: triangle function with peaks at t=1.0, t=3.0
    let x1 = array![1.0, 3.0];
    let y1 = array![4.0, 7.0];
    let f1 = Distribution::function(x1, y1).unwrap();

    // f2: three-point function with values at t=0.0, t=2.0, t=4.0
    let x2 = array![0.0, 2.0, 4.0];
    let y2 = array![1.0, 8.0, 2.0];
    let f2 = Distribution::function(x2, y2).unwrap();

    let actual = distribution_convolution(&f1, &f2).unwrap();

    // Manual calculation of discrete convolution:
    // f1 has points: (1.0, 4.0), (3.0, 7.0)
    // f2 has points: (0.0, 1.0), (2.0, 8.0), (4.0, 2.0)
    //
    // Convolution produces peaks at all pairwise sums:
    // 1.0+0.0=1.0: 4.0*1.0=4.0
    // 1.0+2.0=3.0: 4.0*8.0=32.0
    // 1.0+4.0=5.0: 4.0*2.0=8.0
    // 3.0+0.0=3.0: 7.0*1.0=7.0
    // 3.0+2.0=5.0: 7.0*8.0=56.0
    // 3.0+4.0=7.0: 7.0*2.0=14.0
    //
    // Combined at unique time points:
    // t=1.0: 4.0
    // t=3.0: 32.0+7.0=39.0
    // t=5.0: 8.0+56.0=64.0
    // t=7.0: 14.0

    match actual {
      Distribution::Function(f) => {
        let expected_x = array![1.0, 3.0, 5.0, 7.0];
        let expected_y = array![4.0, 39.0, 64.0, 14.0];
        pretty_assert_ulps_eq!(f.t(), &expected_x, max_ulps = 10);
        pretty_assert_ulps_eq!(f.y(), &expected_y, max_ulps = 10);
      },
      _ => panic!("Expected Distribution::Function, got: {actual:#?}"),
    }
  }

  #[test]
  fn test_convolution_function_function_symmetric() {
    // Test 2: Symmetric convolution with round numbers
    // Both functions have same structure but different amplitudes
    let x1 = array![0.0, 2.0, 4.0];
    let y1 = array![1.0, 3.0, 1.0];
    let f1 = Distribution::function(x1, y1).unwrap();

    let x2 = array![0.0, 2.0, 4.0];
    let y2 = array![2.0, 4.0, 2.0];
    let f2 = Distribution::function(x2, y2).unwrap();

    let actual = distribution_convolution(&f1, &f2).unwrap();

    // Manual calculation:
    // f1: (0.0, 1.0), (2.0, 3.0), (4.0, 1.0)
    // f2: (0.0, 2.0), (2.0, 4.0), (4.0, 2.0)
    //
    // Convolution:
    // 0.0+0.0=0.0: 1.0*2.0=2.0
    // 0.0+2.0=2.0: 1.0*4.0=4.0
    // 0.0+4.0=4.0: 1.0*2.0=2.0
    // 2.0+0.0=2.0: 3.0*2.0=6.0
    // 2.0+2.0=4.0: 3.0*4.0=12.0
    // 2.0+4.0=6.0: 3.0*2.0=6.0
    // 4.0+0.0=4.0: 1.0*2.0=2.0
    // 4.0+2.0=6.0: 1.0*4.0=4.0
    // 4.0+4.0=8.0: 1.0*2.0=2.0
    //
    // Combined:
    // t=0.0: 2.0
    // t=2.0: 4.0+6.0=10.0
    // t=4.0: 2.0+12.0+2.0=16.0
    // t=6.0: 6.0+4.0=10.0
    // t=8.0: 2.0

    match actual {
      Distribution::Function(f) => {
        let expected_x = array![0.0, 2.0, 4.0, 6.0, 8.0];
        let expected_y = array![2.0, 10.0, 16.0, 10.0, 2.0];
        pretty_assert_ulps_eq!(f.t(), &expected_x, max_ulps = 10);
        pretty_assert_ulps_eq!(f.y(), &expected_y, max_ulps = 10);
      },
      _ => panic!("Expected Distribution::Function, got: {actual:#?}"),
    }
  }

  #[test]
  fn test_convolution_function_function_single_vs_multiple_points() {
    // Test 3: Single-point function convolving with multi-point function
    let x1 = array![7.0];
    let y1 = array![2.0];
    let f1 = Distribution::function(x1, y1).unwrap();

    let x2 = array![1.0, 4.0, 7.0, 10.0];
    let y2 = array![3.0, 6.0, 9.0, 12.0];
    let f2 = Distribution::function(x2, y2).unwrap();

    let actual = distribution_convolution(&f1, &f2).unwrap();

    // Manual calculation:
    // f1: (7.0, 2.0)
    // f2: (1.0, 3.0), (4.0, 6.0), (7.0, 9.0), (10.0, 12.0)
    //
    // Convolution:
    // 7.0+1.0=8.0: 2.0*3.0=6.0
    // 7.0+4.0=11.0: 2.0*6.0=12.0
    // 7.0+7.0=14.0: 2.0*9.0=18.0
    // 7.0+10.0=17.0: 2.0*12.0=24.0

    match actual {
      Distribution::Function(f) => {
        let expected_x = array![8.0, 11.0, 14.0, 17.0];
        let expected_y = array![6.0, 12.0, 18.0, 24.0];
        pretty_assert_ulps_eq!(f.t(), &expected_x, max_ulps = 10);
        pretty_assert_ulps_eq!(f.y(), &expected_y, max_ulps = 10);
      },
      _ => panic!("Expected Distribution::Function, got: {actual:#?}"),
    }
  }
}
