use crate::distribution::distribution::Distribution;
use crate::distribution::distribution_function::DistributionFunction;
use crate::distribution::distribution_point::DistributionPoint;
use crate::distribution::distribution_range::DistributionRange;
use crate::utils::ndarray::zeros;
use approx::ulps_eq;
use eyre::Report;
use ndarray::{array, Array1};

pub fn distribution_convolution(a: &Distribution, b: &Distribution) -> Result<Distribution, Report> {
  match (a, b) {
    (Distribution::Empty, _) | (_, Distribution::Empty) => {
      Ok(Distribution::Empty) //
    }
    (Distribution::Point(a), Distribution::Point(b)) => {
      Ok(convolution_point_point(a, b)) //
    }
    (Distribution::Point(a), Distribution::Range(b)) | (Distribution::Range(b), Distribution::Point(a)) => {
      Ok(convolution_point_range(a, b)) //
    }
    (Distribution::Range(a), Distribution::Range(b)) => {
      convolution_range_range(a, b) //
    }
    (Distribution::Point(a), Distribution::Function(b)) | (Distribution::Function(b), Distribution::Point(a)) => {
      Ok(convolution_point_function(a, b)?) //
    }
    (Distribution::Range(a), Distribution::Function(b)) | (Distribution::Function(b), Distribution::Range(a)) => {
      Ok(convolution_range_function(a, b)) //
    }
    (Distribution::Function(a), Distribution::Function(b)) => {
      Ok(convolution_function_function(a, b)) //
    }
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
  let begin = r.start() - p.t();
  let end = r.end() - p.t();
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

fn convolution_function_function(a: &DistributionFunction<f64>, b: &DistributionFunction<f64>) -> Distribution {
  unimplemented!()
}

#[cfg(test)]
mod tests {
  use super::*;
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
}
