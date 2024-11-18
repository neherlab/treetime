use crate::distribution::distribution::Distribution;
use crate::distribution::distribution_function::DistributionFunction;
use crate::distribution::distribution_point::DistributionPoint;
use crate::distribution::distribution_range::DistributionRange;
use approx::ulps_eq;
use eyre::Report;
use ndarray::array;

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
      Ok(convolution_point_function(a, b)) //
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
    Distribution::general(x, y)
  } else {
    let x = array![start, peak_start, peak_end, end];
    let y = array![0.0, peak_amplitude, peak_amplitude, 0.0];
    Distribution::general(x, y)
  }
}

fn convolution_point_range(p: &DistributionPoint<f64>, r: &DistributionRange<f64>) -> Distribution {
  let begin = r.start() - p.t();
  let end = r.end() - p.t();
  let amplitude = p.amplitude() * r.amplitude();
  Distribution::range((begin, end), amplitude)
}

fn convolution_point_function(p: &DistributionPoint<f64>, f: &DistributionFunction<f64>) -> Distribution {
  unimplemented!()
}

fn convolution_range_function(r: &DistributionRange<f64>, f: &DistributionFunction<f64>) -> Distribution {
  unimplemented!()
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
    let b = Distribution::general(array![], array![]).unwrap();
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
      Distribution::general(x, y).unwrap()
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
      Distribution::general(x, y).unwrap()
    };
    assert_eq!(expected, actual);
  }
}
