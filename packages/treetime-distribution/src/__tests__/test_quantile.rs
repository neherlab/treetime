#[cfg(test)]
mod tests {
  use crate::Distribution;
  use crate::y_axis_policy::Plain;
  use approx::assert_relative_eq;
  use ndarray::array;

  #[test]
  fn test_quantile_invalid_p_returns_none() {
    let dist = Distribution::<Plain>::point(5.0, 1.0);
    assert!(dist.quantile(-0.1).is_none());
    assert!(dist.quantile(1.1).is_none());
    assert!(dist.quantile(f64::NAN).is_none());
  }

  #[test]
  fn test_quantile_empty_returns_none() {
    let empty = Distribution::<Plain>::Empty;
    assert!(empty.quantile(0.5).is_none());
  }

  #[test]
  fn test_quantile_point_returns_point_location() {
    let point = Distribution::<Plain>::point(7.5, 1.0);
    assert_relative_eq!(point.quantile(0.0).unwrap(), 7.5);
    assert_relative_eq!(point.quantile(0.5).unwrap(), 7.5);
    assert_relative_eq!(point.quantile(1.0).unwrap(), 7.5);
  }

  #[test]
  fn test_quantile_range_linear_interpolation() {
    let range = Distribution::<Plain>::range((10.0, 20.0), 1.0);
    assert_relative_eq!(range.quantile(0.0).unwrap(), 10.0);
    assert_relative_eq!(range.quantile(0.5).unwrap(), 15.0);
    assert_relative_eq!(range.quantile(1.0).unwrap(), 20.0);
    assert_relative_eq!(range.quantile(0.25).unwrap(), 12.5);
  }

  #[test]
  fn test_quantile_uniform_function() {
    // Uniform distribution over [0, 4] with constant density
    let t = array![0.0, 1.0, 2.0, 3.0, 4.0];
    let y = array![1.0, 1.0, 1.0, 1.0, 1.0];
    let func = Distribution::<Plain>::function(t, y).unwrap();

    assert_relative_eq!(func.quantile(0.0).unwrap(), 0.0);
    assert_relative_eq!(func.quantile(1.0).unwrap(), 4.0);
    assert_relative_eq!(func.quantile(0.5).unwrap(), 2.0, epsilon = 1e-10);
    assert_relative_eq!(func.quantile(0.25).unwrap(), 1.0, epsilon = 1e-10);
  }

  #[test]
  fn test_quantile_single_point_function() {
    let t = array![5.0];
    let y = array![1.0];
    let func = Distribution::<Plain>::function(t, y).unwrap();
    assert_relative_eq!(func.quantile(0.5).unwrap(), 5.0);
  }

  #[test]
  fn test_quantile_function_median() {
    // Triangular distribution peaking at center
    let t = array![0.0, 1.0, 2.0];
    let y = array![0.0, 2.0, 0.0];
    let func = Distribution::<Plain>::function(t, y).unwrap();

    // Median of symmetric triangle at center
    assert_relative_eq!(func.quantile(0.5).unwrap(), 1.0, epsilon = 1e-10);
  }

  #[test]
  fn test_quantile_function_boundaries() {
    let t = array![0.0, 1.0, 2.0, 3.0];
    let y = array![1.0, 2.0, 2.0, 1.0];
    let func = Distribution::<Plain>::function(t, y).unwrap();

    // p=0 returns first point
    assert_relative_eq!(func.quantile(0.0).unwrap(), 0.0);
    // p=1 returns last point
    assert_relative_eq!(func.quantile(1.0).unwrap(), 3.0);
  }

  #[test]
  fn test_confidence_interval_95_percent() {
    // Uniform distribution - 95% CI should be [0.025*range, 0.975*range]
    let range = Distribution::<Plain>::range((0.0, 100.0), 1.0);
    let ci = range.confidence_interval(0.025, 0.975).unwrap();
    assert_relative_eq!(ci.0, 2.5);
    assert_relative_eq!(ci.1, 97.5);
  }

  #[test]
  fn test_confidence_interval_empty_returns_none() {
    let empty = Distribution::<Plain>::Empty;
    assert!(empty.confidence_interval(0.025, 0.975).is_none());
  }
}
