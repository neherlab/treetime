#[cfg(test)]
mod tests {
  use crate::distribution::distribution::Distribution;
  use crate::distribution::y_axis_policy::Plain;
  use approx::assert_relative_eq;
  use ndarray::array;

  #[test]
  fn test_distribution_max_value() {
    let point = Distribution::<Plain>::point(1.0, 5.0);
    assert_relative_eq!(point.max_value(), 5.0);

    let range = Distribution::<Plain>::range((0.0, 2.0), 3.0);
    assert_relative_eq!(range.max_value(), 3.0);

    let func = Distribution::<Plain>::function(array![0.0, 1.0, 2.0], array![1.0, 4.0, 2.0]).unwrap();
    assert_relative_eq!(func.max_value(), 4.0);

    let empty = Distribution::<Plain>::Empty;
    assert_relative_eq!(empty.max_value(), 0.0);
  }

  #[test]
  fn test_distribution_scale_by() {
    let point = Distribution::<Plain>::point(1.0, 5.0);
    let scaled = point.scale_by(2.0);
    if let Distribution::Point(p) = scaled {
      assert_relative_eq!(p.amplitude(), 10.0);
    } else {
      panic!("Expected Point");
    }

    let range = Distribution::<Plain>::range((0.0, 2.0), 3.0);
    let scaled = range.scale_by(0.5);
    if let Distribution::Range(r) = scaled {
      assert_relative_eq!(r.amplitude(), 1.5);
    } else {
      panic!("Expected Range");
    }
  }

  #[test]
  fn test_distribution_normalize() {
    let func = Distribution::<Plain>::function(array![0.0, 1.0, 2.0], array![2.0, 8.0, 4.0]).unwrap();
    let normalized = func.normalize();

    if let Distribution::Function(f) = normalized {
      assert_relative_eq!(f.y()[0], 0.25);
      assert_relative_eq!(f.y()[1], 1.0);
      assert_relative_eq!(f.y()[2], 0.5);
    } else {
      panic!("Expected Function");
    }
  }

  #[test]
  fn test_distribution_normalize_empty_on_zero() {
    let func = Distribution::<Plain>::function(array![0.0, 1.0, 2.0], array![0.0, 0.0, 0.0]).unwrap();
    let normalized = func.normalize();
    assert!(matches!(normalized, Distribution::Empty));
  }
}
