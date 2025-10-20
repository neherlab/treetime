use eyre::Report;
use ndarray::Array1;
use ndarray_interp::interp1d::{Interp1D, Interp1DBuilder, Linear};
use serde::{Deserialize, Deserializer, Serialize, Serializer};

/// Function represented on a regular grid
#[derive(Debug)]
pub struct GridFn {
  interp: Interp1D<ndarray::OwnedRepr<f64>, ndarray::OwnedRepr<f64>, ndarray::Ix1, Linear>,
}

impl Clone for GridFn {
  fn clone(&self) -> Self {
    // Reconstruct from existing data since Interp1D doesn't implement Clone
    Self::new(self.x().clone(), self.y().clone()).expect("Clone should not fail for valid GridFn")
  }
}

impl Serialize for GridFn {
  fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
  where
    S: Serializer,
  {
    use serde::ser::SerializeStruct;
    let mut state = serializer.serialize_struct("GridFn", 2)?;
    state.serialize_field("x", self.x())?;
    state.serialize_field("y", self.y())?;
    state.end()
  }
}

impl<'de> Deserialize<'de> for GridFn {
  fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
  where
    D: Deserializer<'de>,
  {
    use serde::de::{self, MapAccess, Visitor};
    use std::fmt;

    #[derive(Deserialize)]
    #[serde(field_identifier, rename_all = "lowercase")]
    enum Field {
      X,
      Y,
    }

    struct GridFnVisitor;

    impl<'de> Visitor<'de> for GridFnVisitor {
      type Value = GridFn;

      fn expecting(&self, formatter: &mut fmt::Formatter) -> fmt::Result {
        formatter.write_str("struct GridFn")
      }

      fn visit_map<V>(self, mut map: V) -> Result<GridFn, V::Error>
      where
        V: MapAccess<'de>,
      {
        let mut x: Option<Array1<f64>> = None;
        let mut y: Option<Array1<f64>> = None;
        while let Some(key) = map.next_key()? {
          match key {
            Field::X => {
              if x.is_some() {
                return Err(de::Error::duplicate_field("x"));
              }
              x = Some(map.next_value()?);
            },
            Field::Y => {
              if y.is_some() {
                return Err(de::Error::duplicate_field("y"));
              }
              y = Some(map.next_value()?);
            },
          }
        }
        let x = x.ok_or_else(|| de::Error::missing_field("x"))?;
        let y = y.ok_or_else(|| de::Error::missing_field("y"))?;

        GridFn::new(x, y).map_err(|e| de::Error::custom(format!("Failed to create GridFn: {e}")))
      }
    }

    const FIELDS: &[&str] = &["x", "y"];
    deserializer.deserialize_struct("GridFn", FIELDS, GridFnVisitor)
  }
}

impl GridFn {
  /// Create new grid function
  pub fn new(x: Array1<f64>, y: Array1<f64>) -> Result<Self, Report> {
    assert_eq!(x.len(), y.len(), "x and y arrays must have same length");
    assert!(x.len() >= 2, "Arrays must have at least 2 points");

    let interp = Interp1DBuilder::new(y).x(x).build()?;
    Ok(Self { interp })
  }

  /// Create grid function from domain and y-calculation function
  pub fn from_grid<F>((x_min, x_max): (f64, f64), dx: f64, y_fn: F) -> Result<Self, Report>
  where
    F: Fn(f64) -> f64,
  {
    assert!(dx > 0.0, "dx must be positive");
    assert!(x_max > x_min, "domain max must be greater than min");
    assert!(x_min.is_finite(), "x_min must be finite");
    assert!(x_max.is_finite(), "x_max must be finite");
    assert!(dx.is_finite(), "dx must be finite");

    let n_points = ((x_max - x_min) / dx + 1.0).round() as usize;
    let x = Array1::linspace(x_min, x_max, n_points);
    let y: Array1<f64> = x.mapv(y_fn);
    Self::new(x, y)
  }

  /// Get x values
  pub fn x(&self) -> &Array1<f64> {
    self.interp.x()
  }

  /// Get y values
  pub fn y(&self) -> &Array1<f64> {
    self.interp.data()
  }

  /// Interpolate value at given x using linear interpolation
  pub fn interp(&self, xi: f64) -> Result<f64, Report> {
    let result = self.interp.interp_scalar(xi)?;
    Ok(result)
  }

  /// Interpolate at multiple points
  pub fn interp_many(&self, x_vals: &Array1<f64>) -> Result<Array1<f64>, Report> {
    let mut result = Array1::zeros(x_vals.len());
    for (i, &x) in x_vals.iter().enumerate() {
      result[i] = self.interp.interp_scalar(x)?;
    }
    Ok(result)
  }

  /// Evaluate function on array of x values (like Python f(x_array))
  pub fn eval_array(&self, x_vals: &Array1<f64>) -> Result<Array1<f64>, Report> {
    let mut result = Array1::zeros(x_vals.len());
    for (i, &x) in x_vals.iter().enumerate() {
      result[i] = self.interp(x).unwrap_or(0.0); // Use 0 for out-of-domain
    }
    Ok(result)
  }

  /// Create function from a shared grid using a computation function
  /// This avoids interpolation by directly computing on the same grid
  pub fn from_shared_grid(x_grid: &Array1<f64>, compute_fn: impl Fn(f64) -> f64) -> Result<Self, Report> {
    let y_vals = x_grid.mapv(compute_fn);
    Self::new(x_grid.clone(), y_vals)
  }

  /// Evaluate on same grid as this function (no interpolation needed)
  pub fn eval_on_same_grid(&self, compute_fn: impl Fn(f64) -> f64) -> Result<GridFn, Report> {
    let y_vals = self.x().mapv(compute_fn);
    Self::new(self.x().clone(), y_vals)
  }

  /// Check if two GridFn have identical grids
  pub fn has_same_grid(&self, other: &GridFn) -> bool {
    self.x().len() == other.x().len()
      && self
        .x()
        .iter()
        .zip(other.x().iter())
        .all(|(a, b)| (a - b).abs() < 1e-15)
  }

  /// Get value at grid index (no interpolation)
  pub fn value_at_index(&self, index: usize) -> Option<f64> {
    self.y().get(index).copied()
  }

  /// Get x value at grid index
  pub fn x_at_index(&self, index: usize) -> Option<f64> {
    self.x().get(index).copied()
  }

  /// Get grid size
  pub fn grid_size(&self) -> usize {
    self.x().len()
  }

  /// Get maximum value
  pub fn max_value(&self) -> f64 {
    self.y().fold(f64::NEG_INFINITY, |acc, &x| acc.max(x))
  }

  /// Get x position of maximum value
  pub fn max_position(&self) -> f64 {
    let max_idx = self
      .y()
      .iter()
      .enumerate()
      .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
      .map_or(0, |(idx, _)| idx);
    self.x()[max_idx]
  }
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::pretty_assert_ulps_eq;
  use approx::assert_ulps_eq;
  use ndarray::array;

  #[test]
  fn test_gridfn_from_grid_basic() -> Result<(), Report> {
    let domain = (0.0, 2.0);
    let dx = 1.0;
    let result = GridFn::from_grid(domain, dx, |x| x * x)?;
    pretty_assert_ulps_eq!(result.x(), &array![0.0, 1.0, 2.0], max_ulps = 1);
    pretty_assert_ulps_eq!(result.y(), &array![0.0, 1.0, 4.0], max_ulps = 1);
    Ok(())
  }

  #[test]
  fn test_gridfn_interp_basic() -> Result<(), Report> {
    let x = Array1::from_vec(vec![0.0, 1.0, 2.0]);
    let y = Array1::from_vec(vec![0.0, 1.0, 4.0]);
    let result = GridFn::new(x, y)?;

    // Test exact points
    assert_ulps_eq!(result.interp(0.0)?, 0.0, max_ulps = 1);
    assert_ulps_eq!(result.interp(1.0)?, 1.0, max_ulps = 1);
    assert_ulps_eq!(result.interp(2.0)?, 4.0, max_ulps = 1);

    // Test interpolation
    assert_ulps_eq!(result.interp(0.5)?, 0.5, max_ulps = 1);
    assert_ulps_eq!(result.interp(1.5)?, 2.5, max_ulps = 1);

    Ok(())
  }

  #[test]
  fn test_gridfn_max_value_and_position() -> Result<(), Report> {
    let x = Array1::from_vec(vec![0.0, 1.0, 2.0, 3.0]);
    let y = Array1::from_vec(vec![1.0, 3.0, 5.0, 2.0]);
    let result = GridFn::new(x, y)?;
    assert_ulps_eq!(result.max_value(), 5.0, max_ulps = 1);
    assert_ulps_eq!(result.max_position(), 2.0, max_ulps = 1);
    Ok(())
  }
}
