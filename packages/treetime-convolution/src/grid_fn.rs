use crate::InterpElem;
use eyre::Report;
use ndarray::{Array1, Ix1, OwnedRepr};
use ndarray_interp::interp1d::{Interp1D, Interp1DBuilder, Linear};
use serde::{Deserialize, Deserializer, Serialize, Serializer};

/// Function represented on a regular grid
#[derive(Debug)]
pub struct GridFn<T: InterpElem> {
  interp: Interp1D<OwnedRepr<T>, OwnedRepr<T>, Ix1, Linear>,
}

impl<T: InterpElem> Clone for GridFn<T> {
  fn clone(&self) -> Self {
    // Reconstruct from existing data since Interp1D doesn't implement Clone
    Self::new(self.x().clone(), self.y().clone()).expect("Clone should not fail for valid GridFn")
  }
}

impl<T: InterpElem> Serialize for GridFn<T> {
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

impl<'de> Deserialize<'de> for GridFn<f64> {
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
      type Value = GridFn<f64>;

      fn expecting(&self, formatter: &mut fmt::Formatter) -> fmt::Result {
        formatter.write_str("struct GridFn")
      }

      fn visit_map<V>(self, mut map: V) -> Result<GridFn<f64>, V::Error>
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

impl<T: InterpElem> GridFn<T> {
  /// Create new grid function
  pub fn new(x: Array1<T>, y: Array1<T>) -> Result<Self, Report> {
    assert_eq!(x.len(), y.len(), "x and y arrays must have same length");
    assert!(x.len() >= 2, "Arrays must have at least 2 points");

    let interp = Interp1DBuilder::new(y).x(x).build()?;
    Ok(Self { interp })
  }

  /// Get x values
  pub fn x(&self) -> &Array1<T> {
    self.interp.x()
  }

  /// Get y values
  pub fn y(&self) -> &Array1<T> {
    self.interp.data()
  }

  /// Interpolate value at given x using linear interpolation
  pub fn interp(&self, xi: T) -> Result<T, Report> {
    let result = self.interp.interp_scalar(xi)?;
    Ok(result)
  }

  /// Interpolate at multiple points
  pub fn interp_many(&self, x_vals: &Array1<T>) -> Result<Array1<T>, Report> {
    let mut result = Array1::zeros(x_vals.len());
    for (i, &x) in x_vals.iter().enumerate() {
      result[i] = self.interp.interp_scalar(x)?;
    }
    Ok(result)
  }
}

impl GridFn<f64> {
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
}
