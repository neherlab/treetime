use crate::distribution::distribution::Distribution;
use eyre::Report;
use ndarray::Array1;

/// Applies a function element-wise to the y-values of a distribution.
///
/// Creates a new Function distribution with the same time points but transformed y-values.
pub fn distribution_map<F>(dist: &Distribution, f: F) -> Result<Distribution, Report>
where
  F: Fn(f64) -> f64,
{
  match dist {
    Distribution::Function(func) => {
      let new_y = func.y().mapv(&f);
      Distribution::function(func.t().clone(), new_y)
    }
    Distribution::Point(_) | Distribution::Range(_) => {
      treetime_utils::make_error!("Map operation not supported for Point or Range distributions")
    }
    Distribution::Empty => Ok(Distribution::empty()),
  }
}

/// Creates a distribution from time points and a function that generates y-values.
pub fn distribution_from_fn<F>(time_points: Array1<f64>, f: F) -> Result<Distribution, Report>
where
  F: Fn(f64) -> f64,
{
  let y_values = time_points.mapv(&f);
  Distribution::function(time_points, y_values)
}
