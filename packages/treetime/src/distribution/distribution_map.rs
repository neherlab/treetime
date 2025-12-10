use crate::distribution::distribution::Distribution;
use crate::distribution::y_axis_policy::YAxisPolicy;
use eyre::Report;
use ndarray::Array1;

pub fn distribution_map<Y: YAxisPolicy, F>(dist: &Distribution<Y>, f: F) -> Result<Distribution<Y>, Report>
where
  F: Fn(f64) -> f64,
{
  match dist {
    Distribution::Function(func) => {
      let new_y = func.y().mapv(&f);
      Distribution::function(func.t(), new_y)
    },
    Distribution::Point(point) => {
      let amplitude = f(point.amplitude());
      Ok(Distribution::point(point.t(), amplitude))
    },
    Distribution::Range(range) => {
      let amplitude = f(range.amplitude());
      Ok(Distribution::range((range.start(), range.end()), amplitude))
    },
    Distribution::Empty => Ok(Distribution::empty()),
    Distribution::Formula(_) => panic!("Map not implemented for Formula distributions"),
  }
}

pub fn distribution_from_fn<Y: YAxisPolicy, F>(time_points: Array1<f64>, f: F) -> Result<Distribution<Y>, Report>
where
  F: Fn(f64) -> f64,
{
  let y_values = time_points.mapv(&f);
  Distribution::function(time_points, y_values)
}
