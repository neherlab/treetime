use crate::Distribution;
use crate::y_axis_policy::YAxisPolicy;
use eyre::Report;
use ndarray::Array1;

pub fn distribution_map<Y: YAxisPolicy, F>(dist: &Distribution<Y>, f: F) -> Result<Distribution<Y>, Report>
where
  F: Fn(f64) -> f64,
{
  match dist {
    Distribution::Function(func) => Distribution::function(func.t(), func.y().mapv(&f)),
    Distribution::Point(point) => Ok(Distribution::point(point.t(), f(point.amplitude()))),
    Distribution::Range(range) => Ok(Distribution::range((range.start(), range.end()), f(range.amplitude()))),
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
