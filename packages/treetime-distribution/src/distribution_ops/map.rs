use crate::Distribution;
use crate::policy::YAxisPolicy;
use eyre::Report;

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
