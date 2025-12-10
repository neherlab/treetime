use crate::distribution::distribution::Distribution;
use crate::distribution::y_axis_policy::YAxisPolicy;
use eyre::Report;
use ndarray::Array1;

pub fn distribution_resample<Y: YAxisPolicy>(
  dist: &Distribution<Y>,
  time_points: &Array1<f64>,
) -> Result<Distribution<Y>, Report> {
  let values = dist.eval_many(time_points)?;
  Distribution::function(time_points.clone(), values)
}
