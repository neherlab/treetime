use crate::distribution::distribution::Distribution;
use eyre::Report;
use ndarray::Array1;

/// Resamples a distribution to specific time points by evaluating it at those points.
///
/// This creates a new Function distribution with the specified time points
/// and corresponding values obtained by evaluating the input distribution.
pub fn distribution_resample(dist: &Distribution, time_points: &Array1<f64>) -> Result<Distribution, Report> {
  let values = dist.eval_many(time_points)?;
  Distribution::function(time_points.clone(), values)
}
