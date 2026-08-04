#[cfg(test)]
mod tests {
  use crate::distribution_core::function::DistributionFunction;
  use crate::policy::{NegLog, Plain};
  use eyre::Report;
  use ndarray::array;
  use treetime_utils::pretty_assert_ulps_eq;

  type DistFn = DistributionFunction<f64, Plain>;
  type DistFnNegLog = DistributionFunction<f64, NegLog>;

  #[test]
  fn test_distribution_function_interpolate() -> Result<(), Report> {
    let f: DistFn = DistributionFunction::from_range_values((0.0, 3.0), array![0.0, 10.0, 20.0, 30.0])?;
    pretty_assert_ulps_eq!(f.interp(1.5)?, 15.0);
    Ok(())
  }

  #[test]
  fn test_distribution_function_interpolate_many() -> Result<(), Report> {
    let f: DistFn = DistributionFunction::from_range_values((0.0, 3.0), array![0.0, 10.0, 20.0, 30.0])?;

    let query = array![1.5, 2.25, 2.7];
    let expected = array![15.0, 22.5, 27.0];

    pretty_assert_ulps_eq!(f.interp_many(&query)?, expected);
    Ok(())
  }

  // Oracle: the two functions represent the same distribution on grid [0, 1, 2, 3]. Its mode is
  // at t = 1.0 (probability 5.0, the largest). Under `Plain` the ordinate is probability, so the
  // mode is the maximum ordinate. Under `NegLog` the ordinate is -ln(probability), so the mode is
  // the minimum ordinate. Selecting the opposite extremum returns t = 0.0 (probability 1.0), so
  // the asymmetric peak position discriminates the extremum direction.

  #[test]
  fn test_distribution_function_likely_time_plain_selects_max_ordinate() -> Result<(), Report> {
    let f: DistFn = DistributionFunction::from_range_values((0.0, 3.0), array![1.0, 5.0, 2.0, 4.0])?;
    pretty_assert_ulps_eq!(f.likely_time().unwrap(), 1.0);
    Ok(())
  }

  #[test]
  fn test_distribution_function_likely_time_neglog_selects_min_ordinate() -> Result<(), Report> {
    let neg_log = array![1.0, 5.0, 2.0, 4.0].mapv(|p: f64| -p.ln());
    let f: DistFnNegLog = DistributionFunction::from_range_values((0.0, 3.0), neg_log)?;
    pretty_assert_ulps_eq!(f.likely_time().unwrap(), 1.0);
    Ok(())
  }
}
