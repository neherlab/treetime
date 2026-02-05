use crate::distribution::distribution_function::DistributionFunction;
use crate::distribution::y_axis_policy::Plain;
use crate::pretty_assert_ulps_eq;
use eyre::Report;
use ndarray::array;

type DistFn = DistributionFunction<f64, Plain>;

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
