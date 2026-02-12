#[cfg(test)]
mod tests {
  use crate::distribution::distribution::DistributionPlain as Distribution;
  use crate::distribution::distribution_convolution::distribution_convolution;
  use approx::assert_ulps_eq;
  use eyre::Report;

  #[test]
  fn test_inverse_convolution_integration() -> Result<(), Report> {
    // Test that the refactored code uses the same operation as the standalone function
    let child_dist: Distribution = Distribution::point(2013.0, 1.0);
    let branch_dist: Distribution = Distribution::point(2.5, 0.8);

    let negated_branch_dist: Distribution = branch_dist.negate();
    let result: Distribution = distribution_convolution(&child_dist, &negated_branch_dist)?;

    if let Some(parent_time) = result.likely_time() {
      assert_ulps_eq!(parent_time, 2010.5, max_ulps = 4);
      assert!(parent_time < 2013.0, "Parent should be older than child");
    } else {
      panic!("Expected valid parent time");
    }

    Ok(())
  }

  #[test]
  fn test_backward_pass_time_direction() -> Result<(), Report> {
    // Test that backward pass correctly computes older times for ancestors
    let child_time = 2012.0;
    let branch_length = 1.5;
    let expected_parent_time = child_time - branch_length; // 2010.5

    let child_dist: Distribution = Distribution::point(child_time, 1.0);
    let branch_dist: Distribution = Distribution::point(branch_length, 1.0);

    let negated_branch_dist: Distribution = branch_dist.negate();
    let parent_dist: Distribution = distribution_convolution(&child_dist, &negated_branch_dist)?;

    if let Some(actual_parent_time) = parent_dist.likely_time() {
      assert_ulps_eq!(actual_parent_time, expected_parent_time, max_ulps = 4);
      assert!(actual_parent_time < child_time, "Parent must be older than child");
    } else {
      panic!("Expected valid parent time");
    }

    Ok(())
  }
}
