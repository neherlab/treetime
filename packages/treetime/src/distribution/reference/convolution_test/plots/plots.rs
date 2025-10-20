use crate::distribution::reference::convolution_test::framework::results::TestRunOutcome;
use crate::distribution::reference::convolution_test::framework::test_case::TestCase;
use crate::distribution::reference::convolution_test::plots::error_plots::{
  plot_absolute_error, plot_error_histogram, plot_tolerance_metrics,
};
use crate::distribution::reference::convolution_test::plots::functions_and_convolution::plot_functions_and_convolution;
use crate::distribution::reference::convolution_test::plots::pointwise_plots::{
  plot_derivative_errors, plot_pointwise_error_profiles,
};
use crate::distribution::reference::convolution_test::plots::spatial_plots::plot_spatial_profiles;
use eyre::Report;
use std::fs;

pub fn generate_plot_outputs<T>(output_dir: &str, outcomes: &[TestRunOutcome<T>]) -> Result<(), Report>
where
  T: TestCase,
{
  for outcome in outcomes {
    if let TestRunOutcome::Success(result) = outcome {
      let algorithm_dir = format!(
        "{}/{}/{}",
        output_dir,
        result.test_case.name(),
        result.algorithm.to_lowercase()
      );
      fs::create_dir_all(&algorithm_dir)?;

      plot_functions_and_convolution(result, &algorithm_dir)?;
      plot_absolute_error(result, &algorithm_dir)?;
      plot_tolerance_metrics(result, &algorithm_dir)?;
      plot_pointwise_error_profiles(result, &algorithm_dir)?;
      plot_derivative_errors(result, &algorithm_dir)?;
      plot_spatial_profiles(result, &algorithm_dir)?;
      plot_error_histogram(result, &algorithm_dir)?;
    }
  }
  Ok(())
}
