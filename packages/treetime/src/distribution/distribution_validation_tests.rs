use crate::distribution::distribution::Distribution;
use crate::distribution::distribution_convolution::distribution_convolution;
use crate::distribution::distribution_function::DistributionFunction;
use crate::distribution::y_axis_policy::Plain;
use ndarray::Array1;
use treetime_convolution::testing::metrics::metrics::ConvolutionMetrics;
use treetime_convolution::testing::test_suites::gaussian::GaussianTestSuite;
use treetime_convolution::testing::test_suites::test_suites::TestSuite;

#[test]
fn test_gaussian_convolution_validation() {
  let suite = GaussianTestSuite;
  let test_cases = suite.create_test_cases();

  for case in test_cases {
    println!("Running validation test case: {}", case.name);

    let (grid_min, grid_max) = case.input_grid_domain;
    let n_points = case.input_grid_n_points;
    let input_grid = Array1::linspace(grid_min, grid_max, n_points);
    let dx = input_grid[1] - input_grid[0];

    // Create inputs
    let f_vals = suite.create_f(&case, &input_grid).unwrap();
    let g_vals = suite.create_g(&case, &input_grid).unwrap();

    let dist_f = Distribution::<Plain>::Function(
      DistributionFunction::from_start_dx_values(grid_min, dx, f_vals).unwrap(),
    );
    let dist_g = Distribution::<Plain>::Function(
      DistributionFunction::from_start_dx_values(grid_min, dx, g_vals).unwrap(),
    );

    // Run convolution
    let start = std::time::Instant::now();
    let result = distribution_convolution(&dist_f, &dist_g).unwrap();
    let duration = start.elapsed().as_secs_f64() * 1000.0;

    // Create evaluation grid (2x extent to cover convolution result)
    // The convolution result size in Distribution is dynamic, but for Gaussian comparison
    // we use the standard evaluation grid from the suite logic or similar.
    // The suite usually evaluates on 2x domain.
    let eval_min = grid_min * 2.0;
    let eval_max = grid_max * 2.0;
    // Expected result len is (N + N - 1), but let's stick to what the suite usually does for metrics.
    // However, analytical_convolution can evaluate anywhere.
    // Let's use the result's own grid for evaluation to be fair, or resampling?
    
    // In treetime-convolution tests:
    // let (evaluation_grid_min, evaluation_grid_max) = (input_grid_min * 2.0, input_grid_max * 2.0);
    // let evaluation_grid_n_points = 2 * input_grid_n_points - 1;
    let eval_n_points = 2 * n_points - 1;
    let eval_grid = Array1::linspace(eval_min, eval_max, eval_n_points);

    // Analytical expectation
    let expected = suite.analytical_convolution(&case, &eval_grid).unwrap();

    // Actual values interpolated on evaluation grid
    let actual = result.eval_many(&eval_grid).unwrap();

    // Compute metrics
    let metrics = ConvolutionMetrics::new(&eval_grid, &actual, &expected, duration).unwrap();

    println!("  R²: {:.8}", metrics.aggregate.domain_agreement.quality_metrics.r_squared);
    println!("  RMSE: {:.8e}", metrics.aggregate.domain_agreement.quality_metrics.rmse);
    println!("  Max Error: {:.8e}", metrics.pointwise.errors.summary.abs_max);

    // Assertions
    // We expect very high accuracy for Gaussians
    assert!(
      metrics.aggregate.domain_agreement.quality_metrics.r_squared > 0.99, 
      "R² too low for case {}: {}", case.name, metrics.aggregate.domain_agreement.quality_metrics.r_squared
    );
    
    // Some cases might be stressful (e.g. truncation), so we can be slightly lenient or check specific cases
    if case.stress_type.contains("truncation") || case.stress_type.contains("underflow") {
       // Allow lower precision for stress tests
       assert!(metrics.aggregate.domain_agreement.quality_metrics.r_squared > 0.90);
    }
  }
}
