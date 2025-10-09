use crate::distribution::reference::convolution::{ndarray_convolve, riemann_convolve};
use crate::distribution::reference::domain_agreement_metrics::DomainAgreementMetrics;
use crate::distribution::reference::gaussian::{gaussian_convolution, gaussian_f, gaussian_g};
use crate::io::csv::CsvStructFileWriter;
use crate::io::json::{json_write_file, JsonPretty};
use crate::utils::float_fmt::float_to_significant_digits;
use eyre::Report;
use ndarray::Array1;
use serde::{Deserialize, Serialize};
use std::fmt;

/// Comprehensive test framework for convolution implementations
#[derive(Debug, Clone)]
pub struct ConvolutionTestFramework {
  pub test_cases: Vec<GaussianTestCase>,
  pub algorithms: Vec<ConvolutionAlgorithm>,
  pub output_dir: String,
}

/// Available convolution algorithms for testing
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum ConvolutionAlgorithm {
  Riemann,
  Ndarray,
}

impl fmt::Display for ConvolutionAlgorithm {
  fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
    match self {
      ConvolutionAlgorithm::Riemann => write!(f, "riemann"),
      ConvolutionAlgorithm::Ndarray => write!(f, "ndarray"),
    }
  }
}

/// Test case parameters for Gaussian convolution
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GaussianTestCase {
  pub name: String,
  pub description: String,
  pub sigma_f: f64,
  pub sigma_g: f64,
  pub mu: f64,
  pub f_domain: (f64, f64),
  pub g_domain: (f64, f64),
  pub eval_domain: (f64, f64),
  pub dx: f64,
  pub stress_type: String,
  pub analytical_caution: String,
}

/// Results for a single test case and algorithm combination
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TestResult {
  pub test_case_name: String,
  pub algorithm: ConvolutionAlgorithm,
  pub parameters: GaussianTestCase,
  pub metrics: DomainAgreementMetrics,
  pub execution_time_ms: f64,
  pub grid_points: usize,
  pub peak_error_location: f64,
  pub peak_error_value: f64,
}

/// Aggregate summary across all tests
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TestSummary {
  pub total_tests: usize,
  pub total_algorithms: usize,
  pub execution_time_total_ms: f64,
  pub algorithm_summaries: Vec<AlgorithmSummary>,
  pub overall_assessment: String,
}

/// Summary for a specific algorithm across all test cases
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AlgorithmSummary {
  pub algorithm: ConvolutionAlgorithm,
  pub test_cases_count: usize,
  pub execution_time_total_ms: f64,
  pub execution_time_avg_ms: f64,
  pub r2_min: f64,
  pub r2_max: f64,
  pub r2_mean: f64,
  pub max_abs_error_overall: f64,
  pub max_rel_error_overall: f64,
  pub passed_tests: usize,
  pub failed_tests: usize,
  pub success_rate: f64,
}

/// Serializable test result for TSV/CSV output
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FlatTestResult {
  pub test_case_name: String,
  pub algorithm: String,
  pub sigma_f: f64,
  pub sigma_g: f64,
  pub mu: f64,
  pub dx: f64,
  pub grid_points: usize,
  pub execution_time_ms: f64,
  pub r_squared: f64,
  pub max_abs_error: f64,
  pub mean_abs_error: f64,
  pub max_rel_error_percent: f64,
  pub mean_rel_error_percent: f64,
  pub rmse: f64,
  pub correlation: f64,
  pub mass_conservation_error: f64,
  pub peak_error_location: f64,
  pub peak_error_value: f64,
  pub tolerance_strict_abs: f64,
  pub tolerance_moderate_abs: f64,
  pub tolerance_loose_abs: f64,
  pub tolerance_strict_rel: f64,
  pub tolerance_moderate_rel: f64,
  pub tolerance_loose_rel: f64,
  pub stress_type: String,
  pub overall_assessment: String,
}

impl ConvolutionTestFramework {
  /// Create new test framework with default comprehensive test cases
  pub fn new(output_dir: String) -> Self {
    Self {
      test_cases: create_comprehensive_test_cases(),
      algorithms: vec![ConvolutionAlgorithm::Riemann, ConvolutionAlgorithm::Ndarray],
      output_dir,
    }
  }

  /// Create framework with custom test cases
  pub fn with_test_cases(test_cases: Vec<GaussianTestCase>, output_dir: String) -> Self {
    Self {
      test_cases,
      algorithms: vec![ConvolutionAlgorithm::Riemann, ConvolutionAlgorithm::Ndarray],
      output_dir,
    }
  }

  /// Add additional test case
  pub fn add_test_case(&mut self, test_case: GaussianTestCase) {
    self.test_cases.push(test_case);
  }

  /// Set algorithms to test
  pub fn set_algorithms(&mut self, algorithms: Vec<ConvolutionAlgorithm>) {
    self.algorithms = algorithms;
  }

  /// Run all test cases for all algorithms
  pub fn run_all_tests(&self) -> Result<Vec<TestResult>, Report> {
    let mut results = Vec::new();

    println!("=== Convolution Test Framework ===");
    println!(
      "Running {} test cases with {} algorithms ({} total tests)\n",
      self.test_cases.len(),
      self.algorithms.len(),
      self.test_cases.len() * self.algorithms.len()
    );

    for test_case in &self.test_cases {
      println!("Test Case: {}", test_case.name);
      println!("  Description: {}", test_case.description);
      println!("  Parameters: σ_f={}, σ_g={}, μ={}, dx={}",
        test_case.sigma_f, test_case.sigma_g, test_case.mu, test_case.dx);

      for &algorithm in &self.algorithms {
        print!("  Running {} algorithm... ", algorithm);
        let start_time = std::time::Instant::now();

        match self.run_single_test(test_case, algorithm) {
          Ok(result) => {
            let elapsed = start_time.elapsed().as_secs_f64() * 1000.0;
            println!("✓ ({:.1}ms, R²={:.6})", elapsed, result.metrics.quality_metrics.r_squared);
            results.push(result);
          }
          Err(e) => {
            println!("✗ Error: {}", e);
          }
        }
      }
      println!();
    }

    Ok(results)
  }

  /// Run single test case with specific algorithm
  pub fn run_single_test(
    &self,
    test_case: &GaussianTestCase,
    algorithm: ConvolutionAlgorithm,
  ) -> Result<TestResult, Report> {
    let start_time = std::time::Instant::now();

    // Create input functions
    let f = gaussian_f(test_case.sigma_f, test_case.f_domain, test_case.dx)?;
    let g = gaussian_g(test_case.sigma_g, test_case.mu, test_case.g_domain, test_case.dx)?;

    // Create evaluation grid
    let (eval_min, eval_max) = test_case.eval_domain;
    let n_eval_points = ((eval_max - eval_min) / test_case.dx + 1.0).round() as usize;
    let eval_grid = Array1::from_iter((0..n_eval_points).map(|i| eval_min + i as f64 * test_case.dx));

    // Run convolution
    let actual_result = match algorithm {
      ConvolutionAlgorithm::Riemann => riemann_convolve(&f, &g, &eval_grid)?,
      ConvolutionAlgorithm::Ndarray => ndarray_convolve(&f, &g, &eval_grid)?,
    };

    // Compute analytical expected result
    let expected_result = gaussian_convolution(
      test_case.sigma_f,
      test_case.sigma_g,
      test_case.mu,
      test_case.eval_domain,
      test_case.dx,
    )?;

    // Compute metrics
    let metrics = DomainAgreementMetrics::new(
      actual_result.x(),
      actual_result.y(),
      expected_result.y(),
    )?;

    let execution_time = start_time.elapsed().as_secs_f64() * 1000.0;

    // Find peak error location
    let abs_errors: Array1<f64> = actual_result
      .y()
      .iter()
      .zip(expected_result.y().iter())
      .map(|(&a, &e)| (a - e).abs())
      .collect();

    let max_error_idx = abs_errors
      .iter()
      .enumerate()
      .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
      .unwrap()
      .0;

    let peak_error_location = actual_result.x()[max_error_idx];
    let peak_error_value = abs_errors[max_error_idx];

    Ok(TestResult {
      test_case_name: test_case.name.clone(),
      algorithm,
      parameters: test_case.clone(),
      metrics,
      execution_time_ms: execution_time,
      grid_points: eval_grid.len(),
      peak_error_location,
      peak_error_value,
    })
  }

  /// Generate comprehensive summary from test results
  pub fn generate_summary(&self, results: &[TestResult]) -> TestSummary {
    let total_execution_time: f64 = results.iter().map(|r| r.execution_time_ms).sum();

    let mut algorithm_summaries = Vec::new();

    for &algorithm in &self.algorithms {
      let algo_results: Vec<_> = results.iter().filter(|r| r.algorithm == algorithm).collect();

      if algo_results.is_empty() {
        continue;
      }

      let r2_values: Vec<f64> = algo_results.iter().map(|r| r.metrics.quality_metrics.r_squared).collect();
      let execution_times: Vec<f64> = algo_results.iter().map(|r| r.execution_time_ms).collect();

      let r2_min = r2_values.iter().fold(f64::INFINITY, |a, &b| a.min(b));
      let r2_max = r2_values.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b));
      let r2_mean = r2_values.iter().sum::<f64>() / r2_values.len() as f64;

      let execution_time_total = execution_times.iter().sum::<f64>();
      let execution_time_avg = execution_time_total / execution_times.len() as f64;

      let max_abs_error_overall = algo_results
        .iter()
        .map(|r| r.metrics.abs_error_stats.max)
        .fold(0.0f64, |a, b| a.max(b));

      let max_rel_error_overall = algo_results
        .iter()
        .map(|r| r.metrics.rel_error_stats.max)
        .fold(0.0f64, |a, b| a.max(b));

      // Simple pass/fail based on R² > 0.95
      let passed_tests = algo_results.iter().filter(|r| r.metrics.quality_metrics.r_squared > 0.95).count();
      let failed_tests = algo_results.len() - passed_tests;
      let success_rate = passed_tests as f64 / algo_results.len() as f64;

      algorithm_summaries.push(AlgorithmSummary {
        algorithm,
        test_cases_count: algo_results.len(),
        execution_time_total_ms: execution_time_total,
        execution_time_avg_ms: execution_time_avg,
        r2_min,
        r2_max,
        r2_mean,
        max_abs_error_overall,
        max_rel_error_overall,
        passed_tests,
        failed_tests,
        success_rate,
      });
    }

    let overall_assessment = if algorithm_summaries.iter().all(|s| s.success_rate > 0.9) {
      "Excellent - All algorithms perform well".to_owned()
    } else if algorithm_summaries.iter().any(|s| s.success_rate > 0.8) {
      "Good - Some algorithms perform well".to_owned()
    } else {
      "Poor - Significant issues detected".to_owned()
    };

    TestSummary {
      total_tests: results.len(),
      total_algorithms: self.algorithms.len(),
      execution_time_total_ms: total_execution_time,
      algorithm_summaries,
      overall_assessment,
    }
  }

  /// Save results to JSON file
  pub fn save_results_json(&self, results: &[TestResult], summary: &TestSummary) -> Result<(), Report> {
    std::fs::create_dir_all(&self.output_dir)?;

    #[derive(Serialize)]
    struct FullResults<'a> {
      summary: &'a TestSummary,
      results: &'a [TestResult],
    }

    let full_results = FullResults { summary, results };
    let json_path = format!("{}/convolution_test_results.json", self.output_dir);
    json_write_file(&json_path, &full_results, JsonPretty(true))?;
    println!("Saved detailed JSON results to: {}", json_path);
    Ok(())
  }

  /// Save results to TSV file for analysis in Python/R
  pub fn save_results_tsv(&self, results: &[TestResult]) -> Result<(), Report> {
    std::fs::create_dir_all(&self.output_dir)?;

    let tsv_path = format!("{}/convolution_test_results.tsv", self.output_dir);
    let mut writer = CsvStructFileWriter::new(&tsv_path, b'\t')?;

    for result in results {
      let flat_result = FlatTestResult {
        test_case_name: result.test_case_name.clone(),
        algorithm: result.algorithm.to_string(),
        sigma_f: result.parameters.sigma_f,
        sigma_g: result.parameters.sigma_g,
        mu: result.parameters.mu,
        dx: result.parameters.dx,
        grid_points: result.grid_points,
        execution_time_ms: result.execution_time_ms,
        r_squared: result.metrics.quality_metrics.r_squared,
        max_abs_error: result.metrics.abs_error_stats.max,
        mean_abs_error: result.metrics.abs_error_stats.mean,
        max_rel_error_percent: result.metrics.rel_error_stats.max * 100.0,
        mean_rel_error_percent: result.metrics.rel_error_stats.mean * 100.0,
        rmse: result.metrics.quality_metrics.rmse,
        correlation: result.metrics.quality_metrics.correlation,
        mass_conservation_error: result.metrics.quality_metrics.mass_error,
        peak_error_location: result.peak_error_location,
        peak_error_value: result.peak_error_value,
        tolerance_strict_abs: result.metrics.abs_tolerance_fraction(0),
        tolerance_moderate_abs: result.metrics.abs_tolerance_fraction(1),
        tolerance_loose_abs: result.metrics.abs_tolerance_fraction(2),
        tolerance_strict_rel: result.metrics.rel_tolerance_fraction(0),
        tolerance_moderate_rel: result.metrics.rel_tolerance_fraction(1),
        tolerance_loose_rel: result.metrics.rel_tolerance_fraction(2),
        stress_type: result.parameters.stress_type.clone(),
        overall_assessment: format!("{}", result.metrics.overall_assessment()),
      };

      writer.write(&flat_result)?;
    }

    println!("Saved TSV results to: {}", tsv_path);
    Ok(())
  }

  /// Print comprehensive summary to console
  pub fn print_summary(&self, summary: &TestSummary) {
    println!("=== CONVOLUTION TEST FRAMEWORK SUMMARY ===\n");

    println!("Overall Statistics:");
    println!("  Total tests run: {}", summary.total_tests);
    println!("  Total algorithms: {}", summary.total_algorithms);
    println!("  Total execution time: {:.1}ms", summary.execution_time_total_ms);
    println!("  Assessment: {}\n", summary.overall_assessment);

    println!("Algorithm Performance:");
    println!("{:>10} {:>8} {:>10} {:>8} {:>8} {:>8} {:>12} {:>12} {:>8}",
      "Algorithm", "Tests", "Time(ms)", "R²Min", "R²Max", "R²Mean", "MaxAbsErr", "MaxRelErr%", "Success%");
    println!("{}", "-".repeat(100));

    for algo_summary in &summary.algorithm_summaries {
      println!("{:>10} {:>8} {:>10.1} {:>8.4} {:>8.4} {:>8.4} {:>12} {:>12.1} {:>8.1}",
        algo_summary.algorithm,
        algo_summary.test_cases_count,
        algo_summary.execution_time_total_ms,
        algo_summary.r2_min,
        algo_summary.r2_max,
        algo_summary.r2_mean,
        float_to_significant_digits(algo_summary.max_abs_error_overall, 3),
        algo_summary.max_rel_error_overall * 100.0,
        algo_summary.success_rate * 100.0,
      );
    }
    println!();
  }
}

/// Create comprehensive test cases based on convolution-gaussian.md reference
pub fn create_comprehensive_test_cases() -> Vec<GaussianTestCase> {
  vec![
    // Case 1 — Baseline, moderate widths
    GaussianTestCase {
      name: "baseline_moderate".to_owned(),
      description: "General correctness with asymmetric widths and shift".to_owned(),
      sigma_f: 1.0,
      sigma_g: 2.0,
      mu: 1.0,
      f_domain: (-8.0, 8.0),
      g_domain: (-12.0, 14.0),
      eval_domain: (-10.0, 12.0),
      dx: 0.01,
      stress_type: "none".to_owned(),
      analytical_caution: "none".to_owned(),
    },

    // Case 2 — Centered, unequal widths
    GaussianTestCase {
      name: "centered_unequal_widths".to_owned(),
      description: "Width mixing; peak alignment at x=0".to_owned(),
      sigma_f: 1.5,
      sigma_g: 0.5,
      mu: 0.0,
      f_domain: (-12.0, 12.0),
      g_domain: (-4.0, 4.0),
      eval_domain: (-8.0, 8.0),
      dx: 0.01,
      stress_type: "resolution vs narrow kernel".to_owned(),
      analytical_caution: "none".to_owned(),
    },

    // Case 3 — Sharp g, small shift (delta-like)
    GaussianTestCase {
      name: "delta_like_small_shift".to_owned(),
      description: "Sensitivity to sub-grid shifts; near-delta kernel behavior".to_owned(),
      sigma_f: 1.0,
      sigma_g: 0.05,
      mu: 0.2,
      f_domain: (-8.0, 8.0),
      g_domain: (-0.4, 1.0),
      eval_domain: (-4.0, 4.0),
      dx: 0.002,
      stress_type: "discretization error, normalization".to_owned(),
      analytical_caution: "potential overflow for large |x|/σ_g".to_owned(),
    },

    // Case 4 — Very wide g, negative shift
    GaussianTestCase {
      name: "wide_kernel_negative_shift".to_owned(),
      description: "Robustness to ultra-broad kernel tails; negative translation".to_owned(),
      sigma_f: 0.5,
      sigma_g: 5.0,
      mu: -1.0,
      f_domain: (-6.0, 6.0),
      g_domain: (-41.0, 39.0),
      eval_domain: (-20.0, 20.0),
      dx: 0.02,
      stress_type: "tail truncation, numeric underflow, FFT padding".to_owned(),
      analytical_caution: "underflow in far tails expected".to_owned(),
    },

    // Case 5 — Large positive shift
    GaussianTestCase {
      name: "large_positive_shift".to_owned(),
      description: "Translation invariance for large μ".to_owned(),
      sigma_f: 1.0,
      sigma_g: 1.2,
      mu: 10.0,
      f_domain: (-10.0, 10.0),
      g_domain: (-2.0, 22.0),
      eval_domain: (0.0, 20.0),
      dx: 0.01,
      stress_type: "circular convolution artifacts, insufficient padding".to_owned(),
      analytical_caution: "none".to_owned(),
    },

    // Case 6 — Large negative shift
    GaussianTestCase {
      name: "large_negative_shift".to_owned(),
      description: "Translation invariance for negative μ".to_owned(),
      sigma_f: 2.0,
      sigma_g: 1.0,
      mu: -7.0,
      f_domain: (-16.0, 16.0),
      g_domain: (-15.0, 1.0),
      eval_domain: (-20.0, 0.0),
      dx: 0.01,
      stress_type: "padding and index handling".to_owned(),
      analytical_caution: "none".to_owned(),
    },

    // Case 7 — Extreme ratio σ_g >> σ_f
    GaussianTestCase {
      name: "extreme_ratio_wide_g".to_owned(),
      description: "Convolution where g is nearly constant over f's support".to_owned(),
      sigma_f: 0.1,
      sigma_g: 10.0,
      mu: 3.0,
      f_domain: (-1.0, 1.0),
      g_domain: (-77.0, 83.0),
      eval_domain: (-5.0, 15.0),
      dx: 0.01,
      stress_type: "tail coverage, numerical cancellation".to_owned(),
      analytical_caution: "underflow in extreme tails".to_owned(),
    },

    // Case 8 — Extreme ratio σ_f >> σ_g
    GaussianTestCase {
      name: "extreme_ratio_wide_f".to_owned(),
      description: "Near-delta g translating f".to_owned(),
      sigma_f: 8.0,
      sigma_g: 0.2,
      mu: -2.5,
      f_domain: (-64.0, 64.0),
      g_domain: (-4.1, -0.9),
      eval_domain: (-30.0, 30.0),
      dx: 0.02,
      stress_type: "grid resolution near narrow kernel peak".to_owned(),
      analytical_caution: "potential overflow if domain widened".to_owned(),
    },

    // Case 9 — Tail precision stress
    GaussianTestCase {
      name: "tail_precision_stress".to_owned(),
      description: "Accurate accumulation of tiny tail mass".to_owned(),
      sigma_f: 3.0,
      sigma_g: 4.0,
      mu: 0.0,
      f_domain: (-36.0, 36.0),
      g_domain: (-48.0, 48.0),
      eval_domain: (-25.0, 25.0),
      dx: 0.01,
      stress_type: "truncation strategy, FFT padding".to_owned(),
      analytical_caution: "underflow to zero at extreme |x|".to_owned(),
    },

    // Case 10 — Tight truncation (intentional error exposure)
    GaussianTestCase {
      name: "tight_truncation".to_owned(),
      description: "Sensitivity to insufficient domain coverage".to_owned(),
      sigma_f: 1.0,
      sigma_g: 2.0,
      mu: 0.0,
      f_domain: (-3.0, 3.0),
      g_domain: (-6.0, 6.0),
      eval_domain: (-4.0, 4.0),
      dx: 0.01,
      stress_type: "boundary effects, wrap-around".to_owned(),
      analytical_caution: "none".to_owned(),
    },

    // Case 11 — Coarse grid
    GaussianTestCase {
      name: "coarse_grid".to_owned(),
      description: "Discretization error with sub-grid shift".to_owned(),
      sigma_f: 1.0,
      sigma_g: 1.0,
      mu: 0.5,
      f_domain: (-8.0, 8.0),
      g_domain: (-7.5, 8.5),
      eval_domain: (-6.0, 6.0),
      dx: 0.1,
      stress_type: "aliasing, grid spacing scaling".to_owned(),
      analytical_caution: "none".to_owned(),
    },

    // Case 12 — Very fine grid (reference-quality)
    GaussianTestCase {
      name: "fine_grid_reference".to_owned(),
      description: "High-accuracy baseline; should be closest to analytical target".to_owned(),
      sigma_f: 1.0,
      sigma_g: 1.0,
      mu: 0.0,
      f_domain: (-8.0, 8.0),
      g_domain: (-8.0, 8.0),
      eval_domain: (-6.0, 6.0),
      dx: 0.002,
      stress_type: "performance/memory, summation order".to_owned(),
      analytical_caution: "none".to_owned(),
    },

    // Case 13 — Wide dynamic range
    GaussianTestCase {
      name: "wide_dynamic_range".to_owned(),
      description: "Simultaneous handling of very sharp and very broad scales".to_owned(),
      sigma_f: 0.05,
      sigma_g: 7.5,
      mu: 6.0,
      f_domain: (-0.6, 0.6),
      g_domain: (-54.0, 66.0),
      eval_domain: (-10.0, 20.0),
      dx: 0.005,
      stress_type: "catastrophic cancellation, summation stability".to_owned(),
      analytical_caution: "under/overflow risks".to_owned(),
    },

    // Case 14 — Near-overflow guard (large |x| range)
    GaussianTestCase {
      name: "overflow_guard".to_owned(),
      description: "Numerical stability over very wide grids".to_owned(),
      sigma_f: 2.0,
      sigma_g: 2.5,
      mu: 0.0,
      f_domain: (-80.0, 80.0),
      g_domain: (-100.0, 100.0),
      eval_domain: (-50.0, 50.0),
      dx: 0.02,
      stress_type: "padding, exponent range, accumulation depth".to_owned(),
      analytical_caution: "underflow in far tails expected".to_owned(),
    },

    // Case 15 — Tiny shift, equal widths
    GaussianTestCase {
      name: "tiny_shift_equal_widths".to_owned(),
      description: "Sensitivity to small translations; symmetry properties".to_owned(),
      sigma_f: 1.25,
      sigma_g: 1.25,
      mu: 0.05,
      f_domain: (-10.0, 10.0),
      g_domain: (-9.95, 10.05),
      eval_domain: (-8.0, 8.0),
      dx: 0.01,
      stress_type: "interpolation/rounding at sub-grid peaks".to_owned(),
      analytical_caution: "none".to_owned(),
    },
  ]
}
