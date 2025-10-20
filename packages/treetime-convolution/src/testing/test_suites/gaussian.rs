#![allow(clippy::many_single_char_names)]
use crate::testing::framework::test_case::TestCase;
use crate::testing::test_suites::test_suites::TestSuite;
use eyre::Report;
use ndarray::Array1;
use serde::{Deserialize, Serialize};
use std::f64::consts::PI;

#[derive(Default)]
pub struct GaussianTestSuite;

impl TestSuite for GaussianTestSuite {
  type TestCase = GaussianTestCase;

  fn test_suite_name(&self) -> &'static str {
    "gaussian"
  }

  fn create_f(&self, test_case: &Self::TestCase, grid: &Array1<f64>) -> Result<Array1<f64>, Report> {
    let sigma_f = test_case.sigma_f;
    Ok(grid.mapv(|x| (-(0.5 * (x / sigma_f).powi(2))).exp() / (sigma_f * (2.0 * PI).sqrt())))
  }

  fn create_g(&self, test_case: &Self::TestCase, grid: &Array1<f64>) -> Result<Array1<f64>, Report> {
    let sigma_g = test_case.sigma_g;
    let mu = test_case.mu;
    Ok(grid.mapv(|x| (-(0.5 * ((x - mu) / sigma_g).powi(2))).exp() / (sigma_g * (2.0 * PI).sqrt())))
  }

  fn analytical_convolution(&self, test_case: &Self::TestCase, eval_grid: &Array1<f64>) -> Result<Array1<f64>, Report> {
    let sigma_f = test_case.sigma_f;
    let sigma_g = test_case.sigma_g;
    let mu = test_case.mu;
    let variance_sum = sigma_f.powi(2) + sigma_g.powi(2);
    Ok(eval_grid.mapv(|x| (-(0.5 * (x - mu).powi(2) / variance_sum)).exp() / (2.0 * PI * variance_sum).sqrt()))
  }

  fn create_test_cases(&self) -> Vec<Self::TestCase> {
    vec![
      GaussianTestCase {
        name: "python_notebook_case_1".to_owned(),
        description: "Parameters from conv_gauss_py.ipynb: σ_f=1, σ_g=1, μ=1".to_owned(),
        stress_type: "reference implementation validation".to_owned(),
        analytical_caution: "none".to_owned(),
        sigma_f: 1.0,
        sigma_g: 1.0,
        mu: 1.0,
        input_grid_domain: (-5.0, 5.0),
        input_grid_n_points: 101,
        output_grid_domain: (-9.0, 11.0),
        output_grid_n_points: 101,
      },
      GaussianTestCase {
        name: "python_notebook_case_2".to_owned(),
        description: "Parameters from conv_gauss_py.ipynb: σ_f=2, σ_g=1, μ=1".to_owned(),
        stress_type: "reference implementation validation".to_owned(),
        analytical_caution: "none".to_owned(),
        sigma_f: 2.0,
        sigma_g: 1.0,
        mu: 1.0,
        input_grid_domain: (-5.0, 5.0),
        input_grid_n_points: 101,
        output_grid_domain: (-9.0, 11.0),
        output_grid_n_points: 101,
      },
      GaussianTestCase {
        name: "tight_truncation".to_owned(),
        description: "Tight truncation highlighting sensitivity to insufficient domain coverage.".to_owned(),
        stress_type: "boundary effects, wrap-around".to_owned(),
        analytical_caution: "none".to_owned(),
        sigma_f: 1.0,
        sigma_g: 2.0,
        mu: 0.0,
        input_grid_domain: (-6.0, 6.0),
        input_grid_n_points: 1201,
        output_grid_domain: (-8.0, 8.0),
        output_grid_n_points: 1601,
      },
      GaussianTestCase {
        name: "coarse_grid".to_owned(),
        description: "Coarse grid coverage emphasizing discretization error with sub-grid shift.".to_owned(),
        stress_type: "aliasing, grid spacing scaling".to_owned(),
        analytical_caution: "none".to_owned(),
        sigma_f: 1.0,
        sigma_g: 1.0,
        mu: 0.5,
        input_grid_domain: (-8.0, 8.0),
        input_grid_n_points: 161,
        output_grid_domain: (-10.0, 10.0),
        output_grid_n_points: 201,
      },
      GaussianTestCase {
        name: "baseline_moderate".to_owned(),
        description: "Baseline scenario with moderate widths; general correctness with asymmetric widths and shift."
          .to_owned(),
        stress_type: "none".to_owned(),
        analytical_caution: "none".to_owned(),
        sigma_f: 1.0,
        sigma_g: 2.0,
        mu: 1.0,
        input_grid_domain: (-12.0, 14.0),
        input_grid_n_points: 2601,
        output_grid_domain: (-10.0, 12.0),
        output_grid_n_points: 2201,
      },
      GaussianTestCase {
        name: "centered_unequal_widths".to_owned(),
        description: "Centered configuration with unequal widths; width mixing and peak alignment at x=0.".to_owned(),
        stress_type: "resolution vs narrow kernel".to_owned(),
        analytical_caution: "none".to_owned(),
        sigma_f: 1.5,
        sigma_g: 0.5,
        mu: 0.0,
        input_grid_domain: (-12.0, 12.0),
        input_grid_n_points: 2401,
        output_grid_domain: (-8.0, 8.0),
        output_grid_n_points: 1601,
      },
      GaussianTestCase {
        name: "delta_like_small_shift".to_owned(),
        description: "Sharp g with small shift; sensitivity to sub-grid shifts and near-delta kernel behavior."
          .to_owned(),
        stress_type: "discretization error, normalization".to_owned(),
        analytical_caution: "potential overflow for large |x|/σ_g".to_owned(),
        sigma_f: 1.0,
        sigma_g: 0.05,
        mu: 0.2,
        input_grid_domain: (-8.0, 8.0),
        input_grid_n_points: 8001,
        output_grid_domain: (-4.0, 4.0),
        output_grid_n_points: 4001,
      },
      GaussianTestCase {
        name: "wide_kernel_negative_shift".to_owned(),
        description: "Very wide g with negative shift; robustness to ultra-broad kernel tails and translation."
          .to_owned(),
        stress_type: "tail truncation, numeric underflow, FFT padding".to_owned(),
        analytical_caution: "underflow in far tails expected".to_owned(),
        sigma_f: 0.5,
        sigma_g: 5.0,
        mu: -1.0,
        input_grid_domain: (-41.0, 39.0),
        input_grid_n_points: 4001,
        output_grid_domain: (-20.0, 20.0),
        output_grid_n_points: 2001,
      },
      GaussianTestCase {
        name: "large_positive_shift".to_owned(),
        description: "Large positive shift validating translation invariance for sizeable μ.".to_owned(),
        stress_type: "circular convolution artifacts, insufficient padding".to_owned(),
        analytical_caution: "none".to_owned(),
        sigma_f: 1.0,
        sigma_g: 1.2,
        mu: 10.0,
        input_grid_domain: (-10.0, 22.0),
        input_grid_n_points: 3201,
        output_grid_domain: (0.0, 20.0),
        output_grid_n_points: 2001,
      },
      GaussianTestCase {
        name: "large_negative_shift".to_owned(),
        description: "Large negative shift exercising translation invariance for negative μ.".to_owned(),
        stress_type: "padding and index handling".to_owned(),
        analytical_caution: "none".to_owned(),
        sigma_f: 2.0,
        sigma_g: 1.0,
        mu: -7.0,
        input_grid_domain: (-16.0, 16.0),
        input_grid_n_points: 3201,
        output_grid_domain: (-20.0, 0.0),
        output_grid_n_points: 2001,
      },
      GaussianTestCase {
        name: "extreme_ratio_wide_g".to_owned(),
        description:
          "Extreme width ratio with g much wider than f; convolution where g is nearly constant over f's support."
            .to_owned(),
        stress_type: "tail coverage, numerical cancellation".to_owned(),
        analytical_caution: "underflow in extreme tails".to_owned(),
        sigma_f: 0.1,
        sigma_g: 10.0,
        mu: 3.0,
        input_grid_domain: (-77.0, 83.0),
        input_grid_n_points: 16001,
        output_grid_domain: (-5.0, 15.0),
        output_grid_n_points: 2001,
      },
      GaussianTestCase {
        name: "extreme_ratio_wide_f".to_owned(),
        description: "Extreme width ratio with f much wider than g; near-delta g translating f.".to_owned(),
        stress_type: "grid resolution near narrow kernel peak".to_owned(),
        analytical_caution: "potential overflow if domain widened".to_owned(),
        sigma_f: 8.0,
        sigma_g: 0.2,
        mu: -2.5,
        input_grid_domain: (-64.0, 64.0),
        input_grid_n_points: 6401,
        output_grid_domain: (-30.0, 30.0),
        output_grid_n_points: 3001,
      },
      GaussianTestCase {
        name: "tail_precision_stress".to_owned(),
        description: "Tail precision stress with accurate accumulation of tiny tail mass.".to_owned(),
        stress_type: "truncation strategy, FFT padding".to_owned(),
        analytical_caution: "underflow to zero at extreme |x|".to_owned(),
        sigma_f: 3.0,
        sigma_g: 4.0,
        mu: 0.0,
        input_grid_domain: (-48.0, 48.0),
        input_grid_n_points: 9601,
        output_grid_domain: (-25.0, 25.0),
        output_grid_n_points: 5001,
      },
      GaussianTestCase {
        name: "fine_grid_reference".to_owned(),
        description: "Reference-quality very fine grid; high-accuracy baseline closest to analytical target."
          .to_owned(),
        stress_type: "performance/memory, summation order".to_owned(),
        analytical_caution: "none".to_owned(),
        sigma_f: 1.0,
        sigma_g: 1.0,
        mu: 0.0,
        input_grid_domain: (-8.0, 8.0),
        input_grid_n_points: 8001,
        output_grid_domain: (-6.0, 6.0),
        output_grid_n_points: 6001,
      },
      GaussianTestCase {
        name: "wide_dynamic_range".to_owned(),
        description: "Wide dynamic range with simultaneous handling of very sharp and very broad scales.".to_owned(),
        stress_type: "catastrophic cancellation, summation stability".to_owned(),
        analytical_caution: "under/overflow risks".to_owned(),
        sigma_f: 0.05,
        sigma_g: 7.5,
        mu: 6.0,
        input_grid_domain: (-54.0, 66.0),
        input_grid_n_points: 24001,
        output_grid_domain: (-10.0, 20.0),
        output_grid_n_points: 6001,
      },
      GaussianTestCase {
        name: "overflow_guard".to_owned(),
        description: "Near-overflow guard ensuring numerical stability over very wide grids.".to_owned(),
        stress_type: "padding, exponent range, accumulation depth".to_owned(),
        analytical_caution: "underflow in far tails expected".to_owned(),
        sigma_f: 2.0,
        sigma_g: 2.5,
        mu: 0.0,
        input_grid_domain: (-100.0, 100.0),
        input_grid_n_points: 10001,
        output_grid_domain: (-50.0, 50.0),
        output_grid_n_points: 5001,
      },
      GaussianTestCase {
        name: "tiny_shift_equal_widths".to_owned(),
        description: "Tiny shift with equal widths; sensitivity to small translations and symmetry properties."
          .to_owned(),
        stress_type: "interpolation/rounding at sub-grid peaks".to_owned(),
        analytical_caution: "none".to_owned(),
        sigma_f: 1.25,
        sigma_g: 1.25,
        mu: 0.05,
        input_grid_domain: (-10.0, 10.05),
        input_grid_n_points: 2006,
        output_grid_domain: (-8.0, 8.0),
        output_grid_n_points: 1601,
      },
    ]
  }
}

/// Test case parameters for Gaussian convolution
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GaussianTestCase {
  pub name: String,
  pub description: String,
  pub stress_type: String,
  pub analytical_caution: String,
  pub sigma_f: f64,
  pub sigma_g: f64,
  pub mu: f64,
  pub input_grid_domain: (f64, f64),
  pub input_grid_n_points: usize,
  pub output_grid_domain: (f64, f64),
  pub output_grid_n_points: usize,
}

impl TestCase for GaussianTestCase {
  fn name(&self) -> &str {
    &self.name
  }

  fn description(&self) -> &str {
    &self.description
  }

  fn stress_type(&self) -> &str {
    &self.stress_type
  }

  fn analytical_caution(&self) -> &str {
    &self.analytical_caution
  }

  fn input_grid_domain(&self) -> (f64, f64) {
    self.input_grid_domain
  }

  fn input_grid_n_points(&self) -> usize {
    self.input_grid_n_points
  }

  fn output_grid_domain(&self) -> (f64, f64) {
    self.output_grid_domain
  }

  fn output_grid_n_points(&self) -> usize {
    self.output_grid_n_points
  }
}
