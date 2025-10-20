#![allow(clippy::many_single_char_names)]
use crate::distribution::reference::convolution_test::framework::test_case::TestCase;
use crate::distribution::reference::convolution_test::traits::ConvInput;
use crate::distribution::reference::grid_fn::GridFn;
use eyre::Report;
use ndarray::Array1;
use serde::{Deserialize, Serialize};
use std::f64::consts::PI;

pub struct GaussianConvInput {
  test_cases: Vec<GaussianTestCase>,
}

impl ConvInput for GaussianConvInput {
  type TestCase = GaussianTestCase;

  fn function_type(&self) -> &'static str {
    "gaussian"
  }

  fn new_with_cases(test_cases: Vec<GaussianTestCase>) -> Self {
    Self { test_cases }
  }

  fn create_f(
    &self,
    &GaussianTestCase {
      sigma_f, f_domain, dx, ..
    }: &Self::TestCase,
  ) -> Result<GridFn, Report> {
    // Generate Gaussian function f(x) = (1/√(2π σ_f²)) * exp(-x²/(2σ_f²))
    GridFn::from_grid(f_domain, dx, |x| {
      (-(0.5 * (x / sigma_f).powi(2))).exp() / (sigma_f * (2.0 * PI).sqrt())
    })
  }

  fn create_g(
    &self,
    &GaussianTestCase {
      sigma_g,
      mu,
      g_domain,
      dx,
      ..
    }: &Self::TestCase,
  ) -> Result<GridFn, Report> {
    // Generate Gaussian function g(x) = (1/√(2π σ_g²)) * exp(-(x-μ)²/(2σ_g²))
    GridFn::from_grid(g_domain, dx, |x| {
      (-(0.5 * ((x - mu) / sigma_g).powi(2))).exp() / (sigma_g * (2.0 * PI).sqrt())
    })
  }

  fn eval_domain(&self, test_case: &Self::TestCase) -> (f64, f64) {
    test_case.eval_domain
  }

  fn analytical_convolution(
    &self,
    &GaussianTestCase {
      sigma_f,
      sigma_g,
      mu,
      eval_domain,
      dx,
      ..
    }: &Self::TestCase,
    _eval_grid: &Array1<f64>,
  ) -> Result<GridFn, Report> {
    // Analytical convolution of two Gaussian functions
    //
    // f(x) = (1/√(2π σ_f²)) * exp(-x²/(2σ_f²))
    // g(x) = (1/√(2π σ_g²)) * exp(-(x-μ)²/(2σ_g²))
    //
    // Result: (1/√(2π(σ_f² + σ_g²))) * exp(-(x-μ)²/(2(σ_f² + σ_g²)))
    let variance_sum = sigma_f.powi(2) + sigma_g.powi(2);
    GridFn::from_grid(eval_domain, dx, |x| {
      (-(0.5 * (x - mu).powi(2) / variance_sum)).exp() / (2.0 * PI * variance_sum).sqrt()
    })
  }

  fn test_cases(&self) -> &[Self::TestCase] {
    &self.test_cases
  }

  fn create_test_cases() -> Vec<GaussianTestCase> {
    vec![
      // GaussianTestCase {
      //   name: "baseline_moderate".to_owned(),
      //   description: "Baseline scenario with moderate widths; general correctness with asymmetric widths and shift."
      //     .to_owned(),
      //   stress_type: "none".to_owned(),
      //   analytical_caution: "none".to_owned(),
      //   sigma_f: 1.0,
      //   sigma_g: 2.0,
      //   mu: 1.0,
      //   f_domain: (-8.0, 8.0),
      //   g_domain: (-12.0, 14.0),
      //   eval_domain: (-10.0, 12.0),
      //   dx: 0.01,
      // },
      // GaussianTestCase {
      //   name: "centered_unequal_widths".to_owned(),
      //   description: "Centered configuration with unequal widths; width mixing and peak alignment at x=0.".to_owned(),
      //   stress_type: "resolution vs narrow kernel".to_owned(),
      //   analytical_caution: "none".to_owned(),
      //   sigma_f: 1.5,
      //   sigma_g: 0.5,
      //   mu: 0.0,
      //   f_domain: (-12.0, 12.0),
      //   g_domain: (-4.0, 4.0),
      //   eval_domain: (-8.0, 8.0),
      //   dx: 0.01,
      // },
      // GaussianTestCase {
      //   name: "delta_like_small_shift".to_owned(),
      //   description: "Sharp g with small shift; sensitivity to sub-grid shifts and near-delta kernel behavior."
      //     .to_owned(),
      //   stress_type: "discretization error, normalization".to_owned(),
      //   analytical_caution: "potential overflow for large |x|/σ_g".to_owned(),
      //   sigma_f: 1.0,
      //   sigma_g: 0.05,
      //   mu: 0.2,
      //   f_domain: (-8.0, 8.0),
      //   g_domain: (-0.4, 1.0),
      //   eval_domain: (-4.0, 4.0),
      //   dx: 0.002,
      // },
      // GaussianTestCase {
      //   name: "wide_kernel_negative_shift".to_owned(),
      //   description: "Very wide g with negative shift; robustness to ultra-broad kernel tails and translation."
      //     .to_owned(),
      //   stress_type: "tail truncation, numeric underflow, FFT padding".to_owned(),
      //   analytical_caution: "underflow in far tails expected".to_owned(),
      //   sigma_f: 0.5,
      //   sigma_g: 5.0,
      //   mu: -1.0,
      //   f_domain: (-6.0, 6.0),
      //   g_domain: (-41.0, 39.0),
      //   eval_domain: (-20.0, 20.0),
      //   dx: 0.02,
      // },
      // GaussianTestCase {
      //   name: "large_positive_shift".to_owned(),
      //   description: "Large positive shift validating translation invariance for sizeable μ.".to_owned(),
      //   stress_type: "circular convolution artifacts, insufficient padding".to_owned(),
      //   analytical_caution: "none".to_owned(),
      //   sigma_f: 1.0,
      //   sigma_g: 1.2,
      //   mu: 10.0,
      //   f_domain: (-10.0, 10.0),
      //   g_domain: (-2.0, 22.0),
      //   eval_domain: (0.0, 20.0),
      //   dx: 0.01,
      // },
      // GaussianTestCase {
      //   name: "large_negative_shift".to_owned(),
      //   description: "Large negative shift exercising translation invariance for negative μ.".to_owned(),
      //   stress_type: "padding and index handling".to_owned(),
      //   analytical_caution: "none".to_owned(),
      //   sigma_f: 2.0,
      //   sigma_g: 1.0,
      //   mu: -7.0,
      //   f_domain: (-16.0, 16.0),
      //   g_domain: (-15.0, 1.0),
      //   eval_domain: (-20.0, 0.0),
      //   dx: 0.01,
      // },
      // GaussianTestCase {
      //   name: "extreme_ratio_wide_g".to_owned(),
      //   description:
      //     "Extreme width ratio with g much wider than f; convolution where g is nearly constant over f's support."
      //       .to_owned(),
      //   stress_type: "tail coverage, numerical cancellation".to_owned(),
      //   analytical_caution: "underflow in extreme tails".to_owned(),
      //   sigma_f: 0.1,
      //   sigma_g: 10.0,
      //   mu: 3.0,
      //   f_domain: (-1.0, 1.0),
      //   g_domain: (-77.0, 83.0),
      //   eval_domain: (-5.0, 15.0),
      //   dx: 0.01,
      // },
      // GaussianTestCase {
      //   name: "extreme_ratio_wide_f".to_owned(),
      //   description: "Extreme width ratio with f much wider than g; near-delta g translating f.".to_owned(),
      //   stress_type: "grid resolution near narrow kernel peak".to_owned(),
      //   analytical_caution: "potential overflow if domain widened".to_owned(),
      //   sigma_f: 8.0,
      //   sigma_g: 0.2,
      //   mu: -2.5,
      //   f_domain: (-64.0, 64.0),
      //   g_domain: (-4.1, -0.9),
      //   eval_domain: (-30.0, 30.0),
      //   dx: 0.02,
      // },
      // GaussianTestCase {
      //   name: "tail_precision_stress".to_owned(),
      //   description: "Tail precision stress with accurate accumulation of tiny tail mass.".to_owned(),
      //   stress_type: "truncation strategy, FFT padding".to_owned(),
      //   analytical_caution: "underflow to zero at extreme |x|".to_owned(),
      //   sigma_f: 3.0,
      //   sigma_g: 4.0,
      //   mu: 0.0,
      //   f_domain: (-36.0, 36.0),
      //   g_domain: (-48.0, 48.0),
      //   eval_domain: (-25.0, 25.0),
      //   dx: 0.01,
      // },
      GaussianTestCase {
        name: "tight_truncation".to_owned(),
        description: "Tight truncation highlighting sensitivity to insufficient domain coverage.".to_owned(),
        stress_type: "boundary effects, wrap-around".to_owned(),
        analytical_caution: "none".to_owned(),
        sigma_f: 1.0,
        sigma_g: 2.0,
        mu: 0.0,
        f_domain: (-3.0, 3.0),
        g_domain: (-6.0, 6.0),
        eval_domain: (-4.0, 4.0),
        dx: 0.01,
      },
      GaussianTestCase {
        name: "coarse_grid".to_owned(),
        description: "Coarse grid coverage emphasizing discretization error with sub-grid shift.".to_owned(),
        stress_type: "aliasing, grid spacing scaling".to_owned(),
        analytical_caution: "none".to_owned(),
        sigma_f: 1.0,
        sigma_g: 1.0,
        mu: 0.5,
        f_domain: (-8.0, 8.0),
        g_domain: (-7.5, 8.5),
        eval_domain: (-6.0, 6.0),
        dx: 0.1,
      },
      // GaussianTestCase {
      //   name: "fine_grid_reference".to_owned(),
      //   description: "Reference-quality very fine grid; high-accuracy baseline closest to analytical target.".to_owned(),
      //   stress_type: "performance/memory, summation order".to_owned(),
      //   analytical_caution: "none".to_owned(),
      //   sigma_f: 1.0,
      //   sigma_g: 1.0,
      //   mu: 0.0,
      //   f_domain: (-8.0, 8.0),
      //   g_domain: (-8.0, 8.0),
      //   eval_domain: (-6.0, 6.0),
      //   dx: 0.002,
      // },
      // GaussianTestCase {
      //   name: "wide_dynamic_range".to_owned(),
      //   description: "Wide dynamic range with simultaneous handling of very sharp and very broad scales.".to_owned(),
      //   stress_type: "catastrophic cancellation, summation stability".to_owned(),
      //   analytical_caution: "under/overflow risks".to_owned(),
      //   sigma_f: 0.05,
      //   sigma_g: 7.5,
      //   mu: 6.0,
      //   f_domain: (-0.6, 0.6),
      //   g_domain: (-54.0, 66.0),
      //   eval_domain: (-10.0, 20.0),
      //   dx: 0.005,
      // },
      // GaussianTestCase {
      //   name: "overflow_guard".to_owned(),
      //   description: "Near-overflow guard ensuring numerical stability over very wide grids.".to_owned(),
      //   stress_type: "padding, exponent range, accumulation depth".to_owned(),
      //   analytical_caution: "underflow in far tails expected".to_owned(),
      //   sigma_f: 2.0,
      //   sigma_g: 2.5,
      //   mu: 0.0,
      //   f_domain: (-80.0, 80.0),
      //   g_domain: (-100.0, 100.0),
      //   eval_domain: (-50.0, 50.0),
      //   dx: 0.02,
      // },
      // GaussianTestCase {
      //   name: "tiny_shift_equal_widths".to_owned(),
      //   description: "Tiny shift with equal widths; sensitivity to small translations and symmetry properties."
      //     .to_owned(),
      //   stress_type: "interpolation/rounding at sub-grid peaks".to_owned(),
      //   analytical_caution: "none".to_owned(),
      //   sigma_f: 1.25,
      //   sigma_g: 1.25,
      //   mu: 0.05,
      //   f_domain: (-10.0, 10.0),
      //   g_domain: (-9.95, 10.05),
      //   eval_domain: (-8.0, 8.0),
      //   dx: 0.01,
      // },
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
  pub f_domain: (f64, f64),
  pub g_domain: (f64, f64),
  pub eval_domain: (f64, f64),
  pub dx: f64,
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

  fn dx(&self) -> f64 {
    self.dx
  }
}
