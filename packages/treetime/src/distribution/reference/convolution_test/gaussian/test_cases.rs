use crate::distribution::reference::convolution_test::framework::TestCase;
use serde::{Deserialize, Serialize};

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

/// Create comprehensive Gaussian test cases based on convolution-gaussian.md reference
pub fn create_gaussian_test_cases() -> Vec<GaussianTestCase> {
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
