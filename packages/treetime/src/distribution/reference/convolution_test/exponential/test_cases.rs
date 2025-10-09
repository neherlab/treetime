use crate::distribution::reference::convolution_test::framework::TestCase;
use serde::{Deserialize, Serialize};

/// Test case parameters for Exponential convolution
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ExponentialTestCase {
  pub name: String,
  pub description: String,
  pub a: f64,                  // Decay rate for f(x) = exp(-ax) for x >= 0
  pub b: f64,                  // Decay rate for g(x) = exp(-bx) for x >= 0
  pub f_domain: (f64, f64),    // Should start at 0 for causal support
  pub g_domain: (f64, f64),    // Should start at 0 for causal support
  pub eval_domain: (f64, f64), // Should start at 0 for causal support
  pub dx: f64,
  pub stress_type: String,
  pub analytical_caution: String,
}

impl TestCase for ExponentialTestCase {
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

/// Create comprehensive Exponential test cases based on convolution-exponential.md reference
pub fn create_exponential_test_cases() -> Vec<ExponentialTestCase> {
  vec![
    // Case E1 — Baseline, distinct rates (matches pdf example)
    ExponentialTestCase {
      name: "baseline_distinct_rates".to_owned(),
      description: "General correctness, causal support, finite domains".to_owned(),
      a: 1.0,
      b: 2.0,
      f_domain: (0.0, 10.0),
      g_domain: (0.0, 7.0),
      eval_domain: (0.0, 17.0),
      dx: 0.01,
      stress_type: "none".to_owned(),
      analytical_caution: "none".to_owned(),
    },
    // Case E2 — Equal rates (limit case)
    ExponentialTestCase {
      name: "equal_rates_limit".to_owned(),
      description: "Correct handling of a=b using limit form x*exp(-ax)".to_owned(),
      a: 1.5,
      b: 1.5,
      f_domain: (0.0, 12.0),
      g_domain: (0.0, 12.0),
      eval_domain: (0.0, 24.0),
      dx: 0.01,
      stress_type: "kernel path selection, avoids division by zero".to_owned(),
      analytical_caution: "pdf formula undefined at a=b; use limit form".to_owned(),
    },
    // Case E3 — Near-equal rates (cancellation-sensitive)
    ExponentialTestCase {
      name: "near_equal_rates".to_owned(),
      description: "Numerical stability around a≈b; sensitivity to expm1-style evaluation".to_owned(),
      a: 1.0,
      b: 1.0 + 1e-6,
      f_domain: (0.0, 30.0),
      g_domain: (0.0, 30.0),
      eval_domain: (0.0, 60.0),
      dx: 0.002,
      stress_type: "discretization over long horizons, accumulation error".to_owned(),
      analytical_caution: "use stable expm1 form for (1-exp(-(a-b)x))/(a-b)".to_owned(),
    },
    // Case E4 — Fast + slow decay (wide dynamic range)
    ExponentialTestCase {
      name: "fast_slow_decay".to_owned(),
      description: "Interaction of very fast and very slow decays; long-tail coverage".to_owned(),
      a: 10.0,
      b: 0.1,
      f_domain: (0.0, 2.0),
      g_domain: (0.0, 200.0),
      eval_domain: (0.0, 202.0),
      dx: 0.001,
      stress_type: "truncation strategy for long g tail, FFT zero-padding".to_owned(),
      analytical_caution: "underflow in far tail contributions acceptable".to_owned(),
    },
    // Case E5 — Slow + fast decay (swap roles)
    ExponentialTestCase {
      name: "slow_fast_decay".to_owned(),
      description: "Symmetry under role swap; causal truncation interaction".to_owned(),
      a: 0.1,
      b: 10.0,
      f_domain: (0.0, 200.0),
      g_domain: (0.0, 2.0),
      eval_domain: (0.0, 202.0),
      dx: 0.001,
      stress_type: "padding and index handling at short g support".to_owned(),
      analytical_caution: "none".to_owned(),
    },
    // Case E6 — Moderate rates, coarse grid
    ExponentialTestCase {
      name: "moderate_coarse_grid".to_owned(),
      description: "Discretization error on coarse grids; step size sensitivity".to_owned(),
      a: 1.0,
      b: 0.8,
      f_domain: (0.0, 20.0),
      g_domain: (0.0, 20.0),
      eval_domain: (0.0, 40.0),
      dx: 0.1,
      stress_type: "aliasing/accuracy vs Δx, Δx scaling factor".to_owned(),
      analytical_caution: "none".to_owned(),
    },
    // Case E7 — Fine grid reference
    ExponentialTestCase {
      name: "fine_grid_reference".to_owned(),
      description: "High-accuracy baseline for regression; compare with E6".to_owned(),
      a: 1.0,
      b: 0.8,
      f_domain: (0.0, 20.0),
      g_domain: (0.0, 20.0),
      eval_domain: (0.0, 40.0),
      dx: 0.002,
      stress_type: "performance/memory, summation order".to_owned(),
      analytical_caution: "none".to_owned(),
    },
    // Case E8 — Long horizon tail-precision
    ExponentialTestCase {
      name: "long_horizon_tail_precision".to_owned(),
      description: "Accumulation of small tail mass over long support".to_owned(),
      a: 0.05,
      b: 0.04,
      f_domain: (0.0, 600.0),
      g_domain: (0.0, 750.0),
      eval_domain: (0.0, 1350.0),
      dx: 0.01,
      stress_type: "FFT padding to avoid circular wrap, cumulative round-off".to_owned(),
      analytical_caution: "underflow in far tails expected".to_owned(),
    },
    // Case E9 — Tight truncation (intentional error exposure)
    ExponentialTestCase {
      name: "tight_truncation".to_owned(),
      description: "Sensitivity to insufficient support; boundary effects".to_owned(),
      a: 1.0,
      b: 2.0,
      f_domain: (0.0, 3.0),
      g_domain: (0.0, 2.0),
      eval_domain: (0.0, 5.0),
      dx: 0.01,
      stress_type: "wrap-around artifacts if padding insufficient".to_owned(),
      analytical_caution: "none".to_owned(),
    },
    // Case E10 — Large x range, near-underflow guard
    ExponentialTestCase {
      name: "large_range_underflow_guard".to_owned(),
      description: "Stability across wide index ranges; causal support check".to_owned(),
      a: 0.5,
      b: 0.7,
      f_domain: (0.0, 400.0),
      g_domain: (0.0, 400.0),
      eval_domain: (0.0, 800.0),
      dx: 0.02,
      stress_type: "exponent range, padding, accumulation depth".to_owned(),
      analytical_caution: "far tails may underflow to zero; acceptable".to_owned(),
    },
    // Case E11 — Extreme near-equality (machine epsilon scale)
    ExponentialTestCase {
      name: "extreme_near_equality".to_owned(),
      description: "Robustness using expm1 and stable division".to_owned(),
      a: 1.0,
      b: 1.0 + 1e-12,
      f_domain: (0.0, 80.0),
      g_domain: (0.0, 80.0),
      eval_domain: (0.0, 160.0),
      dx: 0.005,
      stress_type: "long support convolution accuracy, summation stability".to_owned(),
      analytical_caution: "direct formula ill-conditioned; use stable evaluation".to_owned(),
    },
    // Case E12 — Very fast decays, small domains
    ExponentialTestCase {
      name: "very_fast_decays".to_owned(),
      description: "Accuracy on sharply localized causal signals".to_owned(),
      a: 25.0,
      b: 30.0,
      f_domain: (0.0, 0.5),
      g_domain: (0.0, 0.5),
      eval_domain: (0.0, 1.0),
      dx: 0.0005,
      stress_type: "grid resolution and Δx scaling".to_owned(),
      analytical_caution: "none".to_owned(),
    },
  ]
}
