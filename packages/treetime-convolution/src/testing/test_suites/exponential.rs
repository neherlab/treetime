#![allow(clippy::many_single_char_names)]
use crate::testing::framework::test_case::TestCase;
use crate::testing::test_suites::test_suites::TestSuite;
use eyre::Report;
use ndarray::Array1;
use serde::{Deserialize, Serialize};

#[derive(Default)]
pub struct ExponentialTestSuite;

impl TestSuite for ExponentialTestSuite {
  type TestCase = ExponentialTestCase;

  fn test_suite_name(&self) -> &'static str {
    "exponential"
  }

  fn create_f(&self, test_case: &Self::TestCase, grid: &Array1<f64>) -> Result<Array1<f64>, Report> {
    let a = test_case.a;
    Ok(grid.mapv(|x| if x >= 0.0 { a * (-a * x).exp() } else { 0.0 }))
  }

  fn create_g(&self, test_case: &Self::TestCase, grid: &Array1<f64>) -> Result<Array1<f64>, Report> {
    let b = test_case.b;
    Ok(grid.mapv(|x| if x >= 0.0 { b * (-b * x).exp() } else { 0.0 }))
  }

  fn analytical_convolution(&self, test_case: &Self::TestCase, eval_grid: &Array1<f64>) -> Result<Array1<f64>, Report> {
    let a = test_case.a;
    let b = test_case.b;
    Ok(eval_grid.mapv(|x| {
      if x < 0.0 {
        0.0
      } else if (a - b).abs() < 1e-15 {
        a * b * x * (-a * x).exp()
      } else {
        (a * b) / (a - b) * (1.0 - (-(a - b) * x).exp()) * (-b * x).exp()
      }
    }))
  }

  fn create_test_cases(&self) -> Vec<Self::TestCase> {
    vec![
      ExponentialTestCase {
        name: "python_notebook_case".to_owned(),
        description: "Parameters from conv_exp_py.ipynb: a=1, b=2".to_owned(),
        stress_type: "reference implementation validation".to_owned(),
        analytical_caution: "none".to_owned(),
        slowness: 0.0,
        a: 1.0,
        b: 2.0,
        input_grid_domain: (-1.0, 10.0),
        input_grid_n_points: 1101,
      },
      ExponentialTestCase {
        name: "moderate_coarse_grid".to_owned(),
        description: "Moderate rates on a coarse grid. Discretization error on coarse grids; step size sensitivity."
          .to_owned(),
        stress_type: "aliasing/accuracy vs Δx, Δx scaling factor".to_owned(),
        analytical_caution: "none".to_owned(),
        slowness: 0.0,
        a: 1.0,
        b: 0.8,
        input_grid_domain: (0.0, 20.0),
        input_grid_n_points: 201,
      },
      ExponentialTestCase {
        name: "tight_truncation".to_owned(),
        description:
          "Tight truncation to expose intentional errors. Sensitivity to insufficient support; boundary effects."
            .to_owned(),
        stress_type: "wrap-around artifacts if padding insufficient".to_owned(),
        analytical_caution: "none".to_owned(),
        slowness: 0.0,
        a: 1.0,
        b: 2.0,
        input_grid_domain: (0.0, 5.0),
        input_grid_n_points: 501,
      },
      ExponentialTestCase {
        name: "baseline_distinct_rates".to_owned(),
        description: "Baseline with distinct rates (matches pdf example). General correctness, causal support, finite domains."
          .to_owned(),
        stress_type: "none".to_owned(),
        analytical_caution: "none".to_owned(),
        slowness: 0.0,
        a: 1.0,
        b: 2.0,
        input_grid_domain: (0.0, 10.0),
        input_grid_n_points: 1001,
      },
      ExponentialTestCase {
        name: "equal_rates_limit".to_owned(),
        description: "Equal-rate limit scenario. Correct handling of a=b using limit form x*exp(-ax).".to_owned(),
        stress_type: "kernel path selection, avoids division by zero".to_owned(),
        analytical_caution: "pdf formula undefined at a=b; use limit form".to_owned(),
        slowness: 0.0,
        a: 1.5,
        b: 1.5,
        input_grid_domain: (0.0, 12.0),
        input_grid_n_points: 1201,
      },
      ExponentialTestCase {
        name: "near_equal_rates".to_owned(),
        description: "Near-equal rates, cancellation sensitive. Numerical stability around a≈b; sensitivity to expm1-style evaluation."
          .to_owned(),
        stress_type: "discretization over long horizons, accumulation error".to_owned(),
        analytical_caution: "use stable expm1 form for (1-exp(-(a-b)x))/(a-b)".to_owned(),
        slowness: 0.6,
        a: 1.0,
        b: 1.0 + 1e-6,
        input_grid_domain: (0.0, 30.0),
        input_grid_n_points: 15001,
      },
      ExponentialTestCase {
        name: "fast_slow_decay".to_owned(),
        description: "Fast plus slow decay with wide dynamic range. Interaction of very fast and very slow decays; long-tail coverage."
          .to_owned(),
        stress_type: "truncation strategy for long g tail, FFT zero-padding".to_owned(),
        analytical_caution: "underflow in far tail contributions acceptable".to_owned(),
        slowness: 1.0,
        a: 10.0,
        b: 0.1,
        input_grid_domain: (0.0, 200.0),
        input_grid_n_points: 200001,
      },
      ExponentialTestCase {
        name: "slow_fast_decay".to_owned(),
        description: "Slow plus fast decay swapping roles. Symmetry under role swap; causal truncation interaction."
          .to_owned(),
        stress_type: "padding and index handling at short g support".to_owned(),
        analytical_caution: "none".to_owned(),
        slowness: 1.0,
        a: 0.1,
        b: 10.0,
        input_grid_domain: (0.0, 200.0),
        input_grid_n_points: 200001,
      },
      ExponentialTestCase {
        name: "fine_grid_reference".to_owned(),
        description: "High-accuracy fine grid reference baseline. Useful for regression against the coarse grid scenario."
          .to_owned(),
        stress_type: "performance/memory, summation order".to_owned(),
        analytical_caution: "none".to_owned(),
        slowness: 0.9,
        a: 1.0,
        b: 0.8,
        input_grid_domain: (0.0, 20.0),
        input_grid_n_points: 10001,
      },
      ExponentialTestCase {
        name: "long_horizon_tail_precision".to_owned(),
        description: "Long horizon tail precision focus. Accumulation of small tail mass over long support.".to_owned(),
        stress_type: "FFT padding to avoid circular wrap, cumulative round-off".to_owned(),
        analytical_caution: "underflow in far tails expected".to_owned(),
        slowness: 0.95,
        a: 0.05,
        b: 0.04,
        input_grid_domain: (0.0, 750.0),
        input_grid_n_points: 75001,
      },
      ExponentialTestCase {
        name: "large_range_underflow_guard".to_owned(),
        description: "Large x range with near-underflow guard. Stability across wide index ranges; causal support check."
          .to_owned(),
        stress_type: "exponent range, padding, accumulation depth".to_owned(),
        analytical_caution: "far tails may underflow to zero; acceptable".to_owned(),
        slowness: 0.7,
        a: 0.5,
        b: 0.7,
        input_grid_domain: (0.0, 400.0),
        input_grid_n_points: 20001,
      },
      ExponentialTestCase {
        name: "extreme_near_equality".to_owned(),
        description: "Extreme near equality at machine epsilon scale. Robustness using expm1 and stable division."
          .to_owned(),
        stress_type: "long support convolution accuracy, summation stability".to_owned(),
        analytical_caution: "direct formula ill-conditioned; use stable evaluation".to_owned(),
        slowness: 0.65,
        a: 1.0,
        b: 1.0 + 1e-12,
        input_grid_domain: (0.0, 80.0),
        input_grid_n_points: 16001,
        
        
      },
      ExponentialTestCase {
        name: "very_fast_decays".to_owned(),
        description: "Very fast decays on small domains. Accuracy on sharply localized causal signals.".to_owned(),
        stress_type: "grid resolution and Δx scaling".to_owned(),
        analytical_caution: "none".to_owned(),
        slowness: 0.0,
        a: 25.0,
        b: 30.0,
        input_grid_domain: (0.0, 0.5),
        input_grid_n_points: 1001,
        
        
      },
    ]
  }
}

/// Test case parameters for Exponential convolution
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ExponentialTestCase {
  pub name: String,
  pub description: String,
  pub stress_type: String,
  pub analytical_caution: String,
  pub slowness: f64,
  pub a: f64,
  pub b: f64,
  pub input_grid_domain: (f64, f64),
  pub input_grid_n_points: usize,
  
  
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

  fn slowness(&self) -> f64 {
    self.slowness
  }

  fn input_grid_domain(&self) -> (f64, f64) {
    self.input_grid_domain
  }

  fn input_grid_n_points(&self) -> usize {
    self.input_grid_n_points
  }

}
