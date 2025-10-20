use crate::distribution::reference::convolution_test::algorithms::ConvolutionAlgorithm;
use crate::distribution::reference::convolution_test::framework::results::TestResult;
use crate::distribution::reference::convolution_test::framework::test_case::TestCase;
use crate::distribution::reference::convolution_test::metrics::metrics::ConvolutionMetrics;
use crate::distribution::reference::convolution_test::traits::{ConvAlgo, ConvInput};
use eyre::Report;
use ndarray::Array1;
use std::time::Instant;

pub trait ConvolutionTestRunner<T: TestCase>: Send + Sync {
  fn run_test(&self, test_case: &T, algorithm: ConvolutionAlgorithm) -> Result<TestResult<T>, Report>;

  fn test_cases(&self) -> &[T];

  fn function_type(&self) -> &str;
}

pub struct TraitBasedTestRunner<I: ConvInput> {
  input: I,
}

impl<I: ConvInput> TraitBasedTestRunner<I> {
  pub fn new(input: I) -> Self {
    Self { input }
  }
}

impl<I: ConvInput> ConvolutionTestRunner<I::TestCase> for TraitBasedTestRunner<I> {
  fn run_test(
    &self,
    test_case: &I::TestCase,
    algorithm: ConvolutionAlgorithm,
  ) -> Result<TestResult<I::TestCase>, Report> {
    let start_time = Instant::now();

    let f = self.input.create_f(test_case)?;
    let g = self.input.create_g(test_case)?;

    let (eval_min, eval_max) = self.input.eval_domain(test_case);
    let n_eval_points = ((eval_max - eval_min) / test_case.dx() + 1.0).round() as usize;
    let eval_grid = Array1::from_iter((0..n_eval_points).map(|i| eval_min + i as f64 * test_case.dx()));

    let algo: Box<dyn ConvAlgo> = match algorithm {
      ConvolutionAlgorithm::Riemann => {
        Box::new(crate::distribution::reference::convolution_test::algo_impls::RiemannAlgo)
      },
      ConvolutionAlgorithm::Ndarray => {
        Box::new(crate::distribution::reference::convolution_test::algo_impls::NdarrayAlgo)
      },
    };

    let actual_result = algo.convolve(&f, &g, &eval_grid)?;
    let expected_result = self.input.analytical_convolution(test_case, &eval_grid)?;

    let execution_time = start_time.elapsed().as_secs_f64() * 1000.0;

    let metrics = ConvolutionMetrics::new(
      actual_result.x(),
      actual_result.y(),
      expected_result.y(),
      execution_time,
    )?;

    let evaluation_grid = actual_result.x().to_owned();
    let actual_values = actual_result.y().to_owned();
    let expected_values = expected_result.y().to_owned();
    let f_x_values = f.x().to_owned();
    let f_y_values = f.y().to_owned();
    let g_x_values = g.x().to_owned();
    let g_y_values = g.y().to_owned();

    Ok(TestResult {
      algorithm,
      test_case: test_case.clone(),
      execution_time_ms: execution_time,
      f_x_values,
      f_y_values,
      g_x_values,
      g_y_values,
      evaluation_grid,
      actual_values,
      expected_values,
      metrics,
    })
  }

  fn test_cases(&self) -> &[I::TestCase] {
    self.input.test_cases()
  }

  fn function_type(&self) -> &'static str {
    self.input.function_type()
  }
}
