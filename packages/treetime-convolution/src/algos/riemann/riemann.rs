use crate::algos::algo_trait::Algo;
use eyre::Report;
use ndarray::Array1;
use ndarray_interp::interp1d::Interp1DBuilder;

/// Convolution using Riemann sum integration
pub fn convolve_riemann(
  input_grid: &Array1<f64>,
  f_values: &Array1<f64>,
  g_values: &Array1<f64>,
  output_grid: &Array1<f64>,
) -> Result<Array1<f64>, Report> {
  let f_interp = Interp1DBuilder::new(f_values.clone()).x(input_grid.clone()).build()?;
  let g_interp = Interp1DBuilder::new(g_values.clone()).x(input_grid.clone()).build()?;

  let ds = input_grid[1] - input_grid[0];

  let mut result = Array1::zeros(output_grid.len());
  for (i, &x_eval) in output_grid.iter().enumerate() {
    let mut sum = 0.0;
    for &s in input_grid {
      let f_at_s = f_interp.interp_scalar(s).unwrap_or(0.0);
      let g_at_x_minus_s = g_interp.interp_scalar(x_eval - s).unwrap_or(0.0);
      sum += f_at_s * g_at_x_minus_s;
    }
    result[i] = sum * ds;
  }

  Ok(result)
}

pub struct RiemannAlgo;

impl Algo for RiemannAlgo {
  fn name(&self) -> &'static str {
    "riemann"
  }

  fn convolve(
    &self,
    input_grid: &Array1<f64>,
    f_values: &Array1<f64>,
    g_values: &Array1<f64>,
    output_grid: &Array1<f64>,
  ) -> Result<Array1<f64>, Report> {
    convolve_riemann(input_grid, f_values, g_values, output_grid)
  }
}
