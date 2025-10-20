use crate::algos::algos::Algo;
use eyre::Report;
use itertools::izip;
use ndarray::Array1;
use ndarray_interp::interp1d::Interp1DBuilder;

/// Convolution using Riemann sum integration
pub fn convolve_riemann(
  input_grid: &Array1<f64>,
  f_values: &Array1<f64>,
  g_values: &Array1<f64>,
  output_grid: &Array1<f64>,
) -> Result<Array1<f64>, Report> {
  let g_interp = Interp1DBuilder::new(g_values.view()).x(input_grid.view()).build()?;

  let ds = input_grid[1] - input_grid[0];

  let mut result = Array1::zeros(output_grid.len());
  for (output_idx, &x_eval) in output_grid.iter().enumerate() {
    let mut sum = 0.0;
    for (&x_input, &f_at_x_input) in izip!(input_grid, f_values) {
      let g_at_shifted = g_interp.interp_scalar(x_eval - x_input).unwrap_or(0.0);
      sum += f_at_x_input * g_at_shifted;
    }
    result[output_idx] = sum * ds;
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
