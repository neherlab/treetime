use crate::make_report;
use eyre::Report;
use ndarray::{Ix0, LinalgScalar};
use ndarray_einsum_beta::{einsum as einsum_base, ArrayLike};

#[inline]
pub fn einsum_1d<A: LinalgScalar>(input_string: &str, operands: &[&dyn ArrayLike<A>]) -> Result<A, Report> {
  let result = einsum_base(input_string, operands)
    .map_err(|err| make_report!("einsum: {err}"))?
    .into_dimensionality::<Ix0>()?
    .into_scalar();
  Ok(result)
}

#[allow(non_snake_case)]
#[cfg(test)]
mod tests {
  use super::*;
  use approx::assert_ulps_eq;
  use eyre::Report;
  use ndarray::array;
  use rstest::rstest;

  #[rstest]
  fn computes_einsum_simple() -> Result<(), Report> {
    let W = array![
      [0., 1., 1., 1., 1.],
      [1., 0., 1., 1., 1.],
      [1., 1., 0., 1., 1.],
      [1., 1., 1., 0., 1.],
      [1., 1., 1., 1., 0.]
    ];

    let pi = array![0.2, 0.2, 0.2, 0.2, 0.2];

    assert_ulps_eq!(einsum_1d::<f64>("i,ij,j", &[&pi, &W, &pi])?, 0.8000000000000005);

    Ok(())
  }
}
