use eyre::Report;
use ndarray::{Array1, Array2, Axis, Ix1, Ix2};

// Calculates outer product of 2 vectors
pub fn outer(a: &Array1<f32>, b: &Array1<f32>) -> Result<Array2<f32>, Report> {
  let a = a.clone().into_dimensionality::<Ix2>()?.into_shape((a.len(), 1_usize))?;

  let b = b
    .t()
    .to_owned()
    .into_dimensionality::<Ix2>()?
    .into_shape((a.len(), 1_usize))?; // reshape((1, b.len()));

  Ok(a.dot(&b))
}

/// Calculates min sum over given axis
#[inline]
pub fn min_axis(arr: &Array2<f32>, axis: Axis) -> Result<Array1<f32>, Report> {
  Ok(arr.fold_axis(axis, -f32::MIN, |a, b| a.min(*b)))
}

/// Calculates max sum over given axis
#[inline]
pub fn max_axis(arr: &Array2<f32>, axis: Axis) -> Result<Array1<f32>, Report> {
  Ok(arr.fold_axis(axis, -f32::MIN, |a, b| a.max(*b)))
}

/// Calculates cumulative sum over given axis
#[inline]
pub fn cumsum_axis(profile: &Array2<f32>, axis: Axis) -> Result<Array1<f32>, Report> {
  let mut result = profile.to_owned();
  result.accumulate_axis_inplace(axis, |&prev, curr| *curr += prev);
  result.into_dimensionality::<Ix1>()
}
