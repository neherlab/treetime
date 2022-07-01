use eyre::Report;
use ndarray::{Array1, Array2, Axis, Ix1, Ix2};
use num_traits::Bounded;

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
  Ok(arr.fold_axis(axis, f32::MAX, |a, b| a.min(*b)))
}

/// Calculates max sum over given axis
#[inline]
pub fn max_axis(arr: &Array2<f32>, axis: Axis) -> Result<Array1<f32>, Report> {
  Ok(arr.fold_axis(axis, f32::MIN, |a, b| a.max(*b)))
}

/// Finds index of min value over given axis
#[inline]
pub fn argmin_axis<T: Copy + Ord + PartialOrd + Bounded>(arr: &Array1<T>, axis: Axis) -> usize {
  let mut i_curr = 0;
  let res = arr
    .fold_axis(axis, (0_usize, T::max_value()), |(i_min, x_min), x| {
      let res = if x < x_min { (i_curr, *x) } else { (*i_min, *x_min) };
      i_curr += 1;
      res
    })
    .into_raw_vec();
  res[0].0
}

/// Finds index of max value over given axis
#[inline]
pub fn argmax_axis<T: Copy + Ord + PartialOrd + Bounded>(arr: &Array1<T>, axis: Axis) -> usize {
  let mut i_curr = 0;
  let res = arr
    .fold_axis(axis, (0_usize, T::min_value()), |(i_max, x_max), x| {
      let res = if x > x_max { (i_curr, *x) } else { (*i_max, *x_max) };
      i_curr += 1;
      res
    })
    .into_raw_vec();
  res[0].0
}

/// Calculates cumulative sum over given axis
#[inline]
pub fn cumsum_axis(profile: &Array2<f32>, axis: Axis) -> Result<Array1<f32>, Report> {
  let mut result = profile.to_owned();
  result.accumulate_axis_inplace(axis, |&prev, curr| *curr += prev);
  let result = result.into_dimensionality::<Ix1>()?;
  Ok(result)
}
